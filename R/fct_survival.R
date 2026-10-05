#' Shared survival / prognosis layer
#'
#' Omics-agnostic survival analysis used by every pipeline that has clinical
#' follow-up: single-cell (cell-type composition vs outcome), WES (mutation vs
#' outcome), bulk (molecular subtype vs outcome) and multi-omics integration
#' (integrative cluster vs outcome). The normalised clinical table lives on the
#' shared hub as `rv$clinical`, so it is loaded once and reused everywhere.
#'
#' Everything here is a pure function of data frames -- no Seurat, no Shiny --
#' so it is directly testable, and the maths uses `survival`, which ships with R
#' and therefore installs on any machine (no compiled GitHub packages, no
#' `survminer`; the Kaplan-Meier curve is drawn with ggplot2 here).
#'
#' The order the callers must keep (AGENTS.md, section 7): define the analysis
#' set first -- one row per patient, missing times/events/values removed --
#' and only then derive a cut-off with [split_numeric()].
#'
#' @name fct_survival
#' @keywords internal
NULL

# ---- clinical table ---------------------------------------------------------

#' Guess the field separator of a delimited text file
#'
#' Looks at the extension first (.csv / .tsv / .txt), then counts tabs,
#' commas and semicolons on the first line.
#' @param path File path. @param name Original file name (uploads have a temp
#'   path without the user's extension).
#' @return One of "\\t", ",", ";".
#' @keywords internal
guess_sep <- function(path, name = basename(path)) {
  ext <- tolower(tools::file_ext(sub("\\.gz$", "", name)))
  if (ext == "csv") return(",")
  if (ext %in% c("tsv", "tab")) return("\t")
  first <- tryCatch(readLines(path, n = 1, warn = FALSE), error = function(e) "")
  n <- c("\t" = lengths(regmatches(first, gregexpr("\t", first))),
         ","  = lengths(regmatches(first, gregexpr(",", first))),
         ";"  = lengths(regmatches(first, gregexpr(";", first))))
  if (all(n == 0)) return("\t")
  names(n)[which.max(n)]
}

#' Read a clinical table (csv/tsv) with the sample ids as an ordinary column
#'
#' Columns are read as text and then converted one by one, except columns
#' that look like identifiers with leading zeros ("001"), which stay text so
#' they still match the sample names elsewhere.
#' @param path File path.
#' @param sep Field separator, or "auto" to detect it with [guess_sep()].
#' @param name Original file name, used by the separator guess.
#' @keywords internal
read_clinical_table <- function(path, sep = "auto", name = basename(path)) {
  if (identical(sep, "auto")) sep <- guess_sep(path, name)
  df <- utils::read.delim(path, sep = sep, check.names = FALSE,
                          colClasses = "character", na.strings = c("", "NA"))
  for (cl in names(df)) {
    v <- df[[cl]]
    if (any(grepl("^0[0-9]", v))) next
    df[[cl]] <- utils::type.convert(v, as.is = TRUE)
  }
  df
}

#' Normalise a clinical table to the columns the survival code expects
#'
#' Produces `.id` / `.time` / `.event` alongside the untouched original columns,
#' so downstream code never has to care what the user's columns were called.
#' Rows with a missing/negative time or an unusable event value are dropped --
#' `survival` would otherwise fail with a much less obvious message. Survival
#' analysis needs one row per patient, so duplicated ids are an error unless
#' the caller explicitly keeps the first occurrence.
#'
#' @param df Raw clinical data.frame.
#' @param id_col,time_col,event_col Column names chosen by the user.
#' @param event_positive Value of `event_col` that means "event occurred"
#'   (death / progression). Numeric 1 and the usual text codings are detected
#'   automatically when this is NULL.
#' @param time_unit One of "days", "months", "years": the unit of `time_col`.
#'   Times are converted to months, which is what the plots label.
#' @param dedup What to do with repeated ids: "error" (default), "first" (keep
#'   the first row of each id) or "keep" (no check; only for callers that
#'   aggregate afterwards).
#' @return data.frame with `.id`, `.time` (months), `.event` (0/1) plus the
#'   original columns. Attributes: "dropped" (rows removed for unusable
#'   time/event), "duplicates" (number of rows removed as repeated ids).
#' @keywords internal
normalise_clinical <- function(df, id_col, time_col, event_col,
                               event_positive = NULL,
                               time_unit = c("months", "days", "years"),
                               dedup = c("error", "first", "keep")) {
  time_unit <- match.arg(time_unit)
  dedup <- match.arg(dedup)
  stopifnot(is.data.frame(df))
  for (cl in c(id_col, time_col, event_col)) {
    if (!cl %in% names(df)) stop("Column not found in the clinical table: ", cl)
  }

  time <- suppressWarnings(as.numeric(df[[time_col]]))
  time <- switch(time_unit, months = time, days = time / 30.4375, years = time * 12)

  event <- encode_event(df[[event_col]], event_positive)

  out <- df
  out$.id    <- trimws(as.character(df[[id_col]]))
  out$.time  <- time
  out$.event <- event

  keep <- !is.na(out$.id) & nzchar(out$.id) & !is.na(out$.time) &
    is.finite(out$.time) & out$.time >= 0 & !is.na(out$.event)
  dropped <- sum(!keep)
  out <- out[keep, , drop = FALSE]

  dup <- duplicated(out$.id)
  n_dup <- 0L
  if (any(dup) && dedup == "error") {
    ids <- unique(out$.id[dup])
    stop(sprintf(paste("%d patient id(s) occur more than once (e.g. %s).",
                       "Survival analysis needs one row per patient."),
                 length(ids), paste(utils::head(ids, 3), collapse = ", ")))
  }
  if (any(dup) && dedup == "first") {
    n_dup <- sum(dup)
    out <- out[!dup, , drop = FALSE]
  }
  attr(out, "dropped") <- dropped
  attr(out, "duplicates") <- n_dup
  out
}

#' Coerce an event column to 0/1
#'
#' Accepts 0/1, TRUE/FALSE, 1/2 (the `survival` convention), and the common text
#' codings. `positive` forces a specific level to mean "event".
#'
#' @param x The raw event column. @param positive Level meaning "event", or NULL.
#' @return Integer vector of 0/1 (NA where undecidable).
#' @keywords internal
encode_event <- function(x, positive = NULL) {
  if (!is.null(positive) && nzchar(as.character(positive))) {
    return(as.integer(as.character(x) == as.character(positive)))
  }
  if (is.logical(x)) return(as.integer(x))
  if (is.numeric(x)) {
    u <- sort(unique(stats::na.omit(x)))
    # survival's other convention is 1 = censored, 2 = event
    if (length(u) == 2 && all(u == c(1, 2))) return(as.integer(x == 2))
    return(as.integer(x > 0))
  }
  s <- tolower(trimws(as.character(x)))
  dead  <- c("1", "true", "yes", "y", "dead", "deceased", "death", "event",
             "progressed", "progression", "recurrence", "relapse")
  alive <- c("0", "false", "no", "n", "alive", "living", "censored", "censor",
             "no event", "disease free", "disease-free")
  out <- rep(NA_integer_, length(s))
  out[s %in% dead]  <- 1L
  out[s %in% alive] <- 0L
  out
}

# ---- grouping ---------------------------------------------------------------

#' Split a numeric variable into survival groups
#'
#' Call it on the final analysis set only (one row per patient, missing values
#' removed): the cut-off is a property of that set.
#'
#' * "median": above vs at-or-below the median.
#' * "tertile": top vs bottom third, middle third dropped. Fails when ties make
#'   the two tertile cut-offs equal (common with zero-inflated fractions).
#' * "optimal": the cut-point maximising the log-rank statistic within the
#'   central `[min_frac, 1 - min_frac]` range. Selecting the cut-point inflates
#'   the plain log-rank p-value (about 30% false positives at a nominal 5% in
#'   simulation), so the result carries the Lausen-Schumacher adjusted p-value
#'   in `attr(, "p_adjusted")` -- from maxstat when installed, otherwise the
#'   closed-form Lausen & Schumacher (1992) approximation. Report that one.
#'
#' @param x Numeric vector (e.g. a cell-type proportion or an expression score).
#' @param method "median", "tertile" or "optimal".
#' @param time,event Needed only by "optimal".
#' @param min_frac Smallest share of samples allowed in a group by "optimal".
#' @return A factor (levels Low, High) the same length as `x`; NA marks samples
#'   not assigned. Attributes: "cutpoint" (the threshold used, or for tertiles
#'   both thresholds), "method" (the method actually applied), and for
#'   "optimal" also "p_adjusted" and "chisq".
#' @keywords internal
split_numeric <- function(x, method = c("median", "tertile", "optimal"),
                          time = NULL, event = NULL, min_frac = 0.2) {
  method <- match.arg(method)
  x <- as.numeric(x)

  if (method == "median") {
    cut <- stats::median(x, na.rm = TRUE)
    out <- factor(ifelse(x > cut, "High", "Low"), levels = c("Low", "High"))
    attr(out, "cutpoint") <- cut
    attr(out, "method") <- "median"
    return(out)
  }

  if (method == "tertile") {
    q <- stats::quantile(x, probs = c(1 / 3, 2 / 3), na.rm = TRUE)
    if (q[[1]] == q[[2]]) {
      stop("Too many tied values for tertiles (both cut-offs are ", signif(q[[1]], 3),
           "). Use the median split or a continuous Cox model instead.")
    }
    g <- rep(NA_character_, length(x))
    g[x <= q[[1]]] <- "Low"
    g[x >= q[[2]]] <- "High"
    out <- factor(g, levels = c("Low", "High"))
    attr(out, "cutpoint") <- unname(q)
    attr(out, "method") <- "tertile"
    return(out)
  }

  # "optimal": scan candidate cut-points, keep the one with the largest
  # log-rank chi-square, bounded by min_frac so it cannot "find" a split of
  # three samples.
  if (is.null(time) || is.null(event)) {
    stop("The optimal cutpoint needs the survival time and event columns.")
  }
  ok <- !is.na(x) & !is.na(time) & !is.na(event)
  fallback <- function() {
    out <- split_numeric(x, "median")
    attr(out, "method") <- "median (too few distinct values for an optimal cut-point)"
    out
  }
  cand <- sort(unique(x[ok]))
  if (length(cand) < 3) return(fallback())
  lo <- stats::quantile(x[ok], min_frac, na.rm = TRUE)
  hi <- stats::quantile(x[ok], 1 - min_frac, na.rm = TRUE)
  cand <- cand[cand >= lo & cand < hi]
  if (!length(cand)) return(fallback())

  best <- NULL
  best_stat <- -Inf
  for (cut in cand) {
    g <- factor(ifelse(x > cut, "High", "Low"), levels = c("Low", "High"))
    st <- tryCatch(
      survival::survdiff(survival::Surv(time[ok], event[ok]) ~ g[ok])$chisq,
      error = function(e) NA_real_)
    if (!is.na(st) && st > best_stat) {
      best_stat <- st
      best <- cut
    }
  }
  if (is.null(best)) return(fallback())
  out <- factor(ifelse(x > best, "High", "Low"), levels = c("Low", "High"))
  attr(out, "cutpoint") <- best
  attr(out, "method") <- "optimal"
  attr(out, "chisq") <- best_stat
  attr(out, "p_adjusted") <- maxstat_p(x[ok], time[ok], event[ok], best_stat, min_frac)
  out
}

#' Selection-adjusted p-value for an optimal log-rank cut-point
#'
#' Uses `maxstat::maxstat.test(pmethod = "Lau94")` when maxstat is installed;
#' otherwise the Lausen & Schumacher (1992) approximation for the maximum of
#' standardised log-rank statistics over the candidate range:
#' `phi(b) * (b - 1/b) * log(e2 (1 - e1) / ((1 - e2) e1)) + 4 phi(b) / b`,
#' with `b = sqrt(max chi-square)`, `e1 = min_frac`, `e2 = 1 - min_frac`.
#' @param x,time,event Complete-case vectors. @param chisq Maximum chi-square.
#' @param min_frac Lower bound of the candidate range.
#' @keywords internal
maxstat_p <- function(x, time, event, chisq, min_frac = 0.2) {
  # an adjusted p is never below the unadjusted one, and never above 1 (both
  # approximations can overshoot when the statistic is small)
  p_raw <- stats::pchisq(chisq, 1, lower.tail = FALSE)
  bound <- function(p) min(1, max(p, p_raw))
  if (has_pkg("maxstat")) {
    d <- data.frame(x = x, time = time, event = event)
    p <- tryCatch(
      maxstat::maxstat.test(survival::Surv(time, event) ~ x, data = d,
                            smethod = "LogRank", pmethod = "Lau94",
                            minprop = min_frac, maxprop = 1 - min_frac)$p.value,
      error = function(e) NA_real_)
    if (is.finite(p)) return(bound(unname(p)))
  }
  b <- sqrt(chisq)
  if (!is.finite(b) || b <= 0) return(1)
  e1 <- min_frac
  e2 <- 1 - min_frac
  phi <- stats::dnorm(b)
  bound(phi * (b - 1 / b) * log(e2 * (1 - e1) / ((1 - e2) * e1)) + 4 * phi / b)
}

#' Per-sample composition (fraction of each group within each sample)
#'
#' The bridge from a per-observation table to the per-sample table survival
#' analysis needs: for single-cell that is cells -> cell-type fractions per
#' sample; the same function serves any per-observation categorical breakdown.
#'
#' @param meta data.frame with one row per observation (e.g. per cell).
#' @param sample_col Column holding the sample/patient id.
#' @param group_col Column holding the category (cluster, cell type, ...).
#' @return data.frame: `.id`, one column per category (fractions summing to 1),
#'   and `n_obs`.
#' @keywords internal
composition_by_sample <- function(meta, sample_col, group_col) {
  if (!all(c(sample_col, group_col) %in% names(meta))) {
    stop("Sample or group column not found in the object metadata.")
  }
  s <- as.character(meta[[sample_col]])
  g <- as.character(meta[[group_col]])
  ok <- !is.na(s) & !is.na(g)
  tab <- table(s[ok], g[ok])
  frac <- as.data.frame.matrix(tab / rowSums(tab))
  frac$.id <- rownames(frac)
  frac$n_obs <- as.integer(rowSums(tab))
  rownames(frac) <- NULL
  frac[, c(".id", setdiff(names(frac), c(".id", "n_obs")), "n_obs"), drop = FALSE]
}

#' Per-patient composition for outcome analysis
#'
#' The patient-level version of [composition_by_sample()]: optionally keep
#' only some cells first (e.g. tumour tissue at baseline), optionally restrict
#' the denominator to some groups (fraction *of T cells* rather than of all
#' cells), and drop patients with too few cells for a stable fraction.
#'
#' @param meta Cell metadata. @param patient_col Patient id column.
#' @param group_col Category column.
#' @param keep Optional logical vector over rows of `meta`: cells to use.
#' @param denom_levels Optional categories forming the denominator.
#' @param min_cells Patients with fewer cells (after `keep`/`denom_levels`)
#'   are excluded.
#' @return Like [composition_by_sample()] (`.id` = patient); attribute
#'   "excluded" lists the patients dropped by `min_cells`.
#' @keywords internal
composition_by_patient <- function(meta, patient_col, group_col, keep = NULL,
                                   denom_levels = NULL, min_cells = 50) {
  if (!is.null(keep)) meta <- meta[!is.na(keep) & keep, , drop = FALSE]
  if (length(denom_levels)) {
    meta <- meta[as.character(meta[[group_col]]) %in% denom_levels, , drop = FALSE]
  }
  comp <- composition_by_sample(meta, patient_col, group_col)
  low <- comp$n_obs < min_cells
  out <- comp[!low, , drop = FALSE]
  attr(out, "excluded") <- comp$.id[low]
  out
}

# ---- models -----------------------------------------------------------------

#' Kaplan-Meier fit for a normalised clinical table
#' @param df Output of [normalise_clinical()], optionally with a `.group` column.
#' @keywords internal
km_fit <- function(df) {
  if (".group" %in% names(df) && length(unique(stats::na.omit(df$.group))) > 1) {
    survival::survfit(survival::Surv(.time, .event) ~ .group, data = df)
  } else {
    survival::survfit(survival::Surv(.time, .event) ~ 1, data = df)
  }
}

#' Log-rank test across `.group`
#' @param df Normalised clinical table with a `.group` column.
#' @return list(chisq, df, p) or NULL when there is nothing to compare.
#' @keywords internal
logrank_test <- function(df) {
  if (!".group" %in% names(df)) return(NULL)
  d <- df[!is.na(df$.group), , drop = FALSE]
  if (length(unique(d$.group)) < 2) return(NULL)
  sd <- tryCatch(survival::survdiff(survival::Surv(.time, .event) ~ .group, data = d),
                 error = function(e) NULL)
  if (is.null(sd)) return(NULL)
  k <- length(sd$n) - 1
  list(chisq = unname(sd$chisq), df = k,
       p = stats::pchisq(sd$chisq, df = k, lower.tail = FALSE))
}

#' Median survival per group (months)
#' @param fit A [km_fit()] result.
#' @return data.frame: group, n, events, median, lower, upper. `median` is NA
#'   when it is not reached; show it with [format_median()].
#' @keywords internal
km_medians <- function(fit) {
  s <- summary(fit)$table
  if (is.null(dim(s))) s <- t(as.matrix(s))
  gp <- rownames(s)
  if (is.null(gp)) gp <- "All"
  pick <- function(nm) if (nm %in% colnames(s)) unname(s[, nm]) else rep(NA_real_, nrow(s))
  data.frame(
    group  = sub("^\\.group=", "", gp),
    n      = pick("records"),
    events = pick("events"),
    median = pick("median"),
    lower  = pick("0.95LCL"),
    upper  = pick("0.95UCL"),
    stringsAsFactors = FALSE
  )
}

#' Median survival as text: "NR" (not reached) instead of NA
#' @param x Numeric medians. @param digits Rounding.
#' @keywords internal
format_median <- function(x, digits = 1) {
  ifelse(is.na(x), "NR", format(round(x, digits), nsmall = digits))
}

#' Hazard ratio of the second `.group` level against the first
#'
#' The effect size that belongs next to a two-arm Kaplan-Meier plot.
#' Warnings from `coxph` (an infinite coefficient when one arm has no events)
#' are returned in `note` instead of being lost on the console.
#' @param df Normalised clinical table with a two-level `.group`.
#' @return list(hr, lower, upper, p (Wald), n, events, label, note), or NULL.
#' @keywords internal
cox_hr <- function(df) {
  if (!".group" %in% names(df)) return(NULL)
  d <- df[!is.na(df$.group), , drop = FALSE]
  d$.group <- droplevels(factor(d$.group))
  if (nlevels(d$.group) != 2) return(NULL)
  notes <- character(0)
  fit <- withCallingHandlers(
    tryCatch(survival::coxph(survival::Surv(.time, .event) ~ .group, data = d),
             error = function(e) NULL),
    warning = function(w) {
      notes <<- c(notes, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  if (is.null(fit)) return(NULL)
  sm <- summary(fit)
  lv <- levels(d$.group)
  list(hr = unname(sm$conf.int[1, "exp(coef)"]),
       lower = unname(sm$conf.int[1, "lower .95"]),
       upper = unname(sm$conf.int[1, "upper .95"]),
       p = unname(sm$coefficients[1, ncol(sm$coefficients)]),
       n = nrow(d), events = sum(d$.event),
       label = sprintf("%s vs %s", lv[2], lv[1]),
       note = unique(notes))
}

#' Hazard ratio per standard deviation of a continuous variable
#'
#' The primary analysis for a continuous biomarker: no cut-off to choose, so
#' no selection bias; the Kaplan-Meier split is then only an illustration.
#' @param df Normalised clinical table. @param var Numeric column.
#' @return list(hr, lower, upper, p, n, events, sd, label, note), or NULL.
#' @keywords internal
cox_continuous <- function(df, var) {
  d <- df[!is.na(df[[var]]), c(".time", ".event", var), drop = FALSE]
  s <- stats::sd(d[[var]])
  if (nrow(d) < 5 || !is.finite(s) || s == 0) return(NULL)
  d$z <- (d[[var]] - mean(d[[var]])) / s
  notes <- character(0)
  fit <- withCallingHandlers(
    tryCatch(survival::coxph(survival::Surv(.time, .event) ~ z, data = d),
             error = function(e) NULL),
    warning = function(w) {
      notes <<- c(notes, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  if (is.null(fit)) return(NULL)
  sm <- summary(fit)
  list(hr = unname(sm$conf.int[1, "exp(coef)"]),
       lower = unname(sm$conf.int[1, "lower .95"]),
       upper = unname(sm$conf.int[1, "upper .95"]),
       p = unname(sm$coefficients[1, ncol(sm$coefficients)]),
       n = nrow(d), events = sum(d$.event), sd = s,
       label = sprintf("per SD of %s", var), note = unique(notes))
}

#' Proportional-hazards check for the two-arm Cox model
#' @param df Normalised clinical table with `.group`.
#' @return The global `cox.zph` p-value, or NA when it cannot be computed.
#' @keywords internal
cox_ph_p <- function(df) {
  d <- df[!is.na(df$.group), , drop = FALSE]
  if (length(unique(d$.group)) < 2) return(NA_real_)
  fit <- tryCatch(suppressWarnings(
    survival::coxph(survival::Surv(.time, .event) ~ .group, data = d)),
    error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  z <- tryCatch(survival::cox.zph(fit), error = function(e) NULL)
  if (is.null(z)) return(NA_real_)
  unname(z$table[nrow(z$table), "p"])
}

#' Small-sample warnings for a survival comparison
#'
#' Thresholds follow common practice (AGENTS.md section 7): fewer than 10
#' patients in a group, fewer than 10 events overall, or fewer than 5 events in
#' a group make estimates unstable and the analysis exploratory.
#' @param df Normalised clinical table, optionally with `.group`.
#' @param lang "en" or "zh".
#' @return Character vector of warnings, possibly empty.
#' @keywords internal
survival_warnings <- function(df, lang = c("en", "zh")) {
  lang <- match.arg(lang)
  msg <- function(en, zh) if (lang == "zh") zh else en
  out <- character(0)
  d <- if (".group" %in% names(df)) df[!is.na(df$.group), , drop = FALSE] else df
  if (sum(d$.event) < 10) {
    out <- c(out, sprintf(msg("Only %d events in total.", "总事件数只有 %d 个。"), sum(d$.event)))
  }
  if (".group" %in% names(d) && nrow(d)) {
    g <- droplevels(factor(d$.group))
    n <- table(g)
    e <- tapply(d$.event, g, sum)
    small <- names(n)[n < 10]
    if (length(small)) {
      out <- c(out, sprintf(msg("Fewer than 10 patients in: %s.", "少于 10 位患者的组：%s。"),
                            paste(small, collapse = ", ")))
    }
    few <- names(e)[e < 5]
    if (length(few)) {
      out <- c(out, sprintf(msg("Fewer than 5 events in: %s.", "少于 5 个事件的组：%s。"),
                            paste(few, collapse = ", ")))
    }
  }
  out
}

#' Univariable Cox regression, one model per variable
#'
#' Univariable (not multivariable) on purpose: it is the screen users want at
#' this point, and it stays interpretable when the cohort is small. Each
#' variable gets one overall likelihood-ratio p-value (one test even for a
#' multi-level factor) and a BH q-value across the variables screened; rows
#' stay grouped by variable, variables ordered by that overall p.
#'
#' @param df Normalised clinical table.
#' @param vars Column names to test.
#' @return data.frame: variable, level, n, events, HR, lower, upper, p (Wald,
#'   per level), p_overall, q, note.
#' @keywords internal
cox_univariable <- function(df, vars) {
  rows <- lapply(vars, function(v) {
    d <- df[, c(".time", ".event", v), drop = FALSE]
    names(d)[3] <- "x"
    if (is.character(d$x) || is.logical(d$x)) d$x <- factor(d$x)
    d <- d[!is.na(d$x), , drop = FALSE]
    if (nrow(d) < 5 || sum(d$.event) < 2) return(NULL)
    if (is.factor(d$x)) {
      d$x <- droplevels(d$x)
      if (nlevels(d$x) < 2) return(NULL)
    }
    notes <- character(0)
    fit <- withCallingHandlers(
      tryCatch(survival::coxph(survival::Surv(.time, .event) ~ x, data = d),
               error = function(e) NULL),
      warning = function(w) {
        notes <<- c(notes, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    if (is.null(fit)) return(NULL)
    sm <- summary(fit)
    cf <- sm$coefficients
    ci <- sm$conf.int
    data.frame(
      variable  = v,
      level     = sub("^x", "", rownames(cf)),
      n         = nrow(d),
      events    = sum(d$.event),
      HR        = unname(cf[, "exp(coef)"]),
      lower     = unname(ci[, "lower .95"]),
      upper     = unname(ci[, "upper .95"]),
      p         = unname(cf[, ncol(cf)]),
      p_overall = unname(sm$logtest[["pvalue"]]),
      note      = if (length(notes)) paste(unique(notes), collapse = "; ") else "",
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(out)) {
    return(data.frame(variable = character(0), level = character(0), n = integer(0),
                      events = integer(0), HR = numeric(0), lower = numeric(0),
                      upper = numeric(0), p = numeric(0), p_overall = numeric(0),
                      q = numeric(0), note = character(0), stringsAsFactors = FALSE))
  }
  per_var <- out[!duplicated(out$variable), c("variable", "p_overall")]
  per_var$q <- stats::p.adjust(per_var$p_overall, method = "BH")
  out$q <- per_var$q[match(out$variable, per_var$variable)]
  ord <- order(match(out$variable, per_var$variable[order(per_var$p_overall)]), out$p)
  out <- out[ord, c("variable", "level", "n", "events", "HR", "lower", "upper",
                    "p", "p_overall", "q", "note"), drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Multivariable Cox model: the stratifier adjusted for covariates
#'
#' Complete cases only. Returns the per-level table plus the events-per-
#' variable ratio (`attr(, "epv")`): below 10 the estimates are unstable and
#' the caller must say so.
#' @param df Normalised clinical table. @param vars Model terms (column names).
#' @return data.frame: variable, level, HR, lower, upper, p; attributes "n",
#'   "events", "epv", "note". NULL when the model cannot be fitted.
#' @keywords internal
cox_multivariable <- function(df, vars) {
  vars <- unique(vars[vars %in% names(df)])
  if (!length(vars)) return(NULL)
  d <- df[stats::complete.cases(df[, c(".time", ".event", vars), drop = FALSE]),
          c(".time", ".event", vars), drop = FALSE]
  for (v in vars) {
    if (is.character(d[[v]]) || is.logical(d[[v]])) d[[v]] <- factor(d[[v]])
    if (is.factor(d[[v]])) d[[v]] <- droplevels(d[[v]])
  }
  vars <- vars[vapply(vars, function(v) !is.factor(d[[v]]) || nlevels(d[[v]]) > 1,
                      logical(1))]
  if (!length(vars) || nrow(d) < 5) return(NULL)
  notes <- character(0)
  f <- stats::reformulate(sprintf("`%s`", vars), "survival::Surv(.time, .event)")
  fit <- withCallingHandlers(
    tryCatch(survival::coxph(f, data = d), error = function(e) NULL),
    warning = function(w) {
      notes <<- c(notes, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  if (is.null(fit)) return(NULL)
  sm <- summary(fit)
  cf <- sm$coefficients
  ci <- sm$conf.int
  term <- gsub("`", "", rownames(cf))
  var_of <- vapply(term, function(t) {
    hit <- vars[startsWith(t, vars)]
    if (length(hit)) hit[which.max(nchar(hit))] else t
  }, character(1))
  out <- data.frame(variable = unname(var_of),
                    level = unname(substring(term, nchar(var_of) + 1)),
                    HR = unname(cf[, "exp(coef)"]),
                    lower = unname(ci[, "lower .95"]),
                    upper = unname(ci[, "upper .95"]),
                    p = unname(cf[, ncol(cf)]),
                    stringsAsFactors = FALSE)
  attr(out, "n") <- nrow(d)
  attr(out, "events") <- sum(d$.event)
  attr(out, "epv") <- sum(d$.event) / nrow(cf)
  attr(out, "note") <- unique(notes)
  out
}

# ---- plot -------------------------------------------------------------------

#' Tidy a survfit into a step-plottable data.frame (with the t=0 anchor)
#'
#' Groups keep the order of the fit's strata (Low before High), and the 95%
#' confidence limits are carried along for the band.
#' @param fit A [km_fit()] result.
#' @keywords internal
km_tidy <- function(fit) {
  lv <- if (is.null(fit$strata)) "All" else sub("^\\.group=", "", names(fit$strata))
  strata <- if (is.null(fit$strata)) rep("All", length(fit$time))
            else rep(lv, fit$strata)
  lower <- if (is.null(fit$lower)) fit$surv else fit$lower
  upper <- if (is.null(fit$upper)) fit$surv else fit$upper
  d <- data.frame(time = fit$time, surv = fit$surv, n_censor = fit$n.censor,
                  lower = lower, upper = upper, group = strata,
                  stringsAsFactors = FALSE)
  # every curve starts at (0, 1), otherwise the step plot floats
  anchors <- data.frame(time = 0, surv = 1, n_censor = 0, lower = 1, upper = 1,
                        group = lv, stringsAsFactors = FALSE)
  d <- rbind(anchors, d)
  d$group <- factor(d$group, levels = lv)
  d[order(d$group, d$time), , drop = FALSE]
}

#' Step-shaped confidence band from [km_tidy()] output
#' @keywords internal
km_band <- function(d) {
  do.call(rbind, lapply(split(d, d$group), function(x) {
    n <- nrow(x)
    if (n < 2) return(NULL)
    i <- c(1, rep(seq_len(n)[-1], each = 2))
    j <- c(rep(seq_len(n - 1), each = 2), n)
    data.frame(time = x$time[i], lower = x$lower[j], upper = x$upper[j],
               group = x$group[1])
  }))
}

#' Number at risk per group at the plot's time breaks
#' @keywords internal
km_risk <- function(fit, breaks) {
  sm <- summary(fit, times = breaks, extend = TRUE)
  grp <- if (is.null(sm$strata)) rep("All", length(sm$time))
         else sub("^\\.group=", "", as.character(sm$strata))
  lv <- if (is.null(fit$strata)) "All" else sub("^\\.group=", "", names(fit$strata))
  data.frame(time = sm$time, n = sm$n.risk, group = factor(grp, levels = lv))
}

#' Kaplan-Meier plot (ggplot2, no survminer)
#'
#' Curves with their 95% confidence bands, censoring ticks, the log-rank test
#' and -- when given -- the Cox hazard ratio with its interval, plus a
#' number-at-risk table under the panel (needs patchwork; skipped without it).
#'
#' @param fit A [km_fit()] result.
#' @param lr Optional [logrank_test()] result, annotated on the panel.
#' @param title Plot title (name the endpoint: "Overall survival", "PFS"...).
#' @param time_label X axis label.
#' @param hr Optional [cox_hr()] result, annotated under the log-rank p.
#' @param ci Draw the 95% confidence bands.
#' @param risk_table Add the number-at-risk table.
#' @param note Optional subtitle (e.g. "Exploratory: optimal cut-point").
#' @return A ggplot (or patchwork) object.
#' @keywords internal
km_plot <- function(fit, lr = NULL, title = "Survival", time_label = "Months",
                    hr = NULL, ci = TRUE, risk_table = TRUE, note = NULL, colors = NULL) {
  d <- km_tidy(fit)
  groups <- levels(d$group)
  cols <- if (!is.null(colors) && all(groups %in% names(colors))) colors[groups] else km_group_colors(groups)
  cens <- d[d$n_censor > 0, , drop = FALSE]
  breaks <- pretty(c(0, max(d$time, na.rm = TRUE)), n = 5)
  breaks <- breaks[breaks <= max(d$time, na.rm = TRUE)]

  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$time, y = .data$surv,
                                       colour = .data$group))
  if (isTRUE(ci)) {
    band <- km_band(d)
    if (!is.null(band) && nrow(band)) {
      p <- p + ggplot2::geom_ribbon(
        data = band, inherit.aes = FALSE, alpha = 0.15, show.legend = FALSE,
        ggplot2::aes(x = .data$time, ymin = .data$lower, ymax = .data$upper,
                     fill = .data$group))
    }
  }
  p <- p +
    ggplot2::geom_step(linewidth = 0.9) +
    ggplot2::scale_colour_manual(values = cols, name = NULL) +
    ggplot2::scale_fill_manual(values = cols, guide = "none") +
    ggplot2::scale_x_continuous(breaks = breaks, limits = range(c(0, d$time))) +
    ggplot2::scale_y_continuous(limits = c(0, 1), labels = function(v) paste0(v * 100, "%")) +
    ggplot2::labs(x = time_label, y = "Survival probability", title = title,
                  subtitle = note) +
    omicone_theme()

  if (nrow(cens)) {
    p <- p + ggplot2::geom_point(data = cens, shape = 3, size = 1.8,
                                 show.legend = FALSE)
  }
  if (length(groups) < 2) p <- p + ggplot2::theme(legend.position = "none")

  labs <- character(0)
  if (!is.null(lr) && is.finite(lr$p)) {
    labs <- c(labs, if (lr$p < 0.0001) "log-rank p < 0.0001"
                    else sprintf("log-rank p = %.3g", lr$p))
  }
  if (!is.null(hr) && is.finite(hr$hr)) {
    # ASCII only: the PDF device's Latin-1 encoding has no en dash
    labs <- c(labs, sprintf("HR (%s) %.2f, 95%% CI %.2f-%.2f", hr$label,
                            hr$hr, hr$lower, hr$upper))
  }
  if (length(labs)) {
    p <- p + ggplot2::annotate("text", x = 0, y = 0.03, hjust = 0, vjust = 0,
                               label = paste(labs, collapse = "\n"),
                               size = 3.8, colour = style_tokens()$muted, lineheight = 1.1)
  }

  if (!isTRUE(risk_table) || !length(breaks)) return(p)
  rk <- tryCatch(km_risk(fit, breaks), error = function(e) NULL)
  if (is.null(rk) || !nrow(rk)) return(p)
  tbl <- ggplot2::ggplot(rk, ggplot2::aes(x = .data$time, y = .data$group,
                                          label = .data$n, colour = .data$group)) +
    ggplot2::geom_text(size = 3.6, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = cols, guide = "none") +
    ggplot2::scale_x_continuous(breaks = breaks, limits = range(c(0, d$time))) +
    ggplot2::scale_y_discrete(limits = rev(groups)) +
    ggplot2::labs(x = NULL, y = NULL, title = "Number at risk") +
    omicone_theme() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size = 10, face = "plain"))
  compose_grid(list(p, tbl), nrow = 2, ncol = 1, heights = c(4, 1.1))
}

#' Colours of survival groups: fixed meanings where the labels have them
#'
#' WT / Mutant, Low / High (a split of a continuous variable) and TMB low /
#' high get the same colours in every survival figure; anything else follows
#' [value_colors()].
#' @param groups Group labels in plot order.
#' @keywords internal
km_group_colors <- function(groups) {
  tk <- style_tokens()
  low <- tolower(groups)
  is_low <- low == "wt" | grepl("^(tmb )?low", low)
  is_high <- low == "mutant" | grepl("^(tmb )?high", low)
  if (length(groups) == 2 && sum(is_low) == 1 && sum(is_high) == 1) {
    return(stats::setNames(ifelse(is_low, tk$wt, tk$mutant), groups))
  }
  value_colors(groups, levels = groups)
}
