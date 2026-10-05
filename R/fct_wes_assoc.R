#' TMB versus outcome: survival and response (pure functions)
#'
#' TMB is a continuous, very skewed biomarker. The primary analysis therefore
#' uses it on a log2 scale with no cut-off (AGENTS.md section 7.10): a Cox
#' hazard ratio, or a logistic odds ratio, per doubling of TMB; for a binary
#' response also the ROC AUC. A median split is drawn only as an illustration.
#' A pseudo-count of one mutation (1 / capture Mb) keeps TMB-0 samples on the
#' log scale.
#'
#' @name fct_wes_assoc
#' @keywords internal
NULL

#' Per-sample TMB, including sequenced samples with no variant
#' @param maf A MAF object. @param capture Capture size (Mb).
#' @param universe Clinical ids taken as sequenced (TMB 0, `in_maf` FALSE when absent).
#' @param sequenced Samples of the imported MAF (`rv$wes_sequenced`): those the
#'   variant filters emptied are in the MAF with TMB 0.
#' @param tcga12 Join on 12-character TCGA barcodes.
#' @return data.frame(.id, nonsyn, tmb, in_maf).
#' @keywords internal
wes_tmb_table <- function(maf, capture, universe = NULL, sequenced = NULL, tcga12 = FALSE) {
  cap <- suppressWarnings(as.numeric(capture))
  if (length(cap) != 1 || !is.finite(cap) || cap <= 0) stop("Capture size must be a positive number of megabases.")
  ss <- as.data.frame(maftools::getSampleSummary(maf))
  ids <- wes_norm_id(ss$Tumor_Sample_Barcode, tcga12)
  n <- tapply(as.numeric(ss$total), ids, mean)           # several samples under one id: averaged
  out <- data.frame(.id = names(n), nonsyn = as.numeric(n), in_maf = TRUE, stringsAsFactors = FALSE)
  filtered <- setdiff(unique(wes_norm_id(sequenced, tcga12)), c(out$.id, NA, ""))
  if (length(filtered)) out <- rbind(out, data.frame(.id = filtered, nonsyn = 0, in_maf = TRUE, stringsAsFactors = FALSE))
  extra <- setdiff(unique(wes_norm_id(universe, tcga12)), c(out$.id, NA, ""))
  if (length(extra)) out <- rbind(out, data.frame(.id = extra, nonsyn = 0, in_maf = FALSE, stringsAsFactors = FALSE))
  out$tmb <- out$nonsyn / cap
  attr(out, "capture") <- cap
  out
}

#' Join TMB to the clinical table: the per-patient analysis set
#'
#' @param clin Clinical data.frame. @param tmb [wes_tmb_table()].
#' @param id_col Sample id column of `clin`.
#' @param outcome "survival" (time + event) or "response" (binary).
#' @param time_col,event_col,time_unit Survival columns.
#' @param resp_col,positive Response column and the level meaning "responder".
#' @param unmatched_zero Keep clinical ids without a MAF record (TMB 0).
#' @param tcga12 Join on 12-character TCGA barcodes.
#' @return data.frame with `.id`, `tmb`, `x` (log2 TMB with the pseudo-count)
#'   and `.time` / `.event` or `.resp`; attribute "flow".
#' @keywords internal
wes_tmb_assoc_data <- function(clin, tmb, id_col, outcome = c("survival", "response"),
                               time_col = NULL, event_col = NULL, time_unit = "days",
                               resp_col = NULL, positive = NULL, unmatched_zero = TRUE,
                               tcga12 = FALSE) {
  outcome <- match.arg(outcome)
  clin <- as.data.frame(clin, stringsAsFactors = FALSE)
  for (cl in c(id_col, time_col, event_col, resp_col)) {
    if (!cl %in% names(clin)) stop("Column not found in the clinical table: ", cl)
  }
  cap <- attr(tmb, "capture")
  clin$.sid <- wes_norm_id(clin[[id_col]], tcga12)
  n_rows <- nrow(clin)
  clin <- clin[!is.na(clin$.sid) & nzchar(clin$.sid), , drop = FALSE]
  matched <- clin$.sid %in% tmb$.id[tmb$in_maf]
  if (!isTRUE(unmatched_zero)) clin <- clin[matched, , drop = FALSE]
  clin$tmb <- tmb$tmb[match(clin$.sid, tmb$.id)]
  clin$tmb[is.na(clin$tmb)] <- 0                          # sequenced, no variant record
  if (outcome == "survival") {
    tt <- suppressWarnings(as.numeric(clin[[time_col]]))
    tt[!is.finite(tt)] <- NA_real_
    clin[[time_col]] <- tt
    d <- normalise_clinical(clin, ".sid", time_col, event_col, time_unit = time_unit, dedup = "first")
  } else {
    d <- clin
    d$.id <- d$.sid
    lv <- as.character(d[[resp_col]])
    d$.resp <- ifelse(is.na(lv) | !nzchar(trimws(lv)), NA_integer_, as.integer(lv == as.character(positive)))
    d <- d[!is.na(d$.resp), , drop = FALSE]
    dup <- duplicated(d$.id)
    d <- d[!dup, , drop = FALSE]
  }
  d$x <- log2(d$tmb + 1 / cap)
  attr(d, "flow") <- list(clin_rows = n_rows, matched = sum(matched), unmatched = sum(!matched),
                          unmatched_zero = isTRUE(unmatched_zero), n = nrow(d),
                          events = if (outcome == "survival") sum(d$.event) else sum(d$.resp),
                          capture = cap)
  d
}

#' Cox model of survival on log2 TMB (per doubling), with a PH check
#' @param d [wes_tmb_assoc_data()] for "survival".
#' @return list(hr, lower, upper, p, n, events, zph_p, note), or NULL.
#' @keywords internal
wes_tmb_cox <- function(d) {
  if (nrow(d) < 5 || stats::sd(d$x) == 0) return(NULL)
  notes <- character(0)
  fit <- withCallingHandlers(
    tryCatch(survival::coxph(survival::Surv(.time, .event) ~ x, data = d), error = function(e) NULL),
    warning = function(w) {
      notes <<- c(notes, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  if (is.null(fit)) return(NULL)
  sm <- summary(fit)
  z <- tryCatch(survival::cox.zph(fit), error = function(e) NULL)
  list(hr = unname(sm$conf.int[1, "exp(coef)"]), lower = unname(sm$conf.int[1, "lower .95"]),
       upper = unname(sm$conf.int[1, "upper .95"]), p = unname(sm$coefficients[1, "Pr(>|z|)"]),
       n = nrow(d), events = sum(d$.event),
       zph_p = if (is.null(z)) NA_real_ else unname(z$table[1, "p"]), note = unique(notes))
}

#' AUC with a 95% CI (DeLong via pROC when installed, else Hanley-McNeil)
#' @param score Numeric predictor. @param y 0/1 outcome.
#' @keywords internal
wes_auc <- function(score, y) {
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  r <- rank(score)
  auc <- (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  if (has_pkg("pROC")) {
    roc <- pROC::roc(y, score, levels = c(0, 1), direction = "<", quiet = TRUE)
    ci <- as.numeric(pROC::ci.auc(roc, method = "delong"))
    return(list(auc = ci[2], lower = ci[1], upper = ci[3], method = "DeLong (pROC)"))
  }
  q1 <- auc / (2 - auc)
  q2 <- 2 * auc^2 / (1 + auc)
  se <- sqrt((auc * (1 - auc) + (n1 - 1) * (q1 - auc^2) + (n0 - 1) * (q2 - auc^2)) / (n1 * n0))
  list(auc = auc, lower = max(0, auc - 1.96 * se), upper = min(1, auc + 1.96 * se), method = "Hanley-McNeil")
}

#' Logistic regression of response on log2 TMB, plus the ROC AUC and a rank test
#' @param d [wes_tmb_assoc_data()] for "response".
#' @return list(or, lower, upper, p, auc, auc_lower, auc_upper, auc_method,
#'   wilcox_p, n, n_resp, n_non, median_resp, median_non, roc = data.frame(fpr, tpr)).
#' @keywords internal
wes_tmb_logit <- function(d) {
  y <- d$.resp
  if (length(unique(y)) < 2) stop("The response column has a single level in the analysis set.")
  fit <- suppressWarnings(stats::glm(.resp ~ x, data = d, family = stats::binomial()))
  co <- summary(fit)$coefficients
  ci <- suppressMessages(stats::confint.default(fit)["x", ])
  a <- wes_auc(d$x, y)
  thr <- sort(unique(d$x), decreasing = TRUE)
  roc <- rbind(data.frame(fpr = 0, tpr = 0),
               data.frame(fpr = vapply(thr, function(t) mean(d$x[y == 0] >= t), 0),
                          tpr = vapply(thr, function(t) mean(d$x[y == 1] >= t), 0)))
  list(or = exp(co["x", "Estimate"]), lower = exp(ci[[1]]), upper = exp(ci[[2]]),
       p = co["x", "Pr(>|z|)"], auc = a$auc, auc_lower = a$lower, auc_upper = a$upper,
       auc_method = a$method,
       wilcox_p = suppressWarnings(stats::wilcox.test(d$tmb[y == 1], d$tmb[y == 0])$p.value),
       n = nrow(d), n_resp = sum(y == 1), n_non = sum(y == 0),
       median_resp = stats::median(d$tmb[y == 1]), median_non = stats::median(d$tmb[y == 0]),
       roc = roc)
}

#' Small-sample warnings for the TMB association
#' @keywords internal
wes_tmb_assoc_warnings <- function(d, outcome) {
  out <- character(0)
  if (outcome == "survival") {
    if (sum(d$.event) < 10) out <- c(out, sprintf("only %d events", sum(d$.event)))
  } else {
    k <- min(sum(d$.resp == 1), sum(d$.resp == 0))
    if (k < 10) out <- c(out, sprintf("only %d patients in the smaller response group", k))
  }
  if (mean(d$tmb == 0) > 0.25) out <- c(out, sprintf("%.0f%% of patients have TMB 0", 100 * mean(d$tmb == 0)))
  out
}

#' TMB by response group (log scale), with the rank-test p-value
#' @keywords internal
wes_tmb_resp_plot <- function(d, res, positive) {
  d$group <- factor(ifelse(d$.resp == 1, positive, "other"), levels = c("other", positive))
  floor <- 0.5 / attr(d, "flow")$capture
  ggplot2::ggplot(d, ggplot2::aes(x = .data$group, y = pmax(.data$tmb, floor))) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.5, colour = style_tokens()$faint) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$group),
                        position = ggplot2::position_jitter(width = 0.13, height = 0, seed = 1), size = 1.8, alpha = 0.8) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_colour_manual(values = c(style_tokens()$faint, style_tokens()$accent), guide = "none") +
    ggplot2::labs(x = NULL, y = "TMB (mut/Mb, log; TMB 0 drawn at the floor)",
                  subtitle = sprintf("Wilcoxon rank-sum p = %s; OR per doubling %.2f (%.2f-%.2f)",
                                     format(signif(res$wilcox_p, 3)), res$or, res$lower, res$upper)) +
    omicone_theme()
}

#' ROC curve of log2 TMB for the response
#' @keywords internal
wes_tmb_roc_plot <- function(res) {
  ggplot2::ggplot(res$roc, ggplot2::aes(x = .data$fpr, y = .data$tpr)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = style_tokens()$faint) +
    ggplot2::geom_step(colour = style_tokens()$accent, linewidth = 1, direction = "vh") +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "1 - specificity", y = "sensitivity",
                  subtitle = sprintf("AUC %.2f (95%% CI %.2f-%.2f, %s); %d responders, %d others",
                                     res$auc, res$auc_lower, res$auc_upper, res$auc_method,
                                     res$n_resp, res$n_non)) +
    omicone_theme()
}

#' Runnable R code for the TMB association (after the import log's `maf`)
#' @keywords internal
wes_tmb_assoc_code <- function(p) {
  head <- c(
    sprintf("capture <- %s", r_lit(p$capture)),
    "ss <- maftools::getSampleSummary(maf)",
    sprintf("tmb <- setNames(ss$total / capture, toupper(trimws(ss$Tumor_Sample_Barcode)))"),
    "clin <- clin_raw   # the clinical table read in the WES import step",
    sprintf("clin$.id <- toupper(trimws(clin[[%s]]))", r_lit(p$id_col)),
    if (!isTRUE(p$unmatched_zero)) "clin <- clin[clin$.id %in% names(tmb), ]",
    "clin$tmb <- ifelse(clin$.id %in% names(tmb), tmb[clin$.id], 0)   # sequenced, no variant = 0",
    "clin$x <- log2(clin$tmb + 1 / capture)   # per doubling; + 1 mutation keeps TMB 0")
  if (identical(p$outcome, "survival")) {
    conv <- switch(p$time_unit, days = " / 30.4375", years = " * 12", "")
    c(head,
      sprintf("clin$.time <- suppressWarnings(as.numeric(clin[[%s]]))%s   # months", r_lit(p$time_col), conv),
      p$event_code,
      "clin <- clin[is.finite(clin$.time) & !is.na(clin$.event) & !duplicated(clin$.id), ]",
      "fit <- survival::coxph(survival::Surv(.time, .event) ~ x, data = clin)",
      "summary(fit)   # HR per doubling of TMB",
      "survival::cox.zph(fit)")
  } else {
    c(head,
      sprintf("clin$.resp <- as.integer(as.character(clin[[%s]]) == %s)", r_lit(p$resp_col), r_lit(p$positive)),
      sprintf("clin <- clin[!is.na(clin[[%s]]) & nzchar(trimws(as.character(clin[[%s]]))) & !duplicated(clin$.id), ]",
              r_lit(p$resp_col), r_lit(p$resp_col)),
      "fit <- glm(.resp ~ x, data = clin, family = binomial())",
      "exp(cbind(OR = coef(fit), confint.default(fit)))[\"x\", ]   # OR per doubling of TMB",
      "pROC::ci.auc(pROC::roc(clin$.resp, clin$x, levels = c(0, 1), direction = \"<\"), method = \"delong\")",
      "wilcox.test(tmb ~ .resp, data = clin)")
  }
}
