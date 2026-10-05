#' WES variant filters and hypermutator flags (pure functions)
#'
#' A MAF from a caller still holds calls a careful analysis would not keep:
#' records the caller itself flagged (FILTER != PASS), calls supported by a
#' handful of reads, very low allele fractions, and common germline
#' polymorphisms that leaked through (high population allele frequency).
#' [wes_filter_apply()] applies these filters in a fixed order and counts what
#' each one removed, so the funnel can be reported in a methods section.
#'
#' Policy for missing values: a variant whose filter column is empty is kept
#' (the filter cannot judge it) and counted as "not evaluable"; a missing
#' population frequency means "not seen in the population database" (0).
#'
#' @name fct_wes_filter
#' @keywords internal
NULL

#' Population allele-frequency columns, most specific first
#' @keywords internal
wes_pop_af_cols <- function() {
  c("gnomAD_AF", "gnomADe_AF", "gnomAD_exomes_AF", "gnomAD_exome_AF", "gnomADg_AF",
    "gnomAD_genome_AF", "MAX_AF", "ExAC_AF", "ExAC_AF_Adj", "1000G_AF", "AF_1000G",
    "ESP_AF", "AF")
}

#' Which filter columns a MAF offers
#' @param fields MAF field names ([wes_fields()]).
#' @return list(filter, depth, alt, ref, vaf, pop, normal_alt): a column name or
#'   NULL each. `depth` may be NULL while `ref` + `alt` exist (depth = their sum).
#' @keywords internal
wes_filter_cols <- function(fields) {
  pick <- function(x) {
    hit <- x[x %in% fields]
    if (length(hit)) hit[1] else NULL
  }
  vaf <- wes_pick_vaf_col(fields)
  list(filter = pick(c("FILTER", "filter")),
       depth = pick(c("t_depth", "tumor_depth", "t_dp")),
       alt = pick(c("t_alt_count", "tumor_alt_count", "t_ad_alt")),
       ref = pick(c("t_ref_count", "tumor_ref_count")),
       vaf = vaf,
       pop = pick(wes_pop_af_cols()),
       normal_alt = pick(c("n_alt_count", "normal_alt_count")))
}

# numeric version of a MAF column (".", "" and text become NA)
wes_num_col <- function(dt, col) {
  if (is.null(col) || !col %in% names(dt)) return(rep(NA_real_, nrow(dt)))
  suppressWarnings(as.numeric(as.character(dt[[col]])))
}

#' Per-variant values the filters look at
#' @param dt Variant table (data.frame / data.table). @param cols [wes_filter_cols()].
#' @return data.frame(pass, depth, alt, vaf, pop, normal_alt) aligned to `dt`.
#' @keywords internal
wes_filter_values <- function(dt, cols) {
  n <- nrow(dt)
  alt <- wes_num_col(dt, cols$alt)
  depth <- wes_num_col(dt, cols$depth)
  if (is.null(cols$depth) && !is.null(cols$ref) && !is.null(cols$alt)) depth <- wes_num_col(dt, cols$ref) + alt
  vaf <- wes_num_col(dt, cols$vaf)
  if (!is.null(cols$vaf) && any(vaf > 1, na.rm = TRUE)) vaf <- vaf / 100          # percent
  if (is.null(cols$vaf)) vaf <- ifelse(depth > 0, alt / depth, NA_real_)
  pop <- wes_num_col(dt, cols$pop)
  pop[is.na(pop)] <- 0
  pass <- if (is.null(cols$filter)) rep(NA, n) else {
    f <- as.character(dt[[cols$filter]])
    ifelse(is.na(f) | f == "", NA, f %in% c("PASS", "."))
  }
  data.frame(pass = pass, depth = depth, alt = alt, vaf = vaf, pop = pop,
             normal_alt = wes_num_col(dt, cols$normal_alt))
}

#' Apply the variant filters to a MAF
#'
#' @param maf A MAF object (the imported one).
#' @param pass Keep only FILTER = PASS (or ".").
#' @param min_depth,min_alt,min_vaf,max_pop,max_normal_alt Thresholds; NA = off.
#'   `min_vaf` and `max_pop` are fractions (0.05, 0.001).
#' @param cols Column choice (default [wes_filter_cols()] on the MAF's fields).
#' @return list(maf = filtered MAF, funnel = data.frame(step, filter, column,
#'   removed, removed_nonsyn, not_evaluable, nonsyn_left, all_left,
#'   samples_left), lost = samples left with no variant record at all, cols =,
#'   values = per-variant values with `kept` / `nonsyn` flags, vc = the
#'   non-synonymous classes, unavailable = filters asked for whose column the
#'   MAF lacks (switched off), thresholds = what actually ran).
#' @keywords internal
wes_filter_apply <- function(maf, pass = TRUE, min_depth = NA, min_alt = NA, min_vaf = NA,
                             max_pop = NA, max_normal_alt = NA, cols = NULL) {
  cols <- cols %||% wes_filter_cols(wes_fields(maf))
  # a filter whose column is not in this MAF is switched off (and reported)
  has_depth <- !is.null(cols$depth) || (!is.null(cols$ref) && !is.null(cols$alt))
  avail <- c(pass = !is.null(cols$filter), min_depth = has_depth, min_alt = !is.null(cols$alt),
             min_vaf = !is.null(cols$vaf) || (has_depth && !is.null(cols$alt)),
             max_pop = !is.null(cols$pop), max_normal_alt = !is.null(cols$normal_alt))
  asked <- c(pass = isTRUE(pass), min_depth = is.finite(min_depth), min_alt = is.finite(min_alt),
             min_vaf = is.finite(min_vaf), max_pop = is.finite(max_pop),
             max_normal_alt = is.finite(max_normal_alt))
  unavailable <- names(asked)[asked & !avail]
  if (!avail[["pass"]]) pass <- FALSE
  if (!avail[["min_depth"]]) min_depth <- NA
  if (!avail[["min_alt"]]) min_alt <- NA
  if (!avail[["min_vaf"]]) min_vaf <- NA
  if (!avail[["max_pop"]]) max_pop <- NA
  if (!avail[["max_normal_alt"]]) max_normal_alt <- NA
  ns <- as.data.frame(maf@data)
  sil <- as.data.frame(maf@maf.silent)
  common <- union(names(ns), names(sil))
  for (cl in setdiff(common, names(ns))) ns[[cl]] <- NA
  for (cl in setdiff(common, names(sil))) sil[[cl]] <- NA
  dt <- rbind(ns[, common, drop = FALSE], sil[, common, drop = FALSE])
  nonsyn <- c(rep(TRUE, nrow(ns)), rep(FALSE, nrow(sil)))
  v <- wes_filter_values(dt, cols)
  keep <- rep(TRUE, nrow(dt))
  samples0 <- unique(as.character(dt$Tumor_Sample_Barcode))
  rows <- list()
  add_row <- function(label, col, test) {
    if (is.null(test)) return(invisible(NULL))
    ok <- test
    ne <- is.na(ok) & keep
    ok[is.na(ok)] <- TRUE                                        # cannot judge: keep
    removed <- keep & !ok
    keep <<- keep & ok
    rows[[length(rows) + 1]] <<- data.frame(
      filter = label, column = col, removed = sum(removed), removed_nonsyn = sum(removed & nonsyn),
      not_evaluable = sum(ne), nonsyn_left = sum(keep & nonsyn), all_left = sum(keep),
      samples_left = length(unique(as.character(dt$Tumor_Sample_Barcode[keep]))),
      stringsAsFactors = FALSE)
  }
  rows[[1]] <- data.frame(filter = "imported", column = "", removed = 0L, removed_nonsyn = 0L,
                          not_evaluable = 0L, nonsyn_left = sum(nonsyn), all_left = nrow(dt),
                          samples_left = length(samples0),
                          stringsAsFactors = FALSE)
  if (isTRUE(pass)) add_row("FILTER = PASS", cols$filter, v$pass)
  if (is.finite(min_depth)) add_row(sprintf("tumour depth >= %g", min_depth),
                                    cols$depth %||% paste(cols$ref, "+", cols$alt), v$depth >= min_depth)
  if (is.finite(min_alt)) add_row(sprintf("alt reads >= %g", min_alt), cols$alt %||% "", v$alt >= min_alt)
  if (is.finite(min_vaf)) add_row(sprintf("VAF >= %g", min_vaf), cols$vaf %||% "alt / depth", v$vaf >= min_vaf)
  if (is.finite(max_pop)) add_row(sprintf("population AF <= %g", max_pop), cols$pop, v$pop <= max_pop)
  if (is.finite(max_normal_alt)) add_row(sprintf("normal alt reads <= %g", max_normal_alt),
                                         cols$normal_alt %||% "", v$normal_alt <= max_normal_alt)
  funnel <- do.call(rbind, rows)
  funnel$step <- seq_len(nrow(funnel)) - 1L
  if (!any(keep & nonsyn)) stop("No non-synonymous variant passes these filters; loosen them.")
  vc <- unique(as.character(ns$Variant_Classification))
  out <- maftools::read.maf(dt[keep, , drop = FALSE], clinicalData = maf@clinical.data,
                            vc_nonSyn = vc, verbose = FALSE)
  left <- unique(as.character(dt$Tumor_Sample_Barcode[keep]))
  v$kept <- keep
  v$nonsyn <- nonsyn
  list(maf = out, funnel = funnel[, c("step", setdiff(names(funnel), "step"))],
       lost = setdiff(samples0, left), cols = cols, values = v, vc = vc, unavailable = unavailable,
       thresholds = list(pass = isTRUE(pass), min_depth = min_depth,
                         min_alt = min_alt, min_vaf = min_vaf, max_pop = max_pop,
                         max_normal_alt = max_normal_alt))
}

#' Flag hypermutated samples (Tukey's far-out fence on log10 counts)
#'
#' A sample is flagged when its non-synonymous count lies above Q3 + 3 IQR of
#' the cohort's log10 counts. A flag, not an exclusion: such samples (MMR / POLE
#' deficiency, temozolomide, artefacts) dominate frequency and enrichment
#' results, so analyses are usually repeated without them.
#' @param maf A MAF object.
#' @return data.frame(Tumor_Sample_Barcode, nonsyn, hypermutated) sorted by
#'   count; attribute "fence" (count threshold).
#' @keywords internal
wes_hypermutators <- function(maf) {
  ss <- as.data.frame(maftools::getSampleSummary(maf))
  n <- as.numeric(ss$total)
  lg <- log10(pmax(n, 1))
  q <- stats::quantile(lg, c(0.25, 0.75), names = FALSE)
  fence <- 10^(q[2] + 3 * (q[2] - q[1]))
  out <- data.frame(Tumor_Sample_Barcode = as.character(ss$Tumor_Sample_Barcode), nonsyn = n,
                    hypermutated = n > fence, stringsAsFactors = FALSE)
  out <- out[order(-out$nonsyn), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "fence") <- fence
  out
}

#' Runnable R code for the filters that ran (on `maf` from the import log)
#' @param f Output of [wes_filter_apply()].
#' @keywords internal
wes_filter_code <- function(f) {
  th <- f$thresholds
  cols <- f$cols
  num <- function(col) sprintf("suppressWarnings(as.numeric(dt[[%s]]))", r_lit(col))
  depth <- if (!is.null(cols$depth)) num(cols$depth) else sprintf("%s + %s", num(cols$ref), num(cols$alt))
  vaf <- if (!is.null(cols$vaf)) {
    sprintf("vaf <- %s; if (any(vaf > 1, na.rm = TRUE)) vaf <- vaf / 100", num(cols$vaf))
  } else "vaf <- alt / depth"
  okna <- function(cond) sprintf("keep <- keep & (is.na(%s) | %s)", sub(" .*$", "", cond), cond)
  c("# variant filters; a variant whose filter column is empty is kept",
    "dt <- as.data.frame(data.table::rbindlist(list(maf@data, maf@maf.silent), fill = TRUE))",
    "keep <- rep(TRUE, nrow(dt))",
    if (th$pass) sprintf("keep <- keep & (is.na(dt[[%s]]) | dt[[%s]] %%in%% c(\"PASS\", \".\", \"\"))",
                         r_lit(cols$filter), r_lit(cols$filter)),
    if (is.finite(th$min_depth) || is.finite(th$min_alt) || is.finite(th$min_vaf)) {
      c(if (!is.null(cols$alt)) sprintf("alt <- %s", num(cols$alt)),
        if (!is.null(cols$depth) || !is.null(cols$ref)) sprintf("depth <- %s", depth))
    },
    if (is.finite(th$min_depth)) okna(sprintf("depth >= %s", r_lit(th$min_depth))),
    if (is.finite(th$min_alt)) okna(sprintf("alt >= %s", r_lit(th$min_alt))),
    if (is.finite(th$min_vaf)) c(vaf, okna(sprintf("vaf >= %s", r_lit(th$min_vaf)))),
    if (is.finite(th$max_pop) && !is.null(cols$pop)) {
      c(sprintf("pop <- %s; pop[is.na(pop)] <- 0   # absent from the database = 0", num(cols$pop)),
        sprintf("keep <- keep & pop <= %s", r_lit(th$max_pop)))
    },
    if (is.finite(th$max_normal_alt)) c(sprintf("nalt <- %s", num(cols$normal_alt)),
                                        okna(sprintf("nalt <= %s", r_lit(th$max_normal_alt)))),
    sprintf("maf <- maftools::read.maf(dt[keep, ], clinicalData = maf@clinical.data, vc_nonSyn = %s)",
            r_lit(f$vc)))
}

#' Funnel of variants left after each filter
#' @keywords internal
wes_filter_funnel_plot <- function(funnel) {
  d <- funnel
  d$label <- factor(d$filter, levels = rev(d$filter))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$nonsyn_left, y = .data$label)) +
    ggplot2::geom_col(fill = style_tokens()$accent, width = 0.65) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%s%s", format(.data$nonsyn_left, big.mark = ","),
                                                    ifelse(.data$removed_nonsyn > 0,
                                                           sprintf("  (-%s)", format(.data$removed_nonsyn, big.mark = ",")),
                                                           ""))),
                       hjust = -0.08, size = 3.4, colour = style_tokens()$text) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.28))) +
    ggplot2::labs(x = "non-synonymous variants left", y = NULL) +
    omicone_theme()
}

#' VAF and depth of the variants, kept vs removed
#' @keywords internal
wes_filter_vaf_plot <- function(values, min_vaf = NA) {
  d <- values[values$nonsyn & !is.na(values$vaf), , drop = FALSE]
  if (!nrow(d)) stop("No VAF available in this MAF.")
  d$status <- factor(ifelse(d$kept, "kept", "removed"), levels = c("kept", "removed"))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$vaf, fill = .data$status)) +
    ggplot2::geom_histogram(bins = 50, boundary = 0, colour = NA, alpha = 0.9) +
    ggplot2::scale_fill_manual(values = c(kept = style_tokens()$accent, removed = style_tokens()$removed), name = NULL, drop = FALSE) +
    ggplot2::labs(x = "variant allele fraction (tumour)", y = "non-synonymous variants") +
    omicone_theme()
  if (is.finite(min_vaf)) p <- p + ggplot2::geom_vline(xintercept = min_vaf, linetype = 2, colour = style_tokens()$muted)
  p
}

#' Non-synonymous count per sample (log scale), with the hypermutator fence
#' @keywords internal
wes_hyper_plot <- function(hm) {
  d <- hm[order(hm$nonsyn), , drop = FALSE]
  d$rank <- seq_len(nrow(d))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$rank, y = pmax(.data$nonsyn, 0.8), colour = .data$hypermutated)) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_hline(yintercept = attr(hm, "fence"), linetype = 2, colour = style_tokens()$removed) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_colour_manual(values = c(`FALSE` = style_tokens()$accent, `TRUE` = style_tokens()$removed),
                                 labels = c(`FALSE` = "within range", `TRUE` = "hypermutated (flag)"), name = NULL) +
    ggplot2::labs(x = "samples, sorted", y = "non-synonymous variants (log)",
                  caption = sprintf("dashed: Tukey far-out fence on log10 counts (Q3 + 3 IQR) = %s",
                                    format(round(attr(hm, "fence")), big.mark = ","))) +
    omicone_theme()
}
