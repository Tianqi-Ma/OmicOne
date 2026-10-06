#' Sample-level analyses: pseudobulk differential expression, differential abundance
#'
#' Comparisons between conditions (responders vs non-responders, treated vs
#' untreated) must use the *sample* (patient) as the unit of replication: cells
#' from one patient are not independent replicates, and cell-level tests on
#' them give wildly over-confident p-values (pseudo-replication; Squair et al.
#' 2021, Nat Commun). Both analyses here first aggregate cells per sample, then
#' test across samples:
#'
#' * **pseudobulk DE** sums the counts of one cell type per sample and tests
#'   with edgeR (quasi-likelihood), DESeq2 or limma-voom;
#' * **differential abundance** compares cell-type proportions per sample with
#'   the propeller method (Phipson et al. 2022, Bioinformatics): logit
#'   proportions with a 0.5 pseudo-count, limma linear model, robust empirical
#'   Bayes, BH across cell types.
#'
#' All functions are pure (no Shiny) and take a plain cell metadata frame, so
#' they are testable without Seurat.
#'
#' @name fct_sample
#' @keywords internal
NULL

#' One row per sample, checking that the condition is a property of the sample
#'
#' @param md Cell metadata. @param sample_col,condition_col Column names.
#' @param covariate Optional sample-level column (a pairing id or a batch).
#' @return A data.frame with one row per sample: `sample`, `condition`,
#'   `n_cells` (and `covariate`). Stops with a clear message when a sample
#'   carries two values of the condition or of the covariate.
#' @keywords internal
sample_table <- function(md, sample_col, condition_col, covariate = NULL) {
  for (cl in c(sample_col, condition_col, covariate)) {
    if (!cl %in% names(md)) stop("Column not found in the cell metadata: ", cl)
  }
  if (identical(sample_col, condition_col)) stop("The sample and the condition must be different columns.")
  s <- as.character(md[[sample_col]])
  ok <- !is.na(s) & !is.na(md[[condition_col]])
  per_sample <- function(col) {
    per <- split(as.character(md[[col]])[ok], s[ok])
    bad <- names(per)[vapply(per, function(x) length(unique(x)) > 1, logical(1))]
    if (length(bad)) {
      stop(sprintf(paste("'%s' is not constant within %d sample(s) (e.g. %s): it must describe the",
                         "whole sample. Pick a patient / sample-level column."),
                   col, length(bad), paste(utils::head(bad, 3), collapse = ", ")))
    }
    vapply(per, `[`, "", 1)
  }
  cond <- per_sample(condition_col)
  out <- data.frame(sample = names(cond), condition = unname(cond),
                    n_cells = as.integer(table(s[ok])[names(cond)]), stringsAsFactors = FALSE)
  if (!is.null(covariate)) out$covariate <- unname(per_sample(covariate)[out$sample])
  out
}

#' Add sample-level columns from a sample sheet to the cell metadata (pure)
#'
#' @param md Cell metadata. @param sheet data.frame: the first column (or the
#'   column named like `sample_col`) holds sample ids, every other column is a
#'   sample-level variable. @param sample_col Sample column in `md`.
#' @return list(cols = data.frame of new columns aligned to the rows of `md`,
#'   matched = samples found in the sheet, missing = samples of `md` absent
#'   from it, unused = sheet rows matching no sample).
#' @keywords internal
apply_sample_sheet <- function(md, sheet, sample_col) {
  if (!sample_col %in% names(md)) stop("Column not found in the cell metadata: ", sample_col)
  if (ncol(sheet) < 2) stop("The sample sheet needs a sample column and at least one more column.")
  key <- if (sample_col %in% names(sheet)) sample_col else names(sheet)[1]
  ids <- as.character(sheet[[key]])
  if (anyDuplicated(ids)) stop("Sample ids are repeated in the sheet: ", paste(unique(ids[duplicated(ids)]), collapse = ", "))
  s <- as.character(md[[sample_col]])
  idx <- match(s, ids)
  vars <- setdiff(names(sheet), key)
  cols <- as.data.frame(lapply(sheet[vars], function(v) v[idx]), stringsAsFactors = FALSE,
                        check.names = FALSE)
  rownames(cols) <- rownames(md)
  have <- unique(s[!is.na(s)])
  list(cols = cols, matched = intersect(have, ids), missing = setdiff(have, ids),
       unused = setdiff(ids, have))
}

#' Sum counts per sample x cell type
#'
#' @param counts Genes x cells count matrix (dgCMatrix or dense).
#' @param md Cell metadata (rows in the same order as the columns of counts).
#' @param sample_col,group_col Sample and cell-type columns.
#' @param min_cells Pseudobulk columns built from fewer cells are dropped.
#' @return list(counts = genes x pseudobulk matrix, meta = data.frame with
#'   `sample`, `group`, `n_cells`), columns named "<sample>|<group>".
#' @keywords internal
pseudobulk_counts <- function(counts, md, sample_col, group_col, min_cells = 10) {
  s <- as.character(md[[sample_col]])
  g <- as.character(md[[group_col]])
  ok <- !is.na(s) & !is.na(g)
  key <- paste(s[ok], g[ok], sep = "|")
  f <- factor(key)
  first <- match(levels(f), key)
  ind <- Matrix::fac2sparse(f)                           # levels x cells
  pb <- counts[, ok, drop = FALSE] %*% Matrix::t(ind)    # genes x levels
  colnames(pb) <- levels(f)
  meta <- data.frame(column = levels(f), sample = s[ok][first], group = g[ok][first],
                     n_cells = as.integer(Matrix::rowSums(ind)), stringsAsFactors = FALSE)
  keep <- meta$n_cells >= min_cells
  list(counts = pb[, keep, drop = FALSE], meta = meta[keep, , drop = FALSE],
       dropped = meta[!keep, , drop = FALSE], min_cells = min_cells)
}

#' Pseudobulk differential expression for one cell type
#'
#' @param pb Output of [pseudobulk_counts()].
#' @param samples Output of [sample_table()].
#' @param group Cell type to test.
#' @param ref,alt Condition levels: `alt` versus `ref` (logFC > 0 = higher in `alt`).
#' @param method "edgeR" (quasi-likelihood F test), "DESeq2" (Wald) or "limma" (voom).
#' @param covariate Add the sample table's `covariate` column as an additive
#'   term (a pairing id for paired designs, or a batch).
#' @return data.frame: gene, logFC, aveExpr, stat, p, fdr (sorted by p);
#'   attributes "n" (samples per condition) and "method".
#' @keywords internal
pseudobulk_de <- function(pb, samples, group, ref, alt, method = c("edgeR", "DESeq2", "limma"),
                          covariate = FALSE) {
  method <- match.arg(method)
  m <- pb$meta[pb$meta$group == group, , drop = FALSE]
  st <- samples[match(m$sample, samples$sample), , drop = FALSE]
  use <- !is.na(st$condition) & st$condition %in% c(ref, alt)
  m <- m[use, , drop = FALSE]
  st <- st[use, , drop = FALSE]
  cond <- factor(st$condition, levels = c(ref, alt))
  n <- table(cond)
  if (any(n < 2)) {
    stop(sprintf("need at least 2 samples per condition with >= %d cells (have %s)",
                 pb$min_cells %||% 0, paste(sprintf("%s = %d", names(n), as.integer(n)), collapse = ", ")))
  }
  y <- as.matrix(pb$counts[, m$column, drop = FALSE])
  cov <- if (isTRUE(covariate)) factor(st$covariate) else NULL
  design <- if (is.null(cov)) stats::model.matrix(~ cond) else stats::model.matrix(~ cov + cond)
  if (qr(design)$rank < ncol(design)) stop("the covariate is confounded with the condition")
  if (nrow(design) <= ncol(design)) stop("no residual degrees of freedom (too few samples for this design)")
  coef <- ncol(design)

  if (method == "DESeq2") {
    if (!has_pkg("DESeq2")) stop("Package 'DESeq2' is required.")
    keep <- rowSums(y >= 10) >= min(n)
    cd <- data.frame(cond = cond, row.names = colnames(y))
    if (!is.null(cov)) cd$cov <- cov
    dds <- DESeq2::DESeqDataSetFromMatrix(round(y[keep, , drop = FALSE]), cd,
                                          if (is.null(cov)) ~ cond else ~ cov + cond)
    dds <- suppressMessages(DESeq2::DESeq(dds, quiet = TRUE))
    r <- as.data.frame(DESeq2::results(dds, contrast = c("cond", alt, ref)))
    out <- data.frame(gene = rownames(r), logFC = r$log2FoldChange, aveExpr = log2(r$baseMean + 1),
                      stat = r$stat, p = r$pvalue, fdr = r$padj, stringsAsFactors = FALSE)
  } else {
    if (!has_pkg("edgeR")) stop("Package 'edgeR' is required.")
    dge <- edgeR::DGEList(y)
    keep <- edgeR::filterByExpr(dge, design = design)
    if (sum(keep) < 2) stop("fewer than 2 genes pass the expression filter")
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    dge <- edgeR::calcNormFactors(dge)
    if (method == "edgeR") {
      dge <- edgeR::estimateDisp(dge, design)
      fit <- edgeR::glmQLFit(dge, design, robust = TRUE)
      tt <- edgeR::topTags(edgeR::glmQLFTest(fit, coef = coef), n = Inf, sort.by = "none")$table
      out <- data.frame(gene = rownames(tt), logFC = tt$logFC, aveExpr = tt$logCPM, stat = tt$F,
                        p = tt$PValue, fdr = tt$FDR, stringsAsFactors = FALSE)
    } else {
      if (!has_pkg("limma")) stop("Package 'limma' is required.")
      v <- limma::voom(dge, design)
      fit <- limma::eBayes(limma::lmFit(v, design), robust = TRUE)
      tt <- limma::topTable(fit, coef = coef, number = Inf, sort.by = "none")
      out <- data.frame(gene = rownames(tt), logFC = tt$logFC, aveExpr = tt$AveExpr, stat = tt$t,
                        p = tt$P.Value, fdr = tt$adj.P.Val, stringsAsFactors = FALSE)
    }
  }
  out <- out[order(out$p), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "n") <- n
  attr(out, "method") <- method
  out
}

#' Pseudobulk DE for several cell types
#'
#' Runs [pseudobulk_de()] per cell type. The FDR is computed within each cell
#' type (each is its own experiment, as in muscat / the OSCA book); a cell type
#' that cannot be tested is reported in `skipped` with the reason.
#' @inheritParams pseudobulk_de
#' @param groups Cell types to test (default: all).
#' @return list(table = data.frame with a leading `cell_type` column,
#'   overview = one row per tested cell type, skipped = data.frame(cell_type, reason)).
#' @keywords internal
pseudobulk_de_all <- function(pb, samples, ref, alt, method = "edgeR", covariate = FALSE,
                              groups = NULL, fdr = 0.05) {
  groups <- groups %||% sort(unique(pb$meta$group))
  tabs <- list()
  ns <- list()
  skipped <- data.frame(cell_type = character(0), reason = character(0), stringsAsFactors = FALSE)
  for (g in groups) {
    r <- tryCatch(pseudobulk_de(pb, samples, g, ref, alt, method = method, covariate = covariate),
                  error = function(e) conditionMessage(e))
    if (is.character(r)) {
      skipped[nrow(skipped) + 1, ] <- c(g, r)
      next
    }
    tabs[[g]] <- data.frame(cell_type = g, r, stringsAsFactors = FALSE)
    ns[[g]] <- attr(r, "n")
  }
  tab <- if (length(tabs)) do.call(rbind, unname(tabs)) else NULL
  overview <- if (length(tabs)) {
    do.call(rbind, lapply(names(tabs), function(g) {
      d <- tabs[[g]]
      sig <- !is.na(d$fdr) & d$fdr < fdr
      n <- ns[[g]]
      data.frame(cell_type = g, samples_ref = n[[1]], samples_alt = n[[2]], genes_tested = nrow(d),
                 up = sum(sig & d$logFC > 0), down = sum(sig & d$logFC < 0), stringsAsFactors = FALSE)
    }))
  }
  list(table = tab, overview = overview, skipped = skipped, ref = ref, alt = alt,
       method = method, fdr = fdr)
}

#' Differential abundance of cell types between conditions (propeller)
#'
#' The propeller method of Phipson et al. 2022 (speckle::propeller with
#' `transform = "logit"`), written with limma so speckle is not required:
#' cell-type proportions per sample with a 0.5 pseudo-count, logit transform,
#' a linear model with robust empirical Bayes; a moderated t test for two
#' conditions, a moderated F test for more; BH across cell types.
#'
#' @param md Cell metadata. @param sample_col,group_col,condition_col Columns.
#' @param samples Optional [sample_table()] (computed when NULL).
#' @return data.frame: group, one mean-proportion column per condition
#'   (`prop_<level>`), prop_ratio (two conditions: last / first level), stat,
#'   p, fdr (sorted by p); attributes "n", "props" (samples x groups) and
#'   "condition".
#' @keywords internal
abundance_test <- function(md, sample_col, group_col, condition_col, samples = NULL) {
  if (!has_pkg("limma")) stop("Package 'limma' is required.")
  if (is.null(samples)) samples <- sample_table(md, sample_col, condition_col)
  s <- as.character(md[[sample_col]])
  g <- as.character(md[[group_col]])
  ok <- !is.na(s) & !is.na(g) & s %in% samples$sample
  tab <- unclass(table(factor(s[ok], levels = samples$sample), g[ok]))    # samples x groups
  if (ncol(tab) < 2) stop("Differential abundance needs at least 2 cell types.")
  props <- tab / rowSums(tab)
  pseudo <- (tab + 0.5) / rowSums(tab + 0.5)
  y <- t(log(pseudo / (1 - pseudo)))                                      # groups x samples
  cond <- factor(samples$condition)
  n <- table(cond)
  if (nlevels(cond) < 2 || any(n < 2)) {
    stop("Differential abundance needs at least 2 samples in each of at least 2 conditions.")
  }
  mean_by <- vapply(levels(cond), function(lv) colMeans(props[cond == lv, , drop = FALSE]),
                    numeric(ncol(props)))
  if (nlevels(cond) == 2) {
    fit <- abundance_ebayes(limma::lmFit(y, stats::model.matrix(~ cond)))
    tt <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "none")
    out <- data.frame(group = rownames(tt), mean_by[rownames(tt), , drop = FALSE],
                      prop_ratio = mean_by[rownames(tt), 2] / mean_by[rownames(tt), 1],
                      stat = tt$t, p = tt$P.Value, fdr = tt$adj.P.Val,
                      stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    fit <- abundance_ebayes(limma::lmFit(y, stats::model.matrix(~ cond))[, -1])
    tt <- limma::topTable(fit, number = Inf, sort.by = "none")
    out <- data.frame(group = rownames(tt), mean_by[rownames(tt), , drop = FALSE],
                      stat = tt$F, p = tt$P.Value, fdr = tt$adj.P.Val,
                      stringsAsFactors = FALSE, check.names = FALSE)
  }
  names(out)[seq_len(nlevels(cond)) + 1] <- paste0("prop_", levels(cond))
  attr(out, "robust") <- isTRUE(attr(fit, "robust"))
  out <- out[order(out$p), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "n") <- n
  attr(out, "props") <- props
  attr(out, "condition") <- cond
  out
}

# propeller's robust empirical Bayes; limma cannot fit its robust prior when the
# variances are (nearly) all zero -- proportions that barely vary between
# samples -- so the ordinary moderated statistics are used then (and reported)
abundance_ebayes <- function(fit) {
  out <- tryCatch(limma::eBayes(fit, robust = TRUE), error = function(e) NULL)
  if (!is.null(out)) return(structure(out, robust = TRUE))
  structure(limma::eBayes(fit, robust = FALSE), robust = FALSE)
}

#' Runnable R code for a pseudobulk DE run (Seurat object `obj`)
#' @keywords internal
pseudobulk_log_code <- function(sample_col, group_col, condition_col, ref, alt, method,
                                min_cells, covariate = NULL, groups = NULL) {
  design <- if (is.null(covariate)) "~ cond" else "~ cov + cond"
  fit <- switch(method,
    edgeR = c("    y <- edgeR::DGEList(y)",
              "    y <- y[edgeR::filterByExpr(y, design = design), , keep.lib.sizes = FALSE]",
              "    y <- edgeR::calcNormFactors(y)",
              "    y <- edgeR::estimateDisp(y, design)",
              "    fit <- edgeR::glmQLFit(y, design, robust = TRUE)",
              "    edgeR::topTags(edgeR::glmQLFTest(fit, coef = ncol(design)), n = Inf)$table"),
    limma = c("    y <- edgeR::DGEList(y)",
              "    y <- y[edgeR::filterByExpr(y, design = design), , keep.lib.sizes = FALSE]",
              "    y <- edgeR::calcNormFactors(y)",
              "    fit <- limma::eBayes(limma::lmFit(limma::voom(y, design), design), robust = TRUE)",
              "    limma::topTable(fit, coef = ncol(design), number = Inf)"),
    DESeq2 = c("    cd <- data.frame(cond = cond, row.names = colnames(y))",
               if (!is.null(covariate)) "    cd$cov <- cov",
               "    keep <- rowSums(y >= 10) >= min(table(cond))",
               sprintf("    dds <- DESeq2::DESeqDataSetFromMatrix(round(y[keep, ]), cd, %s)", design),
               "    dds <- DESeq2::DESeq(dds)",
               sprintf("    as.data.frame(DESeq2::results(dds, contrast = c(\"cond\", %s, %s)))", r_lit(alt), r_lit(ref))))
  c("md <- obj[[]]",
    "counts <- SeuratObject::LayerData(obj, assay = \"RNA\", layer = \"counts\")",
    sprintf("smp <- as.character(md[[%s]])", r_lit(sample_col)),
    sprintf("grp <- as.character(md[[%s]])", r_lit(group_col)),
    sprintf("cond_of <- tapply(as.character(md[[%s]]), smp, `[`, 1)   # one condition per sample", r_lit(condition_col)),
    if (!is.null(covariate)) sprintf("cov_of <- tapply(as.character(md[[%s]]), smp, `[`, 1)", r_lit(covariate)),
    sprintf("cell_types <- %s", if (is.null(groups)) "sort(unique(grp))" else r_lit(groups)),
    "res <- lapply(setNames(cell_types, cell_types), function(ct) {",
    "  cells <- grp == ct",
    "  f <- factor(smp[cells])",
    "  y <- counts[, cells] %*% Matrix::t(Matrix::fac2sparse(f))   # genes x samples (summed counts)",
    sprintf("  y <- as.matrix(y[, as.vector(table(f)) >= %s, drop = FALSE])   # drop samples with too few cells", min_cells),
    sprintf("  y <- y[, cond_of[colnames(y)] %%in%% %s, drop = FALSE]", r_lit(c(ref, alt))),
    sprintf("  cond <- factor(cond_of[colnames(y)], levels = %s)", r_lit(c(ref, alt))),
    if (!is.null(covariate)) "  cov <- factor(cov_of[colnames(y)])",
    "  if (any(table(cond) < 2)) return(NULL)   # too few samples for this cell type",
    sprintf("  design <- model.matrix(%s)", design),
    "  tryCatch({",
    fit,
    "  }, error = function(e) NULL)",
    "})",
    "res <- Filter(Negate(is.null), res)   # one table per cell type; FDR within each")
}

#' Runnable R code for a propeller differential-abundance run
#' @keywords internal
abundance_log_code <- function(sample_col, group_col, condition_col, n_levels = 2) {
  c("md <- obj[[]]",
    sprintf("tab <- table(md[[%s]], md[[%s]])   # samples x cell types", r_lit(sample_col), r_lit(group_col)),
    "pseudo <- (tab + 0.5) / rowSums(tab + 0.5)",
    "y <- t(log(pseudo / (1 - pseudo)))          # logit proportions",
    sprintf("cond <- factor(tapply(as.character(md[[%s]]), md[[%s]], `[`, 1)[rownames(tab)])",
            r_lit(condition_col), r_lit(sample_col)),
    "fit <- limma::lmFit(y, model.matrix(~ cond))",
    if (n_levels == 2) {
      c("fit <- tryCatch(limma::eBayes(fit, robust = TRUE), error = function(e) limma::eBayes(fit))",
        "limma::topTable(fit, coef = 2, number = Inf)")
    } else {
      c("fit <- fit[, -1]   # F test across conditions",
        "fit <- tryCatch(limma::eBayes(fit, robust = TRUE), error = function(e) limma::eBayes(fit))",
        "limma::topTable(fit, number = Inf)")
    },
    sprintf("# same method: speckle::propeller(clusters = md[[%s]], sample = md[[%s]], group = md[[%s]], transform = \"logit\")",
            r_lit(group_col), r_lit(sample_col), r_lit(condition_col)))
}

#' Volcano plot of one cell type's pseudobulk DE table
#' @keywords internal
pseudobulk_volcano <- function(df, fdr = 0.05, title = "", n_label = 12) {
  df$sig <- factor(ifelse(!is.na(df$fdr) & df$fdr < fdr, ifelse(df$logFC > 0, "up", "down"), "ns"),
                   levels = c("down", "ns", "up"))
  df$y <- -log10(pmax(df$p, 1e-300))
  lab <- utils::head(df[df$sig != "ns", , drop = FALSE], n_label)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$logFC, y = .data$y)) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$sig), size = 1.2, alpha = 0.75) +
    ggplot2::scale_colour_manual(values = c(down = style_tokens()$down, ns = style_tokens()$ns, up = style_tokens()$up),
                                 labels = c(down = sprintf("down, FDR < %g", fdr), ns = "n.s.",
                                            up = sprintf("up, FDR < %g", fdr)),
                                 drop = FALSE, name = NULL) +
    ggplot2::geom_vline(xintercept = 0, colour = style_tokens()$faint, linewidth = 0.3) +
    ggplot2::labs(x = "log2 fold change", y = expression(-log[10] ~ italic(p)), title = title) +
    omicone_theme()
  if (nrow(lab)) {
    p <- p + ggplot2::geom_text(data = lab, ggplot2::aes(label = .data$gene), size = 3, vjust = -0.7,
                                check_overlap = TRUE, show.legend = FALSE)
  }
  p
}

#' DE genes per cell type (bars up / down)
#' @keywords internal
pseudobulk_overview_plot <- function(overview, fdr = 0.05, ref = "", alt = "") {
  d <- rbind(data.frame(cell_type = overview$cell_type, dir = "up", n = overview$up),
             data.frame(cell_type = overview$cell_type, dir = "down", n = -overview$down))
  d$cell_type <- factor(d$cell_type, levels = rev(overview$cell_type[order(-(overview$up + overview$down))]))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$n, y = .data$cell_type, fill = .data$dir)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_vline(xintercept = 0, colour = style_tokens()$muted, linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = c(up = style_tokens()$up, down = style_tokens()$down),
                               labels = c(up = sprintf("higher in %s", alt), down = sprintf("lower in %s", alt)),
                               name = NULL) +
    ggplot2::scale_x_continuous(labels = abs) +
    ggplot2::labs(x = sprintf("genes with FDR < %g (%s vs %s)", fdr, alt, ref), y = NULL) +
    omicone_theme()
}

#' Per-sample proportions of each cell type, by condition, with the FDR
#' @keywords internal
abundance_plot <- function(res, cond_colors = NULL) {
  props <- attr(res, "props")
  cond <- attr(res, "condition")
  cond_colors <- cond_colors %||% value_colors(cond, levels = levels(cond))
  d <- data.frame(sample = rep(rownames(props), ncol(props)),
                  group = rep(colnames(props), each = nrow(props)),
                  prop = as.vector(props), condition = rep(as.character(cond), ncol(props)))
  d$group <- factor(d$group, levels = res$group)
  d$condition <- factor(d$condition, levels = levels(cond))
  lab <- data.frame(group = factor(res$group, levels = res$group),
                    txt = sprintf("FDR = %s", formatC(res$fdr, format = "g", digits = 2)))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$condition, y = .data$prop)) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.55, colour = style_tokens()$faint) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$condition),
                        position = ggplot2::position_jitter(width = 0.12, height = 0, seed = 1), size = 1.9) +
    ggplot2::geom_text(data = lab, ggplot2::aes(x = (nlevels(cond) + 1) / 2, y = Inf, label = .data$txt),
                       inherit.aes = FALSE, vjust = 1.4, size = 3, colour = style_tokens()$muted) +
    ggplot2::facet_wrap(~ group, scales = "free_y") +
    scale_group(cond_colors[levels(cond)], "colour", guide = "none") +
    ggplot2::scale_y_continuous(labels = function(v) paste0(round(100 * v), "%"),
                                expand = ggplot2::expansion(mult = c(0.05, 0.22))) +
    ggplot2::labs(x = NULL, y = "share of the sample's cells") +
    omicone_theme()
}

#' Stacked composition bar per sample, samples ordered within condition
#' @keywords internal
abundance_stack_plot <- function(res, group_colors = NULL) {
  props <- attr(res, "props")
  cond <- attr(res, "condition")
  gcols <- group_colors %||% value_colors(colnames(props))
  d <- data.frame(sample = rep(rownames(props), ncol(props)),
                  group = rep(colnames(props), each = nrow(props)),
                  prop = as.vector(props), condition = rep(as.character(cond), ncol(props)))
  d$condition <- factor(d$condition, levels = levels(cond))
  d$group <- factor(d$group, levels = names(gcols))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$sample, y = .data$prop, fill = .data$group)) +
    ggplot2::geom_col(width = 0.85) +
    ggplot2::facet_grid(~ condition, scales = "free_x", space = "free_x") +
    scale_group(gcols, "fill", name = NULL) +
    ggplot2::scale_y_continuous(labels = function(v) paste0(round(100 * v), "%"), expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = "share of the sample's cells") +
    omicone_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}
