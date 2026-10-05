#' Compute wrappers around Seurat / Bioconductor
#'
#' Each function performs one analysis step on a Seurat object and returns the
#' updated object (or a result). Kept separate from the UI so the science is
#' testable and the modules stay thin. All heavy dependencies are checked by the
#' caller via [require_pkgs()].
#'
#' Next to every compute function sit (a) the builder of its reproducibility
#' code -- native Seurat / harmony / scDblFinder / SingleR calls with the
#' parameters actually used, never OmicOne internals, so the exported script
#' runs without this package -- and (b) the pure data preparation and figure
#' for its preview, so both are unit-testable.
#'
#' @name fct_compute
#' @keywords internal
NULL

# ---- Code-generation and input helpers -------------------------------------

#' An R literal for a value, for reproducibility-log code
#' @param x An atomic value or vector.
#' @return A single string that parses back to `x`.
#' @keywords internal
r_lit <- function(x) {
  # no "L" suffixes: 30 reads better than 30L and parses to the same number
  ctrl <- c("keepNA", "niceNames", "showAttributes")
  paste(deparse(x, width.cutoff = 500L, control = ctrl), collapse = "")
}

#' A numeric input as a whole number, or NA
#'
#' `numericInput()` hands over whatever was typed: blank (NA), 30.5, -3. A
#' `%d` format or `seq_len()` on such a value errors out of the session, so
#' modules coerce first and `req(!is.na())` the result.
#' @param x Input value. @param lo,hi Allowed range (inclusive).
#' @return An integer, or `NA_integer_` when blank / non-numeric / out of range.
#' @keywords internal
int_input <- function(x, lo = -Inf, hi = Inf) {
  v <- suppressWarnings(as.numeric(x %||% NA))
  if (length(v) != 1 || !is.finite(v)) return(NA_integer_)
  v <- as.integer(round(v))
  if (v < lo || v > hi) return(NA_integer_)
  v
}

#' A numeric input as a finite number, or NA
#' @param x Input value. @param lo,hi Allowed range (inclusive).
#' @keywords internal
num_input <- function(x, lo = -Inf, hi = Inf) {
  v <- suppressWarnings(as.numeric(x %||% NA))
  if (length(v) != 1 || !is.finite(v) || v < lo || v > hi) return(NA_real_)
  v
}

#' Evenly spaced row indices, at most `n_max` of them (no RNG involved)
#' @param n Number of rows. @param n_max Cap.
#' @keywords internal
thin_index <- function(n, n_max = 20000) {
  if (n <= n_max) return(seq_len(n))
  unique(round(seq(1, n, length.out = n_max)))
}

# ---- QC ---------------------------------------------------------------------

#' Gene-name patterns for mitochondrial / ribosomal / haemoglobin genes
#'
#' The haemoglobin pattern lists the globin chains explicitly: the older
#' `^HB[^P]` also matched HBEGF, HBS1L and HBCBP, none of which is a globin.
#' @param species "human" or "mouse".
#' @return list(mt, ribo, hb) of regular expressions.
#' @keywords internal
qc_gene_patterns <- function(species = c("human", "mouse")) {
  species <- match.arg(species)
  if (species == "human") {
    return(list(mt = "^MT-", ribo = "^RP[SL]", hb = "^HB[ABDEGMQZ][0-9]*$"))
  }
  list(mt = "^mt-", ribo = "^Rp[sl]", hb = "^Hb[abdegmqz]([0-9]|-|$)")
}

#' Compute QC metrics (mitochondrial / ribosomal / hemoglobin percentages)
#' @param obj Seurat object.
#' @param species "human" or "mouse" (controls gene-name patterns).
#' @return Seurat object with `percent.mt`, `percent.ribo`, `percent.hb` (and
#'   `percent.diss` when dissociation genes are present) in meta.
#' @keywords internal
qc_add_metrics <- function(obj, species = c("human", "mouse")) {
  species <- match.arg(species)
  p <- qc_gene_patterns(species)
  obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = p$mt, assay = "RNA")
  obj[["percent.ribo"]] <- Seurat::PercentageFeatureSet(obj, pattern = p$ribo, assay = "RNA")
  obj[["percent.hb"]] <- Seurat::PercentageFeatureSet(obj, pattern = p$hb, assay = "RNA")
  # Dissociation/stress-gene percentage (borrowed from scCancer): flags cells
  # stressed during tissue dissociation. Uses genes present in the object.
  diss <- intersect(dissociation_genes(species), rownames(obj))
  if (length(diss) > 0) {
    obj[["percent.diss"]] <- Seurat::PercentageFeatureSet(obj, features = diss,
                                                          assay = "RNA")
  }
  obj
}

#' Are the feature names gene symbols or Ensembl IDs?
#'
#' With Ensembl IDs as row names no gene-name pattern matches, so the
#' mitochondrial percentage silently reads 0% for every cell.
#' @param genes Character vector of feature names.
#' @return "ensembl" or "symbol".
#' @keywords internal
feature_id_type <- function(genes) {
  g <- genes[!is.na(genes)]
  if (!length(g)) return("symbol")
  if (mean(grepl("^ENS[A-Z]*G[0-9]{6,}", g)) > 0.5) "ensembl" else "symbol"
}

#' Guess the species from gene names
#'
#' Human symbols are upper case (`MT-CO1`, `RPS6`), mouse symbols are title case
#' (`mt-Co1`, `Rps6`); Ensembl IDs carry the species in their prefix
#' (ENSG / ENSMUSG). Used so displays that have no species control (the import
#' overview) do not silently report 0% mitochondrial reads on mouse data.
#'
#' @param obj A Seurat object (or anything with gene symbols as rownames).
#' @return "human" or "mouse"; defaults to "human" when undecidable.
#' @keywords internal
guess_species <- function(obj) {
  g <- tryCatch(rownames(obj), error = function(e) NULL)
  if (!length(g)) return("human")
  if (identical(feature_id_type(g), "ensembl")) {
    return(if (mean(grepl("^ENSMUSG", g)) > 0.5) "mouse" else "human")
  }
  if (any(grepl("^mt-", g)) || any(grepl("^Rp[sl]", g))) return("mouse")
  if (any(grepl("^MT-", g)) || any(grepl("^RP[SL]", g))) return("human")
  # Fall back on overall casing: mouse symbols are mostly title case.
  alpha <- g[grepl("^[A-Za-z]{2,}", g)]
  if (length(alpha) && mean(grepl("^[A-Z][a-z]", alpha)) > 0.5) return("mouse")
  "human"
}

#' Dissociation/stress gene set (van den Brink et al.), human or mouse
#' @keywords internal
dissociation_genes <- function(species = c("human", "mouse")) {
  species <- match.arg(species)
  g <- c("FOS", "FOSB", "JUN", "JUNB", "JUND", "EGR1", "ATF3", "HSPA1A",
         "HSPA1B", "HSP90AB1", "HSPB1", "DNAJB1", "DNAJA1", "DUSP1", "IER2",
         "NR4A1", "PPP1R15A", "SOCS3", "ZFP36", "UBC", "HSPA8", "JUN")
  if (species == "mouse") g <- paste0(substr(g, 1, 1),
                                      tolower(substring(g, 2)))
  unique(g)
}

#' Median +/- n MADs of a vector, or (-Inf, Inf) when the MAD is zero
#'
#' A zero MAD happens when most values are identical -- typically
#' `percent.mt` on nuclei or on data whose mitochondrial genes are absent,
#' where more than half the cells read 0%. A median +/- n x 0 band would then
#' flag every cell with any mitochondrial read, so the metric is not
#' thresholded adaptively (the optional hard cap still applies).
#' @param x Numeric vector. @param nmads Multiplier.
#' @return Numeric length 2 with attribute `zero` (TRUE when the MAD was 0/NA).
#' @keywords internal
mad_bounds <- function(x, nmads) {
  med <- stats::median(x, na.rm = TRUE)
  dev <- stats::mad(x, center = med, na.rm = TRUE)
  if (!is.finite(dev) || dev <= 0) return(structure(c(-Inf, Inf), zero = TRUE))
  structure(c(med - nmads * dev, med + nmads * dev), zero = FALSE)
}

#' Adaptive (MAD-based) QC thresholds, globally or per sample (pure)
#'
#' Cells more than `nmads_lib` MADs from the median of log1p(nCount_RNA) or
#' log1p(nFeature_RNA), or more than `nmads_mt` MADs above the median
#' percent.mt, are outliers. With `batch`, medians and MADs are computed
#' within each sample, because depth and mitochondrial content differ between
#' samples and a pooled threshold over-filters the shallow ones. Thresholds
#' are returned on the original scale, rounded to 4 decimals, so the log can
#' print them and reproduce the exact same `keep` vector.
#'
#' @param md data.frame with `nCount_RNA`, `nFeature_RNA`, optional `percent.mt`.
#' @param nmads_lib,nmads_mt MAD multipliers.
#' @param batch Optional metadata column; NULL = all cells together.
#' @param max_mt Optional hard upper cap on percent.mt (NA = none).
#' @return data.frame, one row per sample: batch, n, nCount_lo, nCount_hi,
#'   nFeature_lo, nFeature_hi, mt_hi, mad_zero (metrics whose MAD was 0).
#' @keywords internal
qc_mad_thresholds <- function(md, nmads_lib = 5, nmads_mt = 3, batch = NULL,
                              max_mt = NA) {
  rows <- seq_len(nrow(md))
  groups <- if (is.null(batch)) list(rows) else split(rows, as.character(md[[batch]]))
  has_mt <- !is.null(md$percent.mt)
  out <- lapply(seq_along(groups), function(k) {
    i <- groups[[k]]
    cnt <- mad_bounds(log1p(md$nCount_RNA[i]), nmads_lib)
    fea <- mad_bounds(log1p(md$nFeature_RNA[i]), nmads_lib)
    mt <- structure(c(-Inf, Inf), zero = FALSE)
    if (has_mt) mt <- mad_bounds(md$percent.mt[i], nmads_mt)
    lo_up <- function(b) {
      lo <- if (is.finite(b[1])) round(expm1(b[1]), 4) else -Inf
      hi <- if (is.finite(b[2])) round(expm1(b[2]), 4) else Inf
      c(lo, hi)
    }
    c_b <- lo_up(cnt)
    f_b <- lo_up(fea)
    mt_hi <- if (is.finite(mt[2])) round(mt[2], 4) else Inf
    if (length(max_mt) == 1 && is.finite(max_mt)) mt_hi <- min(mt_hi, max_mt)
    zero <- c(nCount_RNA = attr(cnt, "zero"), nFeature_RNA = attr(fea, "zero"),
              percent.mt = has_mt && attr(mt, "zero"))
    data.frame(batch = if (is.null(batch)) NA_character_ else names(groups)[k],
               n = length(i),
               nCount_lo = c_b[1], nCount_hi = c_b[2],
               nFeature_lo = f_b[1], nFeature_hi = f_b[2],
               mt_hi = mt_hi,
               mad_zero = paste(names(zero)[zero], collapse = ","),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

#' Fixed (manual) QC thresholds in the same shape as [qc_mad_thresholds()]
#' @param min_genes,max_genes Genes-per-cell bounds. @param max_mt Mito cap (%).
#' @keywords internal
qc_manual_thresholds <- function(min_genes, max_genes, max_mt) {
  data.frame(batch = NA_character_, n = NA_integer_,
             nCount_lo = -Inf, nCount_hi = Inf,
             nFeature_lo = min_genes, nFeature_hi = max_genes,
             mt_hi = if (length(max_mt) == 1 && is.finite(max_mt)) max_mt else Inf,
             mad_zero = "", stringsAsFactors = FALSE)
}

#' Apply a QC threshold table to cell metadata (pure)
#'
#' Cells whose sample is not in the table, or with a missing metric, are not
#' kept: a missing value is never a pass.
#' @param md Cell metadata. @param thr Threshold table.
#' @param batch Column matching `thr$batch`, or NULL for a single global row.
#' @return Logical keep vector with the table as attribute `thresholds`.
#' @keywords internal
qc_apply_thresholds <- function(md, thr, batch = NULL) {
  n <- nrow(md)
  idx <- if (is.null(batch)) rep(1L, n) else match(as.character(md[[batch]]), thr$batch)
  keep <- md$nCount_RNA >= thr$nCount_lo[idx] & md$nCount_RNA <= thr$nCount_hi[idx] &
    md$nFeature_RNA >= thr$nFeature_lo[idx] & md$nFeature_RNA <= thr$nFeature_hi[idx]
  if (!is.null(md$percent.mt)) keep <- keep & md$percent.mt <= thr$mt_hi[idx]
  keep[is.na(keep)] <- FALSE
  attr(keep, "thresholds") <- thr
  keep
}

#' Adaptive (MAD-based) outlier flags from a metadata data.frame (pure)
#'
#' The testable core of the QC rule; see [qc_mad_thresholds()].
#' @param md data.frame with `nCount_RNA`, `nFeature_RNA`, optional `percent.mt`.
#' @param nmads_lib,nmads_mt MAD multipliers.
#' @param batch Optional sample column for per-sample thresholds.
#' @param max_mt Optional hard cap on percent.mt.
#' @return Logical vector: TRUE = keep, FALSE = flagged outlier (attribute
#'   `thresholds` holds the table used).
#' @keywords internal
qc_mad_keep_from_meta <- function(md, nmads_lib = 5, nmads_mt = 3, batch = NULL,
                                  max_mt = NA) {
  thr <- qc_mad_thresholds(md, nmads_lib, nmads_mt, batch = batch, max_mt = max_mt)
  qc_apply_thresholds(md, thr, batch)
}

#' Adaptive (MAD-based) outlier flags for a Seurat object
#' @param obj Seurat object with QC metrics.
#' @inheritParams qc_mad_keep_from_meta
#' @return Logical keep vector.
#' @keywords internal
qc_mad_keep <- function(obj, nmads_lib = 5, nmads_mt = 3, batch = NULL, max_mt = NA) {
  qc_mad_keep_from_meta(obj_meta(obj), nmads_lib, nmads_mt, batch = batch,
                        max_mt = max_mt)
}

#' Apply manual QC thresholds to a metadata data.frame (pure)
#' @param md data.frame with `nFeature_RNA`, optional `percent.mt`.
#' @keywords internal
qc_manual_keep_from_meta <- function(md, min_genes, max_genes, max_mt) {
  qc_apply_thresholds(md, qc_manual_thresholds(min_genes, max_genes, max_mt))
}

#' Apply manual QC thresholds to a Seurat object
#' @keywords internal
qc_manual_keep <- function(obj, min_genes, max_genes, max_mt) {
  qc_manual_keep_from_meta(obj_meta(obj), min_genes, max_genes, max_mt)
}

#' Reproducible R code for the QC step
#'
#' Recomputes the metrics, then applies the threshold table *as numbers*, so
#' the script keeps exactly the cells the app kept.
#' @param species "human"/"mouse". @param thr Threshold table used.
#' @param batch Sample column (NULL = global). @param rule One-line comment
#'   describing how the thresholds were derived.
#' @return Character vector of R code lines.
#' @keywords internal
qc_log_code <- function(species, thr, batch = NULL, rule = "") {
  p <- qc_gene_patterns(species)
  c(
    sprintf('obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = %s, assay = "RNA")',
            r_lit(p$mt)),
    sprintf('obj[["percent.ribo"]] <- Seurat::PercentageFeatureSet(obj, pattern = %s, assay = "RNA")',
            r_lit(p$ribo)),
    sprintf('obj[["percent.hb"]] <- Seurat::PercentageFeatureSet(obj, pattern = %s, assay = "RNA")',
            r_lit(p$hb)),
    sprintf("diss_genes <- intersect(%s, rownames(obj))", r_lit(dissociation_genes(species))),
    paste0('if (length(diss_genes)) obj[["percent.diss"]] <- ',
           'Seurat::PercentageFeatureSet(obj, features = diss_genes, assay = "RNA")'),
    if (nzchar(rule)) paste0("# ", rule) else NULL,
    sprintf(paste0("qc_thr <- data.frame(batch = %s, nCount_lo = %s, nCount_hi = %s,",
                   " nFeature_lo = %s, nFeature_hi = %s, mt_hi = %s)"),
            r_lit(thr$batch), r_lit(thr$nCount_lo), r_lit(thr$nCount_hi),
            r_lit(thr$nFeature_lo), r_lit(thr$nFeature_hi), r_lit(thr$mt_hi)),
    "md <- obj[[]]",
    if (is.null(batch)) "i <- rep(1L, nrow(md))"
    else sprintf("i <- match(as.character(md[[%s]]), qc_thr$batch)", r_lit(batch)),
    "keep <- md$nCount_RNA >= qc_thr$nCount_lo[i] & md$nCount_RNA <= qc_thr$nCount_hi[i] &",
    "  md$nFeature_RNA >= qc_thr$nFeature_lo[i] & md$nFeature_RNA <= qc_thr$nFeature_hi[i]",
    "if (!is.null(md$percent.mt)) keep <- keep & md$percent.mt <= qc_thr$mt_hi[i]",
    "keep[is.na(keep)] <- FALSE",
    "obj <- subset(obj, cells = colnames(obj)[keep])"
  )
}

#' QC preview: per-cell metrics with the thresholds that were applied
#'
#' Three panels (genes, UMIs, mitochondrial %) of violins with every cell as a
#' point, coloured kept / flagged, and the threshold of each sample as a
#' dashed line. With a sample column, each sample is its own violin.
#' @param md Pre-filter cell metadata. @param keep Logical keep vector.
#' @param thr Threshold table. @param batch Sample column or NULL.
#' @return A ggplot.
#' @keywords internal
qc_plot <- function(md, keep, thr, batch = NULL) {
  grp <- if (is.null(batch)) rep("all cells", nrow(md)) else as.character(md[[batch]])
  status <- factor(ifelse(as.vector(keep), "kept", "flagged"), levels = c("kept", "flagged"))
  spec <- list(
    list(col = "nFeature_RNA", lab = "Genes per cell (log10)", log = TRUE,
         lo = "nFeature_lo", hi = "nFeature_hi"),
    list(col = "nCount_RNA", lab = "UMIs per cell (log10)", log = TRUE,
         lo = "nCount_lo", hi = "nCount_hi"),
    list(col = "percent.mt", lab = "Mitochondrial reads (%)", log = FALSE,
         lo = NULL, hi = "mt_hi"))
  spec <- Filter(function(s) !is.null(md[[s$col]]), spec)
  tf <- function(v, lg) {
    if (!lg) return(v)
    out <- rep(NA_real_, length(v))
    ok <- !is.na(v) & v > 0
    out[ok] <- log10(v[ok])
    out
  }
  idx <- thin_index(nrow(md), 30000)
  pts <- do.call(rbind, lapply(spec, function(s) {
    data.frame(metric = s$lab, group = grp[idx], status = status[idx],
               value = tf(md[[s$col]][idx], s$log), stringsAsFactors = FALSE)
  }))
  pts <- pts[order(pts$status), , drop = FALSE]
  vio <- do.call(rbind, lapply(spec, function(s) {
    data.frame(metric = s$lab, group = grp, value = tf(md[[s$col]], s$log),
               stringsAsFactors = FALSE)
  }))
  thr_grp <- if (is.null(batch)) rep("all cells", nrow(thr)) else thr$batch
  lines <- do.call(rbind, lapply(spec, function(s) {
    v <- c(if (!is.null(s$lo)) thr[[s$lo]], thr[[s$hi]])
    g <- c(if (!is.null(s$lo)) thr_grp, thr_grp)
    data.frame(metric = s$lab, group = g, y = tf(ifelse(is.finite(v), v, NA), s$log),
               stringsAsFactors = FALSE)
  }))
  lines <- lines[!is.na(lines$y), , drop = FALSE]
  lv <- vapply(spec, function(s) s$lab, character(1))
  pts$metric <- factor(pts$metric, levels = lv)
  vio$metric <- factor(vio$metric, levels = lv)
  if (nrow(lines)) lines$metric <- factor(lines$metric, levels = lv)
  p <- ggplot2::ggplot(vio, ggplot2::aes(x = .data$group, y = .data$value)) +
    ggplot2::geom_violin(fill = "grey93", colour = "grey65", scale = "width",
                         trim = TRUE, na.rm = TRUE) +
    ggplot2::geom_point(data = pts, ggplot2::aes(colour = .data$status),
                        position = ggplot2::position_jitter(width = 0.25, height = 0, seed = 1),
                        size = 0.4, alpha = 0.55, na.rm = TRUE)
  if (nrow(lines)) {
    p <- p + ggplot2::geom_errorbar(data = lines,
                                    ggplot2::aes(x = .data$group, ymin = .data$y, ymax = .data$y),
                                    inherit.aes = FALSE, width = 0.85, linetype = "dashed",
                                    colour = style_tokens()$removed, linewidth = 0.6)
  }
  n_grp <- length(unique(grp))
  p <- p +
    ggplot2::facet_wrap(~metric, scales = "free_y", nrow = 1) +
    ggplot2::scale_colour_manual(values = c(kept = style_tokens()$kept, flagged = style_tokens()$removed),
                                 drop = FALSE, name = NULL) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
    ggplot2::labs(x = if (is.null(batch)) NULL else batch, y = NULL,
                  title = sprintf("QC: %s of %s cells kept",
                                  format(sum(keep), big.mark = ","),
                                  format(length(keep), big.mark = ",")),
                  subtitle = "Dashed lines: thresholds applied (per sample when a sample column is set)") +
    omicone_theme()
  if (n_grp > 6) {
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }
  p
}

# ---- Doublets ---------------------------------------------------------------

#' Expected doublet fraction under scDblFinder's default prior
#'
#' scDblFinder assumes 0.8% doublets per 1,000 recovered cells *per sample*
#' (`dbr.per1k = 0.008`), so 10,000 cells in one sample -> 8%.
#' @param n_cells Cells in the sample.
#' @keywords internal
doublet_expected_rate <- function(n_cells) {
  0.008 * n_cells / 1000
}

#' Final doublet call: the method's own class, or a custom score threshold
#' @param score Doublet scores. @param cls The method's classes.
#' @param threshold NULL / NA for the method's call, else a score cut-off.
#' @return Character vector "singlet" / "doublet".
#' @keywords internal
doublet_call <- function(score, cls, threshold = NULL) {
  if (length(threshold) == 1 && is.finite(threshold)) {
    return(ifelse(score > threshold, "doublet", "singlet"))
  }
  as.character(cls)
}

#' Run scDblFinder per sample and add the call to the metadata
#'
#' scDblFinder's prior grows with the number of cells it is given, and its
#' artificial doublets pair any two cells it sees; run on several samples
#' pooled, it would over-estimate the doublet rate and simulate cross-sample
#' doublets that cannot physically occur. `samples` makes it work within each
#' sample. A minimal SingleCellExperiment is built from the RNA counts (no
#' `as.SingleCellExperiment()`, which drags every assay and layer along).
#'
#' @param obj Seurat object.
#' @param samples Metadata column naming the sample of each cell, or NULL for
#'   a single-sample dataset.
#' @param seed Random seed (artificial-doublet simulation is stochastic).
#' @param threshold Optional score cut-off overriding scDblFinder's call.
#' @return Seurat object with `doublet_score` and `doublet_class` in meta.
#' @keywords internal
run_doublets <- function(obj, samples = NULL, seed = 42, threshold = NULL) {
  counts <- obj_counts(obj, "RNA")
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
  smp <- NULL
  if (!is.null(samples)) {
    smp <- as.character(obj_meta(obj)[[samples]])
    if (!length(smp) || anyNA(smp)) {
      stop("Sample column '", samples, "' is missing or has empty values.")
    }
  }
  set.seed(seed)
  sce <- if (is.null(smp)) scDblFinder::scDblFinder(sce)
         else scDblFinder::scDblFinder(sce, samples = smp)
  cd <- SummarizedExperiment::colData(sce)
  obj$doublet_score <- cd$scDblFinder.score
  obj$doublet_class <- doublet_call(cd$scDblFinder.score, cd$scDblFinder.class, threshold)
  obj
}

#' Reproducible R code for the doublet step
#' @param samples Sample column or NULL. @param seed Seed used.
#' @param threshold Custom cut-off or NULL. @param action "flag" / "remove".
#' @keywords internal
doublet_log_code <- function(samples = NULL, seed = 42, threshold = NULL,
                             action = "remove") {
  c(
    'counts <- SeuratObject::LayerData(obj, assay = "RNA", layer = "counts")',
    "sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))",
    sprintf("set.seed(%s)", r_lit(seed)),
    if (is.null(samples)) "sce <- scDblFinder::scDblFinder(sce)"
    else sprintf("sce <- scDblFinder::scDblFinder(sce, samples = as.character(obj[[]][[%s]]))",
                 r_lit(samples)),
    "obj$doublet_score <- sce$scDblFinder.score",
    "obj$doublet_class <- as.character(sce$scDblFinder.class)",
    if (length(threshold) == 1 && is.finite(threshold)) {
      sprintf('obj$doublet_class <- ifelse(obj$doublet_score > %s, "doublet", "singlet")',
              r_lit(threshold))
    },
    if (identical(action, "remove")) {
      'obj <- subset(obj, cells = colnames(obj)[obj$doublet_class == "singlet"])'
    }
  )
}

#' Tidy data for the doublet preview (pure)
#'
#' scDblFinder returns its class as a factor; `ifelse()` on a factor yields its
#' integer codes, which once lost the singlet/doublet colouring. The class is
#' therefore taken through `as.character()` into a factor with fixed levels.
#' @param md Metadata with `doublet_score`, `doublet_class`.
#' @param sample_col Optional sample column (facets).
#' @return data.frame(score, class, sample).
#' @keywords internal
doublet_plot_data <- function(md, sample_col = NULL) {
  smp <- if (!is.null(sample_col) && sample_col %in% names(md)) {
    as.character(md[[sample_col]])
  } else {
    rep("all cells", nrow(md))
  }
  data.frame(score = as.numeric(md$doublet_score),
             class = factor(as.character(md$doublet_class),
                            levels = c("singlet", "doublet")),
             sample = smp, stringsAsFactors = FALSE)
}

#' Doublet preview: score histogram coloured by the final call
#' @param pd Output of [doublet_plot_data()].
#' @param threshold Custom cut-off (drawn as a vertical line) or NULL.
#' @keywords internal
doublet_plot <- function(pd, threshold = NULL) {
  n_d <- sum(pd$class == "doublet", na.rm = TRUE)
  p <- ggplot2::ggplot(pd, ggplot2::aes(x = .data$score, fill = .data$class)) +
    ggplot2::geom_histogram(bins = 50, alpha = 0.8, position = "identity", na.rm = TRUE) +
    ggplot2::scale_fill_manual(values = c(singlet = style_tokens()$kept, doublet = style_tokens()$removed),
                               drop = FALSE, na.value = "grey60", name = NULL) +
    ggplot2::labs(x = "Doublet score (scDblFinder)", y = "Cells",
                  title = sprintf("Doublet scores: %s of %s cells called doublets",
                                  format(n_d, big.mark = ","),
                                  format(nrow(pd), big.mark = ","))) +
    omicone_theme()
  if (length(threshold) == 1 && is.finite(threshold)) {
    p <- p + ggplot2::geom_vline(xintercept = threshold, linetype = "dashed",
                                 colour = "grey25")
  }
  if (length(unique(pd$sample)) > 1) {
    p <- p + ggplot2::facet_wrap(~sample, scales = "free_y")
  }
  p
}

# ---- Normalization ----------------------------------------------------------

#' Normalize counts (LogNormalize or SCTransform), always from the RNA counts
#'
#' The default assay is set to RNA first: once SCTransform has run, the
#' default is "SCT", and a later LogNormalize would otherwise normalise the
#' SCT assay. Switching back to LogNormalize drops the SCT assay, so later
#' steps cannot pick up stale Pearson residuals.
#' @param obj Seurat object.
#' @param method "LogNormalize" or "SCT".
#' @param scale_factor LogNormalize target total.
#' @param vars_to_regress SCT only: metadata columns to regress out.
#' @keywords internal
normalize_obj <- function(obj, method = c("LogNormalize", "SCT"), scale_factor = 1e4,
                          vars_to_regress = NULL) {
  method <- match.arg(method)
  SeuratObject::DefaultAssay(obj) <- "RNA"
  if (method == "LogNormalize") {
    obj <- Seurat::NormalizeData(obj, assay = "RNA", normalization.method = "LogNormalize",
                                 scale.factor = scale_factor, verbose = FALSE)
    if ("SCT" %in% obj_assays(obj)) obj[["SCT"]] <- NULL
  } else {
    miss <- setdiff(vars_to_regress, obj_meta_cols(obj))
    if (length(miss)) stop("Columns to regress not found: ", paste(miss, collapse = ", "))
    vars <- if (length(vars_to_regress)) vars_to_regress else NULL
    obj <- Seurat::SCTransform(obj, assay = "RNA", new.assay.name = "SCT",
                               vars.to.regress = vars, verbose = FALSE)
  }
  obj
}

#' Reproducible R code for the normalization step
#' @param method,scale_factor,vars_to_regress As used.
#' @param dropped_sct Whether an SCT assay was removed (LogNormalize after SCT).
#' @keywords internal
normalize_log_code <- function(method, scale_factor = 1e4, vars_to_regress = NULL,
                               dropped_sct = FALSE) {
  c('SeuratObject::DefaultAssay(obj) <- "RNA"',
    if (method == "LogNormalize") {
      c(sprintf(paste0("obj <- Seurat::NormalizeData(obj, assay = \"RNA\", ",
                       "normalization.method = \"LogNormalize\", scale.factor = %s)"),
                r_lit(scale_factor)),
        if (isTRUE(dropped_sct)) 'obj[["SCT"]] <- NULL')
    } else {
      sprintf(paste0("obj <- Seurat::SCTransform(obj, assay = \"RNA\", new.assay.name = \"SCT\", ",
                     "vars.to.regress = %s)"),
              r_lit(if (length(vars_to_regress)) vars_to_regress else NULL))
    })
}

#' Depth diagnostic data: expression of the top genes vs library size
#'
#' For each cell: log10 library size against the mean expression of the 500
#' most-expressed genes, before (raw counts) and after normalisation
#' (`expm1()` of the active assay's `data`, i.e. normalised counts). Both are
#' on the linear scale: averaging log values would mix depth with detection
#' rate (shallow cells have more zeros), and show a depth trend that no
#' scaling normalisation can remove. Spearman rho per panel measures how much
#' depth still drives the values.
#' @param obj Normalised Seurat object. @param n_top Genes used.
#' @param max_cells Evenly thinned cells for speed.
#' @return data.frame(log10_lib, mean_expr, stage) with attribute `rho`.
#' @keywords internal
normalize_diag_data <- function(obj, n_top = 500, max_cells = 5000) {
  counts <- obj_counts(obj, "RNA")
  normed <- obj_layer(obj, "data")
  if (is.null(normed)) stop("No normalised data found.")
  tot <- Matrix::rowSums(counts)
  genes <- names(sort(tot[tot > 0], decreasing = TRUE))
  genes <- utils::head(intersect(genes, rownames(normed)), n_top)
  cells <- colnames(counts)[thin_index(ncol(counts), max_cells)]
  lib <- log10(Matrix::colSums(counts[, cells, drop = FALSE]) + 1)
  before <- Matrix::colMeans(counts[genes, cells, drop = FALSE])
  after <- Matrix::colMeans(expm1(normed[genes, cells, drop = FALSE]))
  stage_lv <- c("Before: raw counts", "After: normalised counts")
  dd <- data.frame(log10_lib = c(lib, lib), mean_expr = c(before, after),
                   stage = factor(rep(stage_lv, each = length(cells)), levels = stage_lv))
  rho <- c(before = suppressWarnings(stats::cor(lib, before, method = "spearman")),
           after = suppressWarnings(stats::cor(lib, after, method = "spearman")))
  attr(dd, "rho") <- rho
  attr(dd, "n_genes") <- length(genes)
  dd
}

#' Normalization preview: depth effect before vs after
#' @param dd Output of [normalize_diag_data()]. @param method Method label.
#' @keywords internal
normalize_diag_plot <- function(dd, method = "") {
  rho <- attr(dd, "rho")
  lab <- data.frame(stage = factor(levels(dd$stage), levels = levels(dd$stage)),
                    txt = sprintf("Spearman rho = %.2f", rho))
  dd <- dd[dd$mean_expr > 0, , drop = FALSE]
  ggplot2::ggplot(dd, ggplot2::aes(x = .data$log10_lib, y = .data$mean_expr)) +
    ggplot2::geom_point(size = 0.5, alpha = 0.4, colour = style_tokens()$kept) +
    ggplot2::geom_text(data = lab, ggplot2::aes(label = .data$txt), x = -Inf, y = Inf,
                       hjust = -0.08, vjust = 1.6, size = 4, inherit.aes = FALSE) +
    ggplot2::facet_wrap(~stage, scales = "free_y") +
    ggplot2::scale_y_log10() +
    ggplot2::labs(x = "Library size (log10 UMIs)",
                  y = sprintf("Mean of top %d genes (counts, log10 axis)", attr(dd, "n_genes")),
                  title = sprintf("Depth effect before vs after %s", method)) +
    omicone_theme()
}

# ---- Features + PCA ---------------------------------------------------------

#' HVG selection, scaling, and PCA
#'
#' Every existing reduction is removed first: Harmony, integrated spaces and
#' 2-D maps were computed from the old PCA and would otherwise stay
#' selectable. When SCTransform is the active assay, its variable genes and
#' Pearson residuals (already in `scale.data`) are used as they are --
#' re-running FindVariableFeatures / ScaleData would overwrite the residuals.
#'
#' @param obj Normalised Seurat object.
#' @param n_hvg Number of HVGs (vst / dispersion; ignored by mvp and SCT).
#' @param npcs PCs to compute (capped at features - 1 and cells - 1).
#' @param hvg_method "vst", "mvp" or "dispersion".
#' @param seed PCA seed.
#' @keywords internal
reduce_obj <- function(obj, n_hvg = 2000, npcs = 50, hvg_method = "vst", seed = 42) {
  obj <- drop_reductions(obj, obj_reductions(obj))
  if (identical(obj_default_assay(obj), "SCT")) {
    features <- SeuratObject::VariableFeatures(obj)
    if (!length(features)) stop("The SCT assay has no variable features; re-run Normalize.")
  } else {
    obj <- if (hvg_method == "mvp") {
      Seurat::FindVariableFeatures(obj, selection.method = "mvp", verbose = FALSE)
    } else {
      Seurat::FindVariableFeatures(obj, selection.method = hvg_method,
                                   nfeatures = n_hvg, verbose = FALSE)
    }
    features <- SeuratObject::VariableFeatures(obj)
    obj <- Seurat::ScaleData(obj, features = features, verbose = FALSE)
  }
  npcs <- reduce_npcs(npcs, length(features), ncol(obj))
  obj <- Seurat::RunPCA(obj, features = features, npcs = npcs, seed.use = seed,
                        verbose = FALSE)
  obj@misc$omicone_reduction <- "pca"
  obj
}

#' PCs that can actually be computed
#' @param npcs Requested. @param n_features,n_cells Matrix size.
#' @keywords internal
reduce_npcs <- function(npcs, n_features, n_cells) {
  as.integer(max(2, min(npcs, n_features - 1, n_cells - 1)))
}

#' Reproducible R code for the features / PCA step
#' @param sct Whether the SCT branch ran. @param hvg_method,n_hvg,npcs,seed As
#'   used. @param dropped Reductions removed before PCA.
#' @keywords internal
reduce_log_code <- function(sct, hvg_method, n_hvg, npcs, seed = 42, dropped = character(0)) {
  drop <- setdiff(dropped, "pca")
  c(if (length(drop)) sprintf("for (r in %s) obj[[r]] <- NULL  # built on the old PCA", r_lit(drop)),
    if (sct) {
      "# SCT is the active assay: keep its variable genes and Pearson residuals (no re-scaling)"
    } else if (hvg_method == "mvp") {
      'obj <- Seurat::FindVariableFeatures(obj, selection.method = "mvp")'
    } else {
      sprintf("obj <- Seurat::FindVariableFeatures(obj, selection.method = %s, nfeatures = %s)",
              r_lit(hvg_method), r_lit(as.integer(n_hvg)))
    },
    if (!sct) "obj <- Seurat::ScaleData(obj, features = SeuratObject::VariableFeatures(obj))",
    sprintf(paste0("obj <- Seurat::RunPCA(obj, features = SeuratObject::VariableFeatures(obj), ",
                   "npcs = %s, seed.use = %s)"), r_lit(as.integer(npcs)), r_lit(seed)))
}

#' Variance explained per PC
#'
#' Seurat stores the total variance of the scaled matrix in
#' `obj[["pca"]]@misc$total.variance`; percentages are of that total. When it
#' is absent the share is relative to the PCs computed, and `basis` says so.
#' @param obj Seurat object with a "pca" reduction.
#' @return data.frame(PC, stdev, pct, cum, basis).
#' @keywords internal
pca_variance <- function(obj) {
  red <- tryCatch(obj[["pca"]], error = function(e) NULL)
  if (is.null(red)) stop("No PCA found.")
  sd <- red@stdev
  tv <- tryCatch(red@misc$total.variance, error = function(e) NULL)
  total <- is.numeric(tv) && length(tv) == 1 && is.finite(tv) && tv > 0
  pct <- if (total) 100 * sd^2 / tv else 100 * sd^2 / sum(sd^2)
  data.frame(PC = seq_along(sd), stdev = sd, pct = pct, cum = cumsum(pct),
             basis = if (total) "total" else "computed", stringsAsFactors = FALSE)
}

#' HVG table (mean vs ranking statistic) for the preview
#'
#' Column names differ by method: vst `variance.standardized`, mvp
#' `dispersion.scaled`, dispersion `dispersion`, SCT `residual_variance` (read
#' from the SCT model's feature attributes; `HVFInfo()` refuses SCT assays).
#' @param obj Seurat object after [reduce_obj()]. @param method HVG method or
#'   "SCT".
#' @return data.frame(gene, mean, score, variable) with attributes `xlab`,
#'   `ylab`, or NULL when unavailable.
#' @keywords internal
hvg_plot_data <- function(obj, method = "vst") {
  info <- NULL
  if (identical(method, "SCT")) {
    info <- tryCatch(Seurat::SCTResults(obj[["SCT"]], slot = "feature.attributes"),
                     error = function(e) NULL)
    if (is.list(info) && !is.data.frame(info)) info <- info[[1]]
  } else {
    # SeuratObject 5 names the argument `method`, version 4 `selection.method`
    f <- tryCatch(names(formals(utils::getS3method("HVFInfo", "Seurat",
                                                   envir = asNamespace("SeuratObject")))),
                  error = function(e) character(0))
    arg <- if ("method" %in% f) "method" else "selection.method"
    info <- tryCatch(do.call(SeuratObject::HVFInfo, stats::setNames(list(obj, method), c("", arg))),
                     error = function(e) NULL)
  }
  if (is.null(info)) info <- tryCatch(SeuratObject::HVFInfo(obj), error = function(e) NULL)
  if (is.null(info) || !nrow(info)) return(NULL)
  cols <- colnames(info)
  pick <- function(pat) {
    hit <- cols[endsWith(cols, pat)]
    if (length(hit)) hit[1] else NULL
  }
  ycol <- switch(method,
                 SCT = pick("residual_variance"),
                 mvp = pick("dispersion.scaled"),
                 dispersion = pick("dispersion"),
                 pick("variance.standardized"))
  ycol <- ycol %||% cols[ncol(info)]
  xcol <- pick("gmean") %||% pick("mean") %||% cols[1]
  vf <- SeuratObject::VariableFeatures(obj)
  df <- data.frame(gene = rownames(info), mean = info[[xcol]], score = info[[ycol]],
                   variable = rownames(info) %in% vf, stringsAsFactors = FALSE)
  attr(df, "xlab") <- xcol
  attr(df, "ylab") <- ycol
  attr(df, "logx") <- method %in% c("vst", "SCT") && all(df$mean > 0, na.rm = TRUE)
  df
}

#' Features / PCA preview: HVG plot + variance explained per PC
#' @param hvg Output of [hvg_plot_data()] (or NULL). @param pv [pca_variance()].
#' @keywords internal
reduce_plot <- function(hvg, pv) {
  ylab <- if (identical(pv$basis[1], "total")) "Variance explained (% of total)"
          else "Variance explained (% of computed PCs)"
  p_elb <- ggplot2::ggplot(pv, ggplot2::aes(x = .data$PC, y = .data$pct)) +
    ggplot2::geom_line(colour = style_tokens()$faint, linewidth = 0.4) +
    ggplot2::geom_point(colour = style_tokens()$kept, size = 1.4) +
    ggplot2::labs(x = "Principal component", y = ylab,
                  title = sprintf("PCA: first 10 PCs explain %.1f%%",
                                  pv$cum[min(10, nrow(pv))])) +
    omicone_theme()
  if (is.null(hvg)) return(p_elb)
  hvg <- hvg[order(hvg$variable), , drop = FALSE]
  top <- utils::head(hvg[hvg$variable, , drop = FALSE][order(-hvg$score[hvg$variable]), ], 10)
  p_hvg <- ggplot2::ggplot(hvg, ggplot2::aes(x = .data$mean, y = .data$score,
                                             colour = .data$variable)) +
    ggplot2::geom_point(size = 0.6, alpha = 0.6, na.rm = TRUE) +
    ggplot2::geom_text(data = top, ggplot2::aes(label = .data$gene), size = 3,
                       vjust = -0.6, check_overlap = TRUE, colour = style_tokens()$ink,
                       show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = style_tokens()$removed, `FALSE` = style_tokens()$ns),
                                 labels = c(`TRUE` = "variable", `FALSE` = "other"),
                                 name = NULL) +
    ggplot2::labs(x = attr(hvg, "xlab"), y = attr(hvg, "ylab"),
                  title = sprintf("%s variable genes", format(sum(hvg$variable), big.mark = ","))) +
    omicone_theme() +
    ggplot2::theme(legend.position = "bottom")
  if (isTRUE(attr(hvg, "logx"))) p_hvg <- p_hvg + ggplot2::scale_x_log10()
  compose_grid(list(p_hvg, p_elb + ggplot2::theme(legend.position = "bottom")), nrow = 1, ncol = 2,
               widths = c(1.2, 1))
}

# ---- Integration ------------------------------------------------------------

#' Validate a batch column before integration
#' @param md Cell metadata. @param batch Column name.
#' @keywords internal
check_batch_col <- function(md, batch) {
  if (is.null(batch) || !length(batch) || !batch %in% names(md)) {
    stop("Batch column '", batch %||% "", "' not found.")
  }
  v <- md[[batch]]
  if (anyNA(v)) {
    stop(sprintf("Batch column '%s' has %d missing values; fill or remove them first.",
                 batch, sum(is.na(v))))
  }
  if (length(unique(v)) < 2) stop("Batch column '", batch, "' has a single level.")
  invisible(TRUE)
}

#' The argument harmony uses for the input reduction ("reduction.use" from 1.0)
#' @keywords internal
harmony_reduction_arg <- function() {
  f <- tryCatch(names(formals(utils::getS3method("RunHarmony", "Seurat",
                                                  envir = asNamespace("harmony")))),
                error = function(e) character(0))
  if ("reduction" %in% f && !"reduction.use" %in% f) "reduction" else "reduction.use"
}

#' Batch integration (Harmony by default; CCA/RPCA via Seurat v5)
#'
#' Any earlier integrated reduction is removed first, so choosing "none" (or
#' switching method) never leaves an old Harmony space selectable. The space
#' downstream steps should use is recorded in `obj@misc$omicone_reduction`.
#' CCA / RPCA reuse the variable genes and PCA from the Features/PCA step
#' (they are not recomputed with defaults) and refuse an SCT assay, which
#' needs Seurat's per-sample SCT workflow.
#'
#' @param obj Seurat object (normalized, with PCA).
#' @param batch Metadata column identifying batch/sample.
#' @param method "harmony" (default), "none", "CCA", or "RPCA".
#' @param dims PCs used (capped at the PCs available).
#' @param seed Random seed.
#' @return Seurat object with an integrated reduction ("harmony" or
#'   "integrated.dr").
#' @keywords internal
integrate_obj <- function(obj, batch, method = c("harmony", "none", "CCA", "RPCA"),
                          dims = 30, seed = 42) {
  method <- match.arg(method)
  obj <- drop_reductions(obj, c("harmony", "integrated.dr"))
  obj@misc$omicone_reduction <- "pca"
  if (method == "none") return(obj)
  if (!has_reduction(obj, "pca")) stop("Run Features / PCA before integrating.")
  check_batch_col(obj_meta(obj), batch)
  dims <- min(dims, ncol(SeuratObject::Embeddings(obj, reduction = "pca")))

  if (method == "harmony") {
    args <- list(object = obj, group.by.vars = batch, dims.use = seq_len(dims),
                 verbose = FALSE)
    args[[harmony_reduction_arg()]] <- "pca"
    set.seed(seed)
    obj <- do.call(harmony::RunHarmony, args)
    obj@misc$omicone_reduction <- "harmony"
    return(obj)
  }

  if (!"IntegrateLayers" %in% getNamespaceExports("Seurat")) {
    stop("CCA/RPCA integration requires Seurat v5 (IntegrateLayers). ",
         "Use Harmony, or upgrade Seurat.")
  }
  if (identical(obj_default_assay(obj), "SCT")) {
    stop("CCA/RPCA after SCTransform needs Seurat's per-sample SCT workflow, which ",
         "OmicOne does not run. Use Harmony, or normalise with LogNormalize.")
  }
  hvg <- SeuratObject::VariableFeatures(obj)
  if (!length(hvg)) stop("No variable features; run Features / PCA first.")
  b <- as.character(obj_meta(obj)[[batch]])
  if (!methods::is(obj[["RNA"]], "Assay5")) obj[["RNA"]] <- methods::as(obj[["RNA"]], "Assay5")
  obj[["RNA"]] <- split(obj[["RNA"]], f = b)
  k_weight <- integrate_k_weight(b)
  fun <- getExportedValue("Seurat", if (method == "CCA") "CCAIntegration" else "RPCAIntegration")
  set.seed(seed)
  obj <- Seurat::IntegrateLayers(obj, method = fun, orig.reduction = "pca",
                                 new.reduction = "integrated.dr", features = hvg,
                                 dims = seq_len(dims), k.weight = k_weight,
                                 verbose = FALSE)
  obj[["RNA"]] <- SeuratObject::JoinLayers(obj[["RNA"]])
  SeuratObject::VariableFeatures(obj) <- hvg
  obj@misc$omicone_reduction <- "integrated.dr"
  obj
}

#' Anchor weighting neighbours: Seurat's 100, or fewer for small samples
#' @param b Batch labels per cell.
#' @keywords internal
integrate_k_weight <- function(b) {
  as.integer(max(5, min(100, min(table(b)) - 1)))
}

#' Reproducible R code for the integration step
#' @param method,batch,dims,seed As used. @param k_weight CCA/RPCA only.
#' @keywords internal
integrate_log_code <- function(method, batch = NULL, dims = 30, seed = 42, k_weight = 100) {
  drop <- 'for (r in intersect(c("harmony", "integrated.dr"), names(obj@reductions))) obj[[r]] <- NULL'
  if (method == "none") return(c(drop, "# no batch correction: downstream steps use \"pca\""))
  if (method == "harmony") {
    return(c(drop,
             sprintf("set.seed(%s)", r_lit(seed)),
             sprintf("obj <- harmony::RunHarmony(obj, group.by.vars = %s, %s = \"pca\", dims.use = %s)",
                     r_lit(batch), harmony_reduction_arg(), r_lit(seq_len(dims)))))
  }
  c(drop,
    "hvg <- SeuratObject::VariableFeatures(obj)",
    'if (!methods::is(obj[["RNA"]], "Assay5")) obj[["RNA"]] <- methods::as(obj[["RNA"]], "Assay5")',
    sprintf('obj[["RNA"]] <- split(obj[["RNA"]], f = as.character(obj[[]][[%s]]))', r_lit(batch)),
    sprintf("set.seed(%s)", r_lit(seed)),
    sprintf(paste0("obj <- Seurat::IntegrateLayers(obj, method = Seurat::%s, orig.reduction = \"pca\", ",
                   "new.reduction = \"integrated.dr\", features = hvg, dims = %s, k.weight = %s)"),
            if (method == "CCA") "CCAIntegration" else "RPCAIntegration",
            r_lit(seq_len(dims)), r_lit(k_weight)),
    'obj[["RNA"]] <- SeuratObject::JoinLayers(obj[["RNA"]])',
    "SeuratObject::VariableFeatures(obj) <- hvg")
}

#' Batch-mixing entropy of a kNN graph, normalised to the batch composition
#'
#' For each (thinned) cell, the Shannon entropy of the batch labels among its
#' k nearest neighbours, averaged and divided by the entropy of the global
#' batch composition: 1 = neighbourhoods as mixed as the composition allows,
#' 0 = every neighbourhood from a single batch. High mixing alone does not
#' show that biology was preserved (over-correction also scores high).
#' @param emb Cells x dims embedding. @param batch Batch label per cell.
#' @param k Neighbours. @param max_cells Query cells (evenly thinned).
#' @return A number in [0, ~1], or NA.
#' @keywords internal
batch_mixing_entropy <- function(emb, batch, k = 30, max_cells = 5000) {
  batch <- as.character(batch)
  n <- nrow(emb)
  lv <- unique(batch)
  if (n < 3 || length(lv) < 2) return(NA_real_)
  k <- min(k, n - 1)
  q <- thin_index(n, max_cells)
  nn <- RANN::nn2(data = emb, query = emb[q, , drop = FALSE], k = k + 1)$nn.idx
  nn <- nn[, -1, drop = FALSE]
  ent <- function(lab) {
    p <- tabulate(match(lab, lv), nbins = length(lv)) / length(lab)
    p <- p[p > 0]
    -sum(p * log(p))
  }
  h_glob <- ent(batch)
  h <- apply(nn, 1, function(r) ent(batch[r]))
  mean(h) / h_glob
}

#' Integration preview data: uncorrected PCA next to the integrated space
#' @param obj Seurat object. @param batch Batch column.
#' @param reduction Integrated reduction (or "pca" when none).
#' @param dims Dims used for the mixing metric.
#' @return list(df = plotting data, entropy = named numeric).
#' @keywords internal
integrate_plot_data <- function(obj, batch, reduction = "pca", dims = 30) {
  md <- obj_meta(obj)
  reds <- unique(c("pca", reduction))
  lab <- c(pca = "PCA (uncorrected)")
  lab[reduction] <- if (reduction == "pca") lab[["pca"]] else paste0(reduction, " (integrated)")
  idx <- thin_index(nrow(md), 30000)
  ent <- stats::setNames(rep(NA_real_, length(reds)), reds)
  parts <- lapply(reds, function(r) {
    emb <- SeuratObject::Embeddings(obj, reduction = r)
    d <- seq_len(min(dims, ncol(emb)))
    if (has_pkg("RANN")) ent[[r]] <<- batch_mixing_entropy(emb[, d, drop = FALSE], md[[batch]])
    data.frame(dim1 = emb[idx, 1], dim2 = emb[idx, 2], batch = as.character(md[[batch]][idx]),
               reduction = r, stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, parts)
  panel <- vapply(reds, function(r) {
    if (is.na(ent[[r]])) lab[[r]] else sprintf("%s: mixing %.2f", lab[[r]], ent[[r]])
  }, character(1))
  df$panel <- factor(panel[df$reduction], levels = panel)
  list(df = df, entropy = ent)
}

#' Integration preview: batches on the uncorrected and integrated spaces
#' @param pd Output of [integrate_plot_data()]. @param batch Batch column name.
#' @param colors Named batch colours (default: from the levels).
#' @keywords internal
integrate_plot <- function(pd, batch, colors = NULL) {
  df <- pd$df
  cols <- colors %||% value_colors(df$batch)
  df$batch <- factor(as.character(df$batch), levels = names(cols))
  set.seed(1)
  df <- df[sample.int(nrow(df)), , drop = FALSE]
  ggplot2::ggplot(df, ggplot2::aes(x = .data$dim1, y = .data$dim2, colour = .data$batch)) +
    ggplot2::geom_point(size = dim_point_size(nrow(df) / 2), alpha = 0.7, stroke = 0, shape = 16) +
    ggplot2::facet_wrap(~panel, scales = "free") +
    scale_group(cols, "colour", name = batch) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1))) +
    ggplot2::labs(title = "Batches before and after integration",
                  subtitle = "Mixing: kNN (k = 30) batch entropy / entropy of the batch composition (1 = fully mixed)") +
    omicone_dim_theme() +
    ggplot2::theme(aspect.ratio = NULL, strip.text = ggplot2::element_text(hjust = 0.5))
}

# ---- Clustering -------------------------------------------------------------

#' Which Leiden implementation the installed Seurat can use
#'
#' Seurat >= 5.2 runs Leiden through the R package leidenbase; versions that
#' expose `leiden_method` can fall back to igraph. Older Seurat calls the
#' Python leidenalg through the `leiden` package.
#' @return "leidenbase", "igraph", "leiden-python", or NA when none is usable.
#' @keywords internal
leiden_backend <- function() {
  if (!has_pkg("Seurat")) return(NA_character_)
  fm <- unlist(lapply(c("Seurat", "default"), function(cl) {
    tryCatch(names(formals(utils::getS3method("FindClusters", cl,
                                              envir = asNamespace("Seurat")))),
             error = function(e) character(0))
  }))
  if (utils::packageVersion("Seurat") >= "5.2.0") {
    if (has_pkg("leidenbase")) return("leidenbase")
    if ("leiden_method" %in% fm && has_pkg("igraph")) return("igraph")
    return(NA_character_)
  }
  if (has_pkg("leiden")) "leiden-python" else NA_character_
}

#' Extra FindClusters arguments for the chosen algorithm
#' @param algorithm 1 (Louvain) or 4 (Leiden).
#' @keywords internal
cluster_extra_args <- function(algorithm) {
  if (algorithm != 4) return(list())
  be <- leiden_backend()
  if (is.na(be)) {
    stop("Leiden (algorithm 4) needs the R package 'leidenbase' with Seurat >= 5.2 ",
         "(install.packages(\"leidenbase\")); or choose Louvain.")
  }
  if (be == "igraph") list(leiden_method = "igraph") else list()
}

#' Neighbors + clustering at one or more resolutions
#'
#' Clustering columns of earlier runs (`*_snn_res.*`) are removed first: they
#' were computed on a different graph. One column per resolution is written
#' (`<assay>_snn_res.<r>`); the first resolution becomes the active
#' clustering (see [cluster_set_active()]), and the user can switch.
#'
#' @param obj Seurat object.
#' @param reduction Reduction to build the graph on.
#' @param dims Leading dimensions (capped at those available).
#' @param resolutions Numeric vector.
#' @param algorithm 1 = Louvain, 4 = Leiden.
#' @param k Neighbours (`k.param`).
#' @param seed `random.seed` for FindClusters.
#' @return Seurat object; `obj@misc$omicone_cluster_cols` maps resolution ->
#'   column.
#' @keywords internal
cluster_obj <- function(obj, reduction = "pca", dims = 30, resolutions = 0.5,
                        algorithm = 4, k = 20, seed = 0) {
  extra <- cluster_extra_args(algorithm)
  dims <- min(dims, ncol(SeuratObject::Embeddings(obj, reduction = reduction)))
  assay <- obj_default_assay(obj) %||% "RNA"
  graphs <- paste0(assay, c("_nn", "_snn"))
  md <- obj_meta(obj)
  old <- grep("_snn_res\\.", names(md), value = TRUE)
  for (cl in old) obj[[cl]] <- NULL
  obj <- Seurat::FindNeighbors(obj, reduction = reduction, dims = seq_len(dims),
                               k.param = k, graph.name = graphs, verbose = FALSE)
  cols <- character(0)
  for (res in resolutions) {
    args <- c(list(obj, graph.name = graphs[2], resolution = res, algorithm = algorithm,
                   random.seed = seed, verbose = FALSE), extra)
    obj <- tryCatch(do.call(Seurat::FindClusters, args), error = function(e) {
      if (algorithm != 4) stop(e)
      stop("Leiden clustering failed (", conditionMessage(e), "). Seurat < 5.2 runs Leiden ",
           "through Python leidenalg; install Seurat >= 5.2 with 'leidenbase', or choose Louvain.",
           call. = FALSE)
    })
    cols[[as.character(res)]] <- paste0(graphs[2], "_res.", res)
  }
  obj@misc$omicone_cluster_cols <- cols
  for (cl in cols) obj <- set_group_colors(obj, cl, extend_colors(NULL, level_order(obj_meta(obj)[[cl]])))
  cluster_set_active(obj, cols[[1]])
}

#' Make one clustering column the active clustering
#'
#' Writes `Idents`, `seurat_clusters` (read by later steps) and
#' `obj@misc$omicone_cluster_col`.
#' @param obj Seurat object. @param col Metadata column.
#' @keywords internal
cluster_set_active <- function(obj, col) {
  if (!col %in% obj_meta_cols(obj)) stop("Clustering column '", col, "' not found.")
  SeuratObject::Idents(obj) <- col
  obj$seurat_clusters <- SeuratObject::Idents(obj)
  obj@misc$omicone_cluster_col <- col
  set_group_colors(obj, "seurat_clusters", group_colors(obj, col))
}

#' Reproducible R code for the clustering step
#' @param reduction,dims,resolutions,algorithm,k,seed As used.
#' @param assay Default assay (graph names). @param active Active column.
#' @param leiden_igraph Whether `leiden_method = "igraph"` was passed.
#' @keywords internal
cluster_log_code <- function(reduction, dims, resolutions, algorithm, k = 20, seed = 0,
                             assay = "RNA", active = NULL, leiden_igraph = FALSE) {
  g <- paste0(assay, c("_nn", "_snn"))
  c('obj@meta.data <- obj@meta.data[, !grepl("_snn_res\\\\.", colnames(obj@meta.data)), drop = FALSE]',
    sprintf("obj <- Seurat::FindNeighbors(obj, reduction = %s, dims = %s, k.param = %s, graph.name = %s)",
            r_lit(reduction), r_lit(seq_len(dims)), r_lit(as.integer(k)), r_lit(g)),
    sprintf("for (r in %s) {", r_lit(resolutions)),
    sprintf("  obj <- Seurat::FindClusters(obj, graph.name = %s, resolution = r, algorithm = %s, random.seed = %s%s)",
            r_lit(g[2]), r_lit(as.integer(algorithm)), r_lit(seed),
            if (isTRUE(leiden_igraph)) ', leiden_method = "igraph"' else ""),
    "}",
    if (!is.null(active)) c(sprintf("SeuratObject::Idents(obj) <- %s", r_lit(active)),
                            "obj$seurat_clusters <- SeuratObject::Idents(obj)"))
}

#' Clustering preview: clusters on a 2-D map, or cluster sizes
#' @param obj Seurat object. @param col Active clustering column.
#' @param label Resolution label for the title.
#' @keywords internal
cluster_plot <- function(obj, col, label = "") {
  md <- obj_meta(obj)
  cols <- group_colors(obj, col)
  red <- intersect(c("umap", "tsne"), obj_reductions(obj))[1]
  if (!is.na(red)) {
    return(sc_dim_plot(obj, red, col, colors = cols,
                       title = sprintf("%d clusters on the %s (resolution %s)", length(cols), toupper(red), label)))
  }
  counts <- as.data.frame(table(cluster = factor(as.character(md[[col]]), levels = names(cols))),
                          stringsAsFactors = FALSE)
  names(counts) <- c("cluster", "n")
  counts$cluster <- factor(counts$cluster, levels = names(cols))
  ggplot2::ggplot(counts, ggplot2::aes(x = .data$cluster, y = .data$n, fill = .data$cluster)) +
    ggplot2::geom_col(width = 0.75) +
    scale_group(cols, "fill", guide = "none") +
    ggplot2::labs(x = "Cluster", y = "Cells",
                  title = sprintf("Cluster sizes (resolution %s)", label),
                  caption = "No 2-D map yet: run the Embed step to see clusters on a UMAP.") +
    omicone_theme(grid = "y")
}

# ---- Embedding --------------------------------------------------------------

#' Embedding parameters the data can support
#'
#' Caps dims at the reduction's width, `n_neighbors` below the cell count and
#' t-SNE perplexity below (n - 1) / 3 (Rtsne's limit).
#' @param n_cells,n_dims Data size. @param dims,n_neighbors,perplexity Requested.
#' @return list(dims, n_neighbors, perplexity).
#' @keywords internal
embed_params <- function(n_cells, n_dims, dims, n_neighbors = 30, perplexity = 30) {
  list(dims = as.integer(max(2, min(dims, n_dims))),
       n_neighbors = as.integer(max(2, min(n_neighbors, n_cells - 1))),
       perplexity = max(1, min(perplexity, floor((n_cells - 1) / 3 - 1e-9))))
}

#' Non-linear embedding for visualization (UMAP / t-SNE)
#' @param obj Seurat object. @param method "umap" or "tsne".
#' @param reduction Input reduction. @param dims,n_neighbors,min_dist,perplexity
#'   Parameters (capped by [embed_params()]). @param seed Seed.
#' @keywords internal
embed_obj <- function(obj, method = c("umap", "tsne"), reduction = "pca",
                      dims = 30, n_neighbors = 30, min_dist = 0.3, perplexity = 30,
                      seed = NULL) {
  method <- match.arg(method)
  pr <- embed_params(ncol(obj), ncol(SeuratObject::Embeddings(obj, reduction = reduction)),
                     dims, n_neighbors, perplexity)
  d <- seq_len(pr$dims)
  if (method == "umap") {
    obj <- Seurat::RunUMAP(obj, reduction = reduction, dims = d,
                           n.neighbors = pr$n_neighbors, min.dist = min_dist,
                           seed.use = seed %||% 42, verbose = FALSE)
  } else {
    obj <- Seurat::RunTSNE(obj, reduction = reduction, dims = d,
                           perplexity = pr$perplexity, seed.use = seed %||% 1)
  }
  obj
}

#' Reproducible R code for the embedding step
#' @param method,reduction,min_dist As used. @param pr [embed_params()] output.
#' @keywords internal
embed_log_code <- function(method, reduction, pr, min_dist = 0.3) {
  if (method == "umap") {
    sprintf(paste0("obj <- Seurat::RunUMAP(obj, reduction = %s, dims = %s, n.neighbors = %s, ",
                   "min.dist = %s, seed.use = 42)"),
            r_lit(reduction), r_lit(seq_len(pr$dims)), r_lit(pr$n_neighbors), r_lit(min_dist))
  } else {
    sprintf("obj <- Seurat::RunTSNE(obj, reduction = %s, dims = %s, perplexity = %s, seed.use = 1)",
            r_lit(reduction), r_lit(seq_len(pr$dims)), r_lit(pr$perplexity))
  }
}

#' Embedding preview: the 2-D map coloured by a chosen grouping
#' @param obj Seurat object. @param reduction "umap" / "tsne".
#' @param color_by Metadata column or NULL.
#' @keywords internal
embed_plot <- function(obj, reduction, color_by = NULL) {
  sc_dim_plot(obj, reduction, color_by, label = !is.null(color_by) && length(group_colors(obj, color_by)) <= 30)
}

# ---- Markers ----------------------------------------------------------------

#' Marker genes per group (cell-level marker discovery)
#'
#' `group_by` is applied by setting `Idents` on a local copy, which works on
#' every Seurat version (FindAllMarkers has no `group.by` before v5).
#' @param obj Seurat object. @param group_by Metadata column, or NULL for the
#'   current Idents. @param test,logfc,min_pct,only_pos FindAllMarkers options.
#' @keywords internal
markers_obj <- function(obj, group_by = NULL, test = "wilcox", logfc = 0.25,
                        min_pct = 0.1, only_pos = TRUE) {
  if (!is.null(group_by)) SeuratObject::Idents(obj) <- group_by
  Seurat::FindAllMarkers(obj, test.use = test, logfc.threshold = logfc,
                         min.pct = min_pct, only.pos = only_pos, verbose = FALSE)
}

#' Reproducible R code for the markers step
#' @keywords internal
markers_log_code <- function(group_by = NULL, test = "wilcox", logfc = 0.25,
                             min_pct = 0.1, only_pos = TRUE) {
  call <- sprintf(paste0("markers <- Seurat::FindAllMarkers(%s, test.use = %s, logfc.threshold = %s, ",
                         "min.pct = %s, only.pos = %s)"),
                  if (is.null(group_by)) "obj" else "obj_m",
                  r_lit(test), r_lit(logfc), r_lit(min_pct), r_lit(isTRUE(only_pos)))
  if (is.null(group_by)) return(call)
  c("obj_m <- obj",
    sprintf("SeuratObject::Idents(obj_m) <- %s", r_lit(group_by)),
    call,
    "rm(obj_m)")
}

#' Top-N markers per group: significant first, then by fold change (pure)
#'
#' Genes with BH-adjusted p >= `padj` are dropped before ranking, so a large
#' but non-significant fold change never heads the list. The ROC test has no
#' p-values; there genes are ranked by `power` (or `myAUC`). `label` is
#' `gene___cluster`, unique even when several clusters share a gene.
#' @param df FindAllMarkers output. @param n Genes per group.
#' @param padj Adjusted-p cut-off.
#' @return data.frame with a `label` column (possibly zero rows).
#' @keywords internal
marker_top_n <- function(df, n = 5, padj = 0.05) {
  if (is.null(df) || !nrow(df)) return(df)
  if (is.null(df$gene)) df$gene <- rownames(df)
  if (!is.null(df$p_val_adj)) df <- df[!is.na(df$p_val_adj) & df$p_val_adj < padj, , drop = FALSE]
  rank_col <- if (!is.null(df$avg_log2FC)) "avg_log2FC"
              else if (!is.null(df$power)) "power"
              else if (!is.null(df$myAUC)) "myAUC" else NULL
  parts <- split(df, df$cluster, drop = TRUE)
  picked <- lapply(parts, function(d) {
    ord <- if (is.null(rank_col)) seq_len(nrow(d)) else order(-d[[rank_col]])
    utils::head(d[ord, , drop = FALSE], n)
  })
  out <- do.call(rbind, picked)
  if (is.null(out)) return(df[0, , drop = FALSE])
  rownames(out) <- NULL
  out$label <- paste0(out$gene, "___", out$cluster)
  out
}

#' Markers preview (bar chart): top-N fold changes per group
#' @param top Output of [marker_top_n()].
#' @keywords internal
markers_bar_plot <- function(top, colors = NULL) {
  cols <- colors %||% value_colors(top$cluster)
  top$cluster <- factor(as.character(top$cluster), levels = names(cols))
  y <- if (!is.null(top$avg_log2FC)) "avg_log2FC" else if (!is.null(top$power)) "power" else "myAUC"
  ggplot2::ggplot(top, ggplot2::aes(x = stats::reorder(.data$label, .data[[y]]),
                                    y = .data[[y]], fill = .data$cluster)) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~cluster, scales = "free_y") +
    ggplot2::scale_x_discrete(labels = function(x) sub("___.*$", "", x)) +
    scale_group(cols, "fill", guide = "none") +
    ggplot2::labs(x = NULL, y = if (y == "avg_log2FC") "Average log2 fold change" else y,
                  title = "Top marker genes per group (adjusted p < 0.05)") +
    omicone_theme(grid = "x") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(face = "italic"))
}

#' Markers preview (dot plot) of the top markers, grouped by the group they mark
#' @param obj Seurat object. @param top [marker_top_n()] output.
#' @param group_by Grouping column or NULL (Idents).
#' @keywords internal
markers_dot_plot <- function(obj, top, group_by = NULL) {
  group_by <- group_by %||% "seurat_clusters"
  lv <- names(group_colors(obj, group_by))
  top <- top[order(match(as.character(top$cluster), lv)), , drop = FALSE]
  top <- top[!duplicated(top$gene), , drop = FALSE]
  sc_dot_plot(obj, top$gene, group_by,
              gene_groups = stats::setNames(as.character(top$cluster), top$gene),
              title = "Top markers: dot size = % of cells expressing, colour = scaled mean expression")
}

# ---- Annotation -------------------------------------------------------------

#' Species of a celldex reference
#' @param ref celldex function name.
#' @keywords internal
singler_ref_species <- function(ref) {
  if (ref %in% c("ImmGenData", "MouseRNAseqData")) "mouse" else "human"
}

#' Margin between the best and second-best label score (SingleR delta)
#' @param pred SingleR result.
#' @keywords internal
singler_delta <- function(pred) {
  if ("delta.next" %in% colnames(pred)) return(as.numeric(pred$delta.next))
  sc <- as.matrix(pred$scores)
  apply(sc, 1, function(r) {
    s <- sort(r, decreasing = TRUE)
    if (length(s) > 1) s[1] - s[2] else NA_real_
  })
}

#' Reference-based annotation with SingleR, per cluster and per cell
#'
#' Clusters are labelled from their pooled profile (`clusters =`), which is
#' far more stable than a per-cell vote, and the pruned label is used: a
#' cluster whose best label does not clearly beat the rest gets
#' "Unassigned" instead of a guess. Per-cell pruned labels are kept in
#' `SingleR_cell`. The test matrix is the log-normalised data of the active
#' assay (SingleR correlates ranks, so no SingleCellExperiment is needed).
#'
#' @param obj Seurat object. @param ref A celldex reference (SummarizedExperiment).
#' @param labels Reference labels. @param clusters Cluster label per cell.
#' @return list(obj, table) -- table: cluster, label, raw_label, confidence
#'   (delta), low_conf, n.
#' @keywords internal
annotate_singler <- function(obj, ref, labels, clusters) {
  logm <- obj_layer(obj, "data")
  if (is.null(logm)) stop("No normalised data: run Normalize first.")
  cl <- factor(clusters)
  pred_cl <- SingleR::SingleR(test = logm, ref = ref, labels = labels, clusters = cl)
  pred_cell <- SingleR::SingleR(test = logm, ref = ref, labels = labels)
  lab <- as.character(pred_cl$pruned.labels)
  low <- is.na(lab)
  lab[low] <- "Unassigned"
  names(lab) <- rownames(pred_cl)
  obj$celltype <- unname(lab[as.character(cl)])
  obj$SingleR_cell <- as.character(pred_cell$pruned.labels)
  n <- table(cl)
  tab <- data.frame(cluster = rownames(pred_cl), label = unname(lab),
                    raw_label = as.character(pred_cl$labels),
                    confidence = singler_delta(pred_cl), low_conf = low,
                    n = as.integer(n[rownames(pred_cl)]), stringsAsFactors = FALSE)
  list(obj = obj, table = tab)
}

#' Reproducible R code for SingleR annotation
#' @param ref celldex function name. @param level "main" or "fine".
#' @keywords internal
annotate_singler_log_code <- function(ref, level = "main") {
  lab <- sprintf("ref$label.%s", level)
  c(sprintf("ref <- celldex::%s()", ref),
    'logm <- SeuratObject::LayerData(obj, layer = "data")',
    "clusters <- SeuratObject::Idents(obj)",
    sprintf("pred_cl <- SingleR::SingleR(test = logm, ref = ref, labels = %s, clusters = clusters)", lab),
    sprintf("pred_cell <- SingleR::SingleR(test = logm, ref = ref, labels = %s)", lab),
    "lab <- pred_cl$pruned.labels",
    'lab[is.na(lab)] <- "Unassigned"',
    "names(lab) <- rownames(pred_cl)",
    "obj$celltype <- unname(lab[as.character(clusters)])",
    "obj$SingleR_cell <- pred_cell$pruned.labels")
}

#' The Azimuth label column for an annotation level (pure)
#'
#' References name levels differently (`celltype.l2`, `annotation.l1`,
#' `ann_level_3`); a reference without numbered levels gets its first
#' predicted column.
#' @param cols Metadata column names. @param level "1", "2" or "3".
#' @return Column name, or NULL when Azimuth wrote none.
#' @keywords internal
azimuth_label_col <- function(cols, level = "2") {
  cand <- grep("^predicted\\.", cols, value = TRUE)
  cand <- cand[!grepl("\\.score$", cand)]
  if (!length(cand)) return(NULL)
  hit <- cand[grepl(paste0("(\\.l|_level_)", level, "$"), cand)]
  if (length(hit)) hit[1] else cand[1]
}

#' Per-cluster summary of a per-cell annotation (pure)
#'
#' Majority label of each cluster, the share of its cells carrying it, and
#' the mean per-cell score when one is given.
#' @param clusters Cluster per cell. @param labels Label per cell.
#' @param score Optional score per cell. @param low_score Mean score below
#'   which a cluster is flagged.
#' @keywords internal
annotation_cluster_table <- function(clusters, labels, score = NULL, low_score = 0.5) {
  cl <- factor(clusters)
  rows <- lapply(levels(cl), function(k) {
    i <- which(cl == k)
    tb <- sort(table(as.character(labels[i])), decreasing = TRUE)
    conf <- if (is.null(score)) NA_real_ else mean(score[i], na.rm = TRUE)
    data.frame(cluster = k, label = if (length(tb)) names(tb)[1] else NA_character_,
               raw_label = if (length(tb)) names(tb)[1] else NA_character_,
               confidence = conf,
               low_conf = !is.na(conf) && conf < low_score,
               n = length(i),
               share = if (length(tb)) as.numeric(tb[1]) / length(i) else NA_real_,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' Azimuth reference mapping, keeping the default assay
#' @param obj Seurat object. @param reference Azimuth reference name.
#' @param level "1"/"2"/"3". @param clusters Cluster per cell.
#' @return list(obj, table, col).
#' @keywords internal
annotate_azimuth <- function(obj, reference = "pbmcref", level = "2", clusters) {
  assay0 <- obj_default_assay(obj)
  o <- Azimuth::RunAzimuth(obj, reference = reference)
  if (!is.null(assay0) && assay0 %in% obj_assays(o)) SeuratObject::DefaultAssay(o) <- assay0
  md <- obj_meta(o)
  col <- azimuth_label_col(colnames(md), level)
  if (is.null(col)) stop("Azimuth returned no predicted label column.")
  o$celltype <- as.character(md[[col]])
  score <- md[[paste0(col, ".score")]]
  list(obj = o, table = annotation_cluster_table(clusters, o$celltype, score), col = col)
}

#' Reproducible R code for Azimuth annotation
#' @param reference Reference name. @param col Label column used.
#' @param assay Default assay restored after mapping.
#' @keywords internal
annotate_azimuth_log_code <- function(reference, col, assay = "RNA") {
  c(sprintf("obj <- Azimuth::RunAzimuth(obj, reference = %s)", r_lit(reference)),
    sprintf("SeuratObject::DefaultAssay(obj) <- %s", r_lit(assay)),
    sprintf("obj$celltype <- as.character(obj[[]][[%s]])", r_lit(col)))
}

#' Reproducible R code for manual annotation
#' @param labels Named character vector cluster -> label.
#' @keywords internal
annotate_manual_log_code <- function(labels) {
  c(sprintf("labels <- %s", r_lit(labels)),
    "obj$celltype <- unname(labels[as.character(SeuratObject::Idents(obj))])")
}

#' Annotation preview: labelled map plus per-cluster confidence
#' @param obj Seurat object with `celltype`. @param tab Per-cluster table.
#' @param method "manual", "singler" or "azimuth".
#' @keywords internal
annotate_plot <- function(obj, tab, method = "manual") {
  tk <- style_tokens()
  tab$name <- sprintf("%s: %s", tab$cluster, tab$label)
  tab$name <- factor(tab$name, levels = rev(tab$name))
  has_conf <- any(!is.na(tab$confidence))
  tab$flag <- factor(ifelse(tab$low_conf, "low confidence", "assigned"),
                     levels = c("assigned", "low confidence"))
  xlab <- switch(method,
                 singler = "SingleR delta (best minus next label score)",
                 azimuth = "Mean Azimuth prediction score",
                 "Cells")
  p_bar <- ggplot2::ggplot(tab, ggplot2::aes(y = .data$name,
                                             x = if (has_conf) .data$confidence else .data$n,
                                             fill = .data$flag)) +
    ggplot2::geom_col(width = 0.72) +
    ggplot2::scale_fill_manual(values = c(assigned = tk$kept, `low confidence` = tk$removed),
                               drop = FALSE, name = NULL,
                               guide = if (has_conf) "legend" else "none") +
    ggplot2::labs(x = if (has_conf) xlab else "Cells", y = NULL,
                  title = if (has_conf) "Per-cluster confidence" else "Cells per cluster") +
    omicone_theme(grid = "x") +
    ggplot2::theme(legend.position = "bottom")
  red <- intersect(c("umap", "tsne"), obj_reductions(obj))[1]
  if (is.na(red) || !"celltype" %in% obj_meta_cols(obj)) return(p_bar)
  p_map <- sc_dim_plot(obj, red, "celltype", title = sprintf("%d cell types", length(group_colors(obj, "celltype")))) +
    ggplot2::theme(legend.position = "none")
  compose_grid(list(p_map, p_bar), nrow = 1, ncol = 2, widths = c(1.35, 1))
}

# The scop-engine wrappers (enrichment, trajectory, velocity, communication,
# CNV, signature scoring) live in fct_scop.R.
