#' Pre-QC data overview (survey the dataset before filtering)
#'
#' Helpers for the Import step: a scCancer/Seurat-style survey of the raw data so
#' the user can eyeball quality and structure before any QC. Computation is
#' read-only (metrics are computed for display, never written back / filtered).
#'
#' @name fct_overview
#' @keywords internal
NULL

#' Summary statistics for the header tiles
#'
#' @param obj A Seurat object.
#' @return Named list: cells, genes, meta_cols, median_genes, median_umi.
#' @keywords internal
data_overview <- function(obj) {
  dims <- obj_dims(obj)
  md <- obj_meta(obj)
  mg <- if (!is.null(md$nFeature_RNA)) stats::median(md$nFeature_RNA, na.rm = TRUE) else NA
  mu <- if (!is.null(md$nCount_RNA))   stats::median(md$nCount_RNA,   na.rm = TRUE) else NA
  list(cells = dims$cells, genes = dims$genes, meta_cols = ncol(md),
       median_genes = mg, median_umi = mu)
}

#' A capped dense slice of the counts matrix for a table preview
#'
#' @param obj A Seurat object.
#' @param n_genes,n_cells Max rows/cols to show.
#' @return A data.frame (genes x cells) small enough to render in a table.
#' @keywords internal
counts_preview <- function(obj, n_genes = 50, n_cells = 100) {
  m <- obj_layer(obj, "counts", "RNA")
  if (is.null(m)) stop("Could not access the counts matrix.")
  g <- seq_len(min(n_genes, nrow(m)))
  c <- seq_len(min(n_cells, ncol(m)))
  as.data.frame(as.matrix(m[g, c, drop = FALSE]))
}

#' Pre-QC overview figure: QC violins + top expressed genes + count scatter
#'
#' The Import step draws this once, at import, and keeps the figure: drawn
#' from the live object it would silently show post-QC data under a
#' "pre-filter" title.
#'
#' @param obj A Seurat object.
#' @param species "human"/"mouse" for mito/ribo gene patterns. Defaults to
#'   auto-detection from gene-name casing (this view has no species control, and
#'   guessing wrong makes the mitochondrial panel read a flat 0%).
#' @param when Title suffix saying which object is shown ("at import",
#'   "current object").
#' @return A patchwork/ggplot object (falls back to a single ggplot if patchwork
#'   is unavailable).
#' @keywords internal
overview_plots <- function(obj, species = guess_species(obj), when = "at import") {
  # compute QC metrics for display only (no filtering, not written back)
  o <- tryCatch(qc_add_metrics(obj, species = species), error = function(e) obj)
  md <- obj_meta(o)

  # --- 1. QC metric violins (nCount / nFeature / percent.mt) ---
  metrics <- intersect(c("nCount_RNA", "nFeature_RNA", "percent.mt"), colnames(md))
  labs_map <- c(nCount_RNA = "UMIs / cell", nFeature_RNA = "Genes / cell",
                percent.mt = "Mitochondrial %")
  p_vln <- NULL
  if (length(metrics)) {
    long <- do.call(rbind, lapply(metrics, function(mt) {
      data.frame(metric = labs_map[[mt]], value = md[[mt]], stringsAsFactors = FALSE)
    }))
    p_vln <- ggplot2::ggplot(long, ggplot2::aes(x = .data$metric, y = .data$value,
                                                fill = .data$metric)) +
      ggplot2::geom_violin(scale = "width", trim = TRUE, alpha = 0.85) +
      ggplot2::facet_wrap(~metric, scales = "free", nrow = 1) +
      ggplot2::scale_fill_manual(values = sc_palette(length(metrics)), guide = "none") +
      ggplot2::labs(x = NULL, y = NULL,
                    title = sprintf("Per-cell QC metrics (%s, %s cells)", when,
                                    format(nrow(md), big.mark = ","))) +
      omicone_theme() +
      ggplot2::theme(axis.text.x = ggplot2::element_blank())
  }

  # --- 2. Top-20 highly expressed genes (fraction of total counts) ---
  m <- obj_layer(o, "counts", "RNA")
  p_top <- NULL
  if (!is.null(m)) {
    gene_tot <- Matrix::rowSums(m)
    frac <- sort(gene_tot / sum(gene_tot), decreasing = TRUE)[1:min(20, length(gene_tot))]
    dtop <- data.frame(gene = factor(names(frac), levels = rev(names(frac))),
                       frac = as.numeric(frac) * 100)
    p_top <- ggplot2::ggplot(dtop, ggplot2::aes(x = .data$gene, y = .data$frac)) +
      ggplot2::geom_col(fill = sc_palette(1)) +
      ggplot2::coord_flip() +
      ggplot2::labs(x = NULL, y = "% of total counts", title = "Top expressed genes") +
      omicone_theme()
  }

  # --- 3. nCount vs nFeature scatter ---
  p_sc <- NULL
  if (all(c("nCount_RNA", "nFeature_RNA") %in% colnames(md))) {
    p_sc <- ggplot2::ggplot(md, ggplot2::aes(x = .data$nCount_RNA, y = .data$nFeature_RNA)) +
      ggplot2::geom_point(size = 0.5, alpha = 0.4, colour = sc_palette(1)) +
      ggplot2::scale_x_log10() +
      ggplot2::labs(x = "UMIs / cell (log10)", y = "Genes / cell",
                    title = "Counts vs genes") +
      omicone_theme()
  }

  plots <- Filter(Negate(is.null), list(p_vln, p_top, p_sc))
  if (!length(plots)) stop("Nothing to plot: no QC metrics or counts available.")
  if (patchwork_ok() && length(plots) > 1) {
    patchwork::wrap_plots(plots, ncol = 1)
  } else {
    plots[[1]]
  }
}
