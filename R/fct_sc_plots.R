#' Single-cell figures in the shared style (ggplot2, no Seurat plotting)
#'
#' The five figures every single-cell analysis shows -- the cell map, gene
#' expression on the map, violins, the dot plot and the heatmap -- drawn with
#' ggplot2 in the app's style (fct_style.R), with every grouping coloured by
#' [group_colors()], so a cluster or a cell type has the same colour in every
#' step. Seurat's DimPlot / FeaturePlot / VlnPlot / DotPlot / DoHeatmap are not
#' used.
#'
#' @name fct_sc_plots
#' @keywords internal
NULL

# the 2-D map to draw on: UMAP, else t-SNE, else the default reduction
sc_default_map <- function(obj) {
  red <- intersect(c("umap", "tsne"), obj_reductions(obj))[1]
  if (!is.na(red)) return(red)
  r <- tryCatch(SeuratObject::DefaultDimReduc(obj), error = function(e) NULL)
  if (is.null(r)) stop("No embedding yet: run the Embed step first.")
  r
}

# point size that keeps a map readable from 300 to 300,000 cells
dim_point_size <- function(n) max(0.15, min(2, 64 / sqrt(max(n, 1))))

#' Normalised expression of some genes, as a cells x genes data.frame
#' @param obj Seurat object. @param genes Genes (missing ones are dropped).
#' @return data.frame (rownames = cells); attribute "missing".
#' @keywords internal
sc_expr_df <- function(obj, genes) {
  md <- obj_meta(obj)
  meta <- genes[genes %in% names(md) & vapply(genes, function(g) g %in% names(md) && is.numeric(md[[g]]), logical(1))]
  m <- obj_layer(obj, "data")
  if (is.null(m) && !length(meta)) stop("No normalised expression: run the Normalize step first.")
  keep <- if (is.null(m)) character(0) else setdiff(genes[genes %in% rownames(m)], meta)
  if (!length(keep) && !length(meta)) stop("None of the genes is in the data: ", paste(genes, collapse = ", "))
  out <- if (length(keep)) as.data.frame(as.matrix(Matrix::t(m[keep, , drop = FALSE])), check.names = FALSE)
         else data.frame(row.names = rownames(md))
  for (g in meta) out[[g]] <- md[rownames(out), g]
  out <- out[, intersect(genes, names(out)), drop = FALSE]
  attr(out, "missing") <- setdiff(genes, names(out))
  attr(out, "meta") <- meta
  out
}

#' The cell map: an embedding coloured by a grouping (or a single colour)
#'
#' @param obj Seurat object. @param reduction "umap", "tsne", "pca", ...
#' @param group_by Metadata column, or NULL.
#' @param label Write each group's name at its centre.
#' @param colors Named colours (default [group_colors()]).
#' @param title Plot title (default: what is shown).
#' @param highlight Optional levels to emphasise (others drawn grey).
#' @keywords internal
sc_dim_plot <- function(obj, reduction = NULL, group_by = NULL, label = TRUE, colors = NULL,
                        title = NULL, highlight = NULL) {
  tk <- style_tokens()
  reduction <- reduction %||% sc_default_map(obj)
  df <- embedding_df(obj, reduction = reduction, color_by = group_by)
  n <- nrow(df)
  set.seed(1)
  df <- df[sample.int(n), , drop = FALSE]          # no group systematically drawn on top
  red <- toupper(reduction)
  xl <- paste0(red, " 1")
  yl <- paste0(red, " 2")
  pt <- dim_point_size(n)
  if (is.null(group_by) || is.null(df$color)) {
    p <- ggplot2::ggplot(df, ggplot2::aes(.data$dim1, .data$dim2)) +
      ggplot2::geom_point(size = pt, colour = tk$kept, alpha = 0.75, stroke = 0, shape = 16)
    ttl <- title %||% sprintf("%s of %s cells", red, format(n, big.mark = ","))
  } else {
    cols <- colors %||% group_colors(obj, group_by)
    df$color <- factor(as.character(df$color), levels = names(cols))
    if (!is.null(highlight)) {
      cols[!names(cols) %in% highlight] <- tk$ns
      df <- df[order(df$color %in% highlight), , drop = FALSE]
    }
    p <- ggplot2::ggplot(df, ggplot2::aes(.data$dim1, .data$dim2, colour = .data$color)) +
      ggplot2::geom_point(size = pt, alpha = 0.85, stroke = 0, shape = 16) +
      scale_group(cols, "colour", name = group_by) +
      ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1),
                                                     ncol = if (length(cols) > 16) 2 else 1))
    if (length(cols) > 40) p <- p + ggplot2::theme(legend.position = "none")
    if (isTRUE(label) && length(cols) <= 40) {
      cen <- stats::aggregate(cbind(dim1, dim2) ~ color, data = df, FUN = stats::median)
      p <- p + ggplot2::geom_label(data = cen, ggplot2::aes(label = .data$color), colour = tk$ink,
                                   fill = grDevices::adjustcolor("white", 0.78), label.size = 0,
                                   size = 3.2, fontface = "bold", label.padding = grid::unit(1.6, "pt"),
                                   show.legend = FALSE)
    }
    ttl <- title %||% sprintf("%s of %s cells, by %s", red, format(n, big.mark = ","), group_by)
  }
  p + dim_corner_axes(df$dim1, df$dim2, xl, yl) + ggplot2::labs(title = ttl) + omicone_dim_theme()
}

#' Gene expression on the cell map, one panel per gene
#'
#' Colour is the expression scaled per gene (0 = none, 1 = the gene's 99th
#' percentile), so a weakly expressed marker is as readable as a strong one;
#' expressing cells are drawn last (on top). Numeric metadata columns (module
#' scores, pseudotime) work too, scaled between their 1st and 99th percentile.
#' @param obj Seurat object. @param genes Genes or numeric metadata columns.
#' @param reduction Reduction (default: UMAP, else the default reduction).
#' @param ncol Panels per row (default: up to 3).
#' @keywords internal
sc_feature_plot <- function(obj, genes, reduction = NULL, ncol = NULL) {
  reduction <- reduction %||% sc_default_map(obj)
  ex <- sc_expr_df(obj, genes)
  emb <- embedding_df(obj, reduction = reduction)
  genes <- colnames(ex)
  meta <- attr(ex, "meta")
  long <- do.call(rbind, lapply(genes, function(g) {
    v <- ex[emb$cell, g]
    if (g %in% meta) {                                 # a score (may be negative): 1st-99th percentile
      q <- stats::quantile(v, c(0.01, 0.99), na.rm = TRUE, names = FALSE)
      sc <- if (diff(q) > 0) (v - q[1]) / diff(q) else v * 0
    } else {
      hi <- stats::quantile(v[v > 0], 0.99, names = FALSE)
      if (!is.finite(hi) || hi <= 0) hi <- max(v, 1e-9)
      sc <- v / hi
    }
    d <- data.frame(dim1 = emb$dim1, dim2 = emb$dim2, gene = g, value = pmax(0, pmin(sc, 1)))
    d[order(d$value), , drop = FALSE]
  }))
  long$gene <- factor(long$gene, levels = genes)
  red <- toupper(reduction)
  miss <- attr(ex, "missing")
  ggplot2::ggplot(long, ggplot2::aes(.data$dim1, .data$dim2, colour = .data$value)) +
    ggplot2::geom_point(size = dim_point_size(nrow(emb)) * 0.9, stroke = 0, shape = 16) +
    ggplot2::scale_colour_gradientn(colours = style_seq_colors(), limits = c(0, 1),
                                    breaks = c(0, 1), labels = c("none", "high"),
                                    name = "expression\n(per gene)") +
    ggplot2::facet_wrap(~ gene, ncol = ncol %||% min(3, length(genes))) +
    dim_corner_axes(long$dim1, long$dim2, paste0(red, " 1"), paste0(red, " 2")) +
    ggplot2::labs(title = sprintf("Expression on the %s", red),
                  caption = if (length(miss)) paste("Not in the data:", paste(miss, collapse = ", "))) +
    omicone_dim_theme() +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold.italic", hjust = 0.5))
}

#' Expression distributions per group, one panel per gene
#' @param obj Seurat object. @param genes Genes. @param group_by Grouping.
#' @param colors Named colours (default [group_colors()]).
#' @keywords internal
sc_violin_plot <- function(obj, genes, group_by, colors = NULL) {
  ex <- sc_expr_df(obj, genes)
  md <- obj_meta(obj)
  cols <- colors %||% group_colors(obj, group_by)
  grp <- factor(as.character(md[rownames(ex), group_by]), levels = names(cols))
  long <- do.call(rbind, lapply(colnames(ex), function(g) data.frame(group = grp, gene = g, value = ex[[g]])))
  long$gene <- factor(long$gene, levels = colnames(ex))
  set.seed(1)
  pts <- if (nrow(ex) <= 4000) long else long[sample.int(nrow(long), min(nrow(long), 4000 * ncol(ex))), ]
  ggplot2::ggplot(long, ggplot2::aes(.data$group, .data$value, fill = .data$group)) +
    ggplot2::geom_violin(scale = "width", colour = NA, alpha = 0.92, adjust = 1.2, trim = TRUE) +
    ggplot2::geom_jitter(data = pts, width = 0.18, height = 0, size = 0.2, alpha = 0.14,
                         colour = style_tokens()$ink, show.legend = FALSE) +
    ggplot2::stat_summary(fun = stats::median, geom = "point", size = 1.6, colour = "white", show.legend = FALSE) +
    scale_group(cols, "fill", guide = "none") +
    ggplot2::facet_wrap(~ gene, scales = "free_y", ncol = if (length(genes) > 3) 2 else 1) +
    ggplot2::labs(x = NULL, y = "normalised expression", title = sprintf("Expression by %s", group_by)) +
    omicone_theme(grid = "y") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = if (length(cols) > 6) 45 else 0,
                                                       hjust = if (length(cols) > 6) 1 else 0.5),
                   strip.text = ggplot2::element_text(face = "bold.italic"))
}

#' Mean expression and fraction expressing per group (computed, not Seurat's)
#' @param obj Seurat object. @param genes Genes. @param group_by Grouping.
#' @return data.frame(group, gene, pct, mean, scaled): `mean` is the log of the
#'   mean of expm1(normalised) (Seurat's convention), `scaled` its z-score
#'   across groups, clipped to +-2.5.
#' @keywords internal
sc_dot_data <- function(obj, genes, group_by) {
  ex <- sc_expr_df(obj, genes)
  md <- obj_meta(obj)
  lv <- names(group_colors(obj, group_by))
  grp <- factor(as.character(md[rownames(ex), group_by]), levels = lv)
  rows <- list()
  for (g in colnames(ex)) {
    v <- ex[[g]]
    pct <- tapply(v > 0, grp, mean)
    mn <- log1p(tapply(expm1(v), grp, mean))
    s <- stats::sd(mn, na.rm = TRUE)
    z <- if (is.finite(s) && s > 0) (mn - mean(mn, na.rm = TRUE)) / s else mn * 0
    rows[[g]] <- data.frame(group = factor(names(pct), levels = lv), gene = g, pct = 100 * as.numeric(pct),
                            mean = as.numeric(mn), scaled = pmax(-2.5, pmin(2.5, as.numeric(z))))
  }
  out <- do.call(rbind, unname(rows))
  out$gene <- factor(out$gene, levels = unique(colnames(ex)))
  out
}

#' Dot plot: dot size = % of cells expressing, colour = scaled mean expression
#' @param obj Seurat object. @param genes Genes (in display order).
#' @param group_by Grouping. @param gene_groups Optional named vector gene ->
#'   panel label (e.g. the cluster each marker belongs to).
#' @param title Plot title.
#' @keywords internal
sc_dot_plot <- function(obj, genes, group_by, gene_groups = NULL, title = NULL) {
  d <- sc_dot_data(obj, unique(genes), group_by)
  d$group <- factor(d$group, levels = rev(levels(d$group)))
  p <- ggplot2::ggplot(d, ggplot2::aes(.data$gene, .data$group)) +
    ggplot2::geom_point(ggplot2::aes(size = .data$pct, colour = .data$scaled)) +
    ggplot2::scale_size_area(max_size = 6, limits = c(0, 100), breaks = c(25, 50, 75, 100),
                             name = "% expressing") +
    ggplot2::scale_colour_gradientn(colours = style_div_colors(), limits = c(-2.5, 2.5),
                                    name = "scaled\nmean") +
    ggplot2::labs(x = NULL, y = NULL, title = title %||% sprintf("Expression by %s", group_by)) +
    omicone_theme(grid = "xy") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, face = "italic"),
                   axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
  if (!is.null(gene_groups)) {
    d$panel <- factor(gene_groups[as.character(d$gene)], levels = unique(gene_groups[unique(genes)]))
    p$data <- d
    p <- p + ggplot2::facet_grid(~ panel, scales = "free_x", space = "free_x") +
      ggplot2::theme(strip.text = ggplot2::element_text(hjust = 0.5, angle = 0),
                     panel.spacing.x = grid::unit(4, "pt"))
  }
  p
}

#' Heatmap of scaled expression per cell, cells grouped, with a group colour bar
#'
#' At most `max_per_group` cells per group are drawn (a random, seeded
#' sample), so the figure stays readable and fast on large datasets.
#' @param obj Seurat object. @param genes Genes. @param group_by Grouping.
#' @param max_per_group Cells per group drawn.
#' @return A [compose_grid()] figure.
#' @keywords internal
sc_heatmap <- function(obj, genes, group_by, max_per_group = 150) {
  ex <- sc_expr_df(obj, genes)
  md <- obj_meta(obj)
  cols <- group_colors(obj, group_by)
  grp <- factor(as.character(md[rownames(ex), group_by]), levels = names(cols))
  set.seed(1)
  keep <- unlist(lapply(split(seq_along(grp), grp), function(i) if (length(i) > max_per_group) sample(i, max_per_group) else i))
  ex <- ex[keep, , drop = FALSE]
  grp <- grp[keep]
  ord <- order(grp)
  ex <- ex[ord, , drop = FALSE]
  grp <- grp[ord]
  z <- scale(as.matrix(ex))
  z[!is.finite(z)] <- 0
  z[z > 2.5] <- 2.5
  z[z < -2.5] <- -2.5
  cell <- factor(seq_len(nrow(z)))
  long <- data.frame(cell = rep(cell, ncol(z)), gene = factor(rep(colnames(z), each = nrow(z)), levels = rev(colnames(z))),
                     z = as.vector(z))
  bounds <- cumsum(table(grp))
  mids <- bounds - table(grp) / 2
  bar <- data.frame(cell = cell, group = grp)
  p_bar <- ggplot2::ggplot(bar, ggplot2::aes(.data$cell, 1, fill = .data$group)) +
    ggplot2::geom_raster() +
    scale_group(cols, "fill", name = group_by) +
    ggplot2::annotate("text", x = as.numeric(mids), y = 1, label = names(mids), size = 2.9,
                      colour = "white", fontface = "bold") +
    ggplot2::scale_x_discrete(expand = c(0, 0)) + ggplot2::scale_y_continuous(expand = c(0, 0)) +
    omicone_theme(grid = "none") +
    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.title = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(), axis.line = ggplot2::element_blank(),
                   legend.position = "none", plot.margin = ggplot2::margin(0, 12, 2, 10))
  p_hm <- ggplot2::ggplot(long, ggplot2::aes(.data$cell, .data$gene, fill = .data$z)) +
    ggplot2::geom_raster() +
    ggplot2::geom_vline(xintercept = as.numeric(bounds[-length(bounds)]) + 0.5, colour = "white", linewidth = 0.6) +
    ggplot2::scale_fill_gradientn(colours = style_div_colors(), limits = c(-2.5, 2.5), name = "z-score") +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::labs(x = sprintf("cells (up to %d per group)", max_per_group), y = NULL) +
    omicone_theme(grid = "none") +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                   axis.line = ggplot2::element_blank(), axis.text.y = ggplot2::element_text(face = "italic"),
                   legend.position = "none", plot.margin = ggplot2::margin(0, 12, 8, 10))
  compose_grid(list(p_bar, p_hm), nrow = 2, ncol = 1, heights = c(0.6, max(4, length(genes) * 0.5)),
               legends = list(p_hm + ggplot2::theme(legend.position = "bottom")),
               title = sprintf("Scaled expression by %s", group_by))
}
