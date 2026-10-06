#' Figure style: one theme, one colour registry, one panel composer
#'
#' Every figure of the app is drawn with ggplot2 and follows these rules
#' (AGENTS.md section 5.2):
#'
#' * **One theme**, [omicone_theme()]; embeddings use [omicone_dim_theme()]
#'   (no axes, small corner arrows).
#' * **A level keeps its colour everywhere.** Colours of a grouping (clusters,
#'   cell types, samples, conditions, clinical groups) come from
#'   [group_colors()], looked up by level *name*, never by position. The map of
#'   a grouping is stored in the object (`obj@misc$omicone_colors`) by the step
#'   that creates it ([set_group_colors()]), so it survives subsetting and
#'   export; cell types inherit the colour of the cluster they were named from
#'   ([inherit_group_colors()]).
#' * **Fixed meanings keep fixed colours**: kept / removed, up / down, mutant /
#'   wild-type ([style_tokens()]); expression uses one sequential ramp,
#'   z-scores and fold changes one diverging ramp; mutation classes use the
#'   maftools / COSMIC conventions people already read.
#' * **Multi-panel figures** are built with [compose_grid()], which aligns the
#'   panels itself (no patchwork version dependency).
#'
#' @name fct_style
#' @keywords internal
NULL

#' Semantic colours shared by every figure
#' @keywords internal
style_tokens <- function() {
  list(ink = "#1c2530", text = "#3d4a57", muted = "#66727f", faint = "#9aa5b1",
       grid = "#e6eaef", panel = "#f4f6f9", accent = "#2f81c7", accent2 = "#2669a8",
       kept = "#3b6ea5", removed = "#c1476b", up = "#c1476b", down = "#2f81c7",
       ns = "#c3cad3", mutant = "#c1476b", wt = "#8492a0", highlight = "#f4a261",
       ok = "#2f9e6d", warn = "#b7791f", unassigned = "#c3cad3")
}

#' Sequential ramp for expression / density / counts (low = near-white)
#' @keywords internal
style_seq_colors <- function() c("#eef1f5", "#c6dbef", "#6baed6", "#2171b5", "#08306b")

#' Diverging ramp for z-scores and fold changes (blue = low, red = high)
#' @keywords internal
style_div_colors <- function() c("#2166ac", "#67a9cf", "#f5f6f7", "#ef8a62", "#b2182b")

#' Shared ggplot theme
#' @param base_size Base font size (pt).
#' @param grid "y", "x", "xy" or "none": which major grid lines to keep.
#' @keywords internal
omicone_theme <- function(base_size = 12, grid = "xy") {
  tk <- style_tokens()
  gl <- ggplot2::element_line(colour = tk$grid, linewidth = 0.3)
  bl <- ggplot2::element_blank()
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = tk$ink),
      plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.08), hjust = 0,
                                         margin = ggplot2::margin(b = 4)),
      plot.title.position = "plot",
      plot.subtitle = ggplot2::element_text(colour = tk$muted, size = ggplot2::rel(0.86),
                                            margin = ggplot2::margin(b = 8)),
      plot.caption = ggplot2::element_text(colour = tk$muted, size = ggplot2::rel(0.76), hjust = 0),
      plot.caption.position = "plot",
      axis.title = ggplot2::element_text(colour = tk$muted, size = ggplot2::rel(0.88)),
      axis.text = ggplot2::element_text(colour = tk$text, size = ggplot2::rel(0.82)),
      axis.line = ggplot2::element_line(colour = tk$faint, linewidth = 0.35),
      axis.ticks = ggplot2::element_line(colour = tk$faint, linewidth = 0.3),
      axis.ticks.length = grid::unit(2.5, "pt"),
      panel.grid.major.x = if (grid %in% c("x", "xy")) gl else bl,
      panel.grid.major.y = if (grid %in% c("y", "xy")) gl else bl,
      panel.grid.minor = bl,
      strip.text = ggplot2::element_text(face = "bold", hjust = 0, size = ggplot2::rel(0.88),
                                         colour = tk$ink, margin = ggplot2::margin(4, 2, 4, 2)),
      legend.title = ggplot2::element_text(size = ggplot2::rel(0.82), colour = tk$muted),
      legend.text = ggplot2::element_text(size = ggplot2::rel(0.82), colour = tk$text),
      legend.key.size = grid::unit(11, "pt"),
      legend.position = "right",
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(10, 12, 8, 10))
}

#' Theme for embeddings (UMAP / t-SNE / PCA maps): no axes, square panel
#' @keywords internal
omicone_dim_theme <- function(base_size = 12) {
  bl <- ggplot2::element_blank()
  omicone_theme(base_size, grid = "none") +
    ggplot2::theme(axis.line = bl, axis.ticks = bl, axis.text = bl, axis.title = bl,
                   aspect.ratio = 1)
}

#' Small corner arrows naming the two embedding axes (instead of full axes)
#' @param x,y Coordinates drawn. @param xlab,ylab Axis names.
#' @return A list of ggplot layers.
#' @keywords internal
dim_corner_axes <- function(x, y, xlab, ylab) {
  tk <- style_tokens()
  xr <- range(x, na.rm = TRUE)
  yr <- range(y, na.rm = TRUE)
  x0 <- xr[1] - 0.03 * diff(xr)
  y0 <- yr[1] - 0.03 * diff(yr)
  lx <- 0.2 * diff(xr)
  ly <- 0.2 * diff(yr)
  arr <- grid::arrow(length = grid::unit(5, "pt"), type = "closed")
  list(
    ggplot2::annotate("segment", x = x0, y = y0, xend = x0 + lx, yend = y0, colour = tk$muted,
                      linewidth = 0.45, arrow = arr),
    ggplot2::annotate("segment", x = x0, y = y0, xend = x0, yend = y0 + ly, colour = tk$muted,
                      linewidth = 0.45, arrow = arr),
    ggplot2::annotate("text", x = x0, y = y0 - 0.025 * diff(yr), label = xlab, hjust = 0, vjust = 1,
                      size = 3, colour = tk$muted),
    ggplot2::annotate("text", x = x0 - 0.025 * diff(xr), y = y0, label = ylab, hjust = 0, vjust = 0,
                      angle = 90, size = 3, colour = tk$muted))
}

# ---- colour registry -----------------------------------------------------------

#' Levels of a grouping in their natural order
#'
#' Factor levels (present ones) when the column is a factor; otherwise a
#' natural sort ("2" before "10", "C2" before "C10"), with "Unassigned"-type
#' levels last.
#' @param x A vector.
#' @keywords internal
level_order <- function(x) {
  if (is.factor(x)) return(levels(droplevels(x)))
  u <- unique(as.character(x[!is.na(x)]))
  if (!length(u)) return(u)
  key <- vapply(u, function(v) {
    m <- gregexpr("[0-9]+", v)[[1]]
    if (m[1] == -1) return(tolower(v))
    num <- regmatches(v, list(m))[[1]]
    regmatches(v, list(m)) <- list(formatC(as.numeric(num), width = 15, format = "d", flag = "0"))
    tolower(v)
  }, character(1))
  u[order(is_unassigned_level(u), key)]
}

# levels that mean "no label": drawn in neutral grey, never take a palette slot
is_unassigned_level <- function(x) {
  tolower(trimws(x)) %in% c("unassigned", "unknown", "na", "other", "others", "unclassified",
                            "low quality", "doublet", "none", "")
}

#' Complete a stored colour map for a set of levels
#'
#' Levels already in `stored` keep their colour; new levels take palette
#' colours not used yet, in palette order; "Unassigned"-type levels are grey.
#' @param stored Named character vector (level -> colour) or NULL.
#' @param levels Levels to colour.
#' @keywords internal
extend_colors <- function(stored, levels) {
  out <- stats::setNames(rep(NA_character_, length(levels)), levels)
  hit <- intersect(levels, names(stored))
  out[hit] <- stored[hit]
  grey <- is.na(out) & is_unassigned_level(levels)
  out[grey] <- style_tokens()$unassigned
  need <- sum(is.na(out))
  if (need) {
    pool <- sc_palette(max(24, length(levels) + length(stored)))
    pool <- setdiff(pool, out[!is.na(out)])
    pool <- setdiff(pool, unname(stored))
    if (length(pool) < need) pool <- c(pool, grDevices::colorRampPalette(sc_palette(24))(need))
    out[is.na(out)] <- pool[seq_len(need)]
  }
  out
}

#' Fixed colours of groupings whose levels have a meaning
#'
#' Cell-cycle phase, malignancy calls, doublet calls and QC flags look the
#' same in every figure; a stored map still wins over these.
#' @keywords internal
style_known_groups <- function() {
  tk <- style_tokens()
  list(Phase = c(G1 = "#7aa6d6", S = "#f4a261", G2M = "#c1476b"),
       malignant = c(normal = tk$kept, malignant = tk$removed),
       scDblFinder.class = c(singlet = tk$kept, doublet = tk$removed),
       doublet_class = c(singlet = tk$kept, doublet = tk$removed))
}

#' Colours of a grouping, by level name
#'
#' @param x A Seurat object (uses its stored map) or a data.frame / metadata
#'   (deterministic default from the level order).
#' @param col Column name.
#' @param levels Levels to return (default: those present, natural order).
#' @return Named character vector level -> colour.
#' @keywords internal
group_colors <- function(x, col, levels = NULL) {
  md <- if (is.data.frame(x)) x else obj_meta(x)
  if (is.null(levels)) levels <- if (!is.null(md) && col %in% names(md)) level_order(md[[col]]) else character(0)
  stored <- if (!is.data.frame(x)) tryCatch(obj_misc(x, "omicone_colors")[[col]], error = function(e) NULL)
  extend_colors(stored %||% style_known_groups()[[col]], levels)
}

#' Colours of a plain vector of values (no object; deterministic)
#' @param values Values (or a factor). @param levels Optional level order.
#' @keywords internal
value_colors <- function(values, levels = NULL) {
  levels <- levels %||% level_order(values)
  low <- tolower(levels)
  if (length(levels) == 2 && all(sort(low) %in% c("mutant", "wt"))) {
    tk <- style_tokens()
    return(stats::setNames(ifelse(low == "mutant", tk$mutant, tk$wt), levels))
  }
  extend_colors(NULL, levels)
}

#' Store the colour map of a grouping in the object
#' @param obj Seurat object. @param col Column.
#' @param colors Named colours (default: [group_colors()], i.e. keep stored,
#'   extend for new levels).
#' @return The object.
#' @keywords internal
set_group_colors <- function(obj, col, colors = NULL) {
  cols <- colors %||% group_colors(obj, col)
  m <- obj@misc[["omicone_colors"]] %||% list()
  m[[col]] <- cols
  obj@misc[["omicone_colors"]] <- m
  obj
}

#' Colours for a grouping derived from another (cell types from clusters)
#'
#' Each level of `col` takes the colour of the `from` level that contributes
#' most of its cells (a cell type named from cluster 3 is drawn in cluster 3's
#' colour), unless another level already took that colour; the rest get unused
#' palette colours.
#' @param obj Seurat object. @param col New grouping. @param from Source grouping.
#' @return Named colours for `col`.
#' @keywords internal
inherit_group_colors <- function(obj, col, from) {
  md <- obj_meta(obj)
  lv <- level_order(md[[col]])
  src <- group_colors(obj, from)
  tab <- table(factor(as.character(md[[col]]), levels = lv), as.character(md[[from]]))
  out <- stats::setNames(rep(NA_character_, length(lv)), lv)
  used <- character(0)
  sizes <- rowSums(tab)
  for (l in lv[order(-sizes)]) {                       # big populations choose first
    if (is_unassigned_level(l)) next
    cand <- colnames(tab)[order(-tab[l, ])]
    cand <- cand[tab[l, cand] > 0]
    cc <- setdiff(unname(src[cand]), used)
    if (length(cc)) {
      out[l] <- cc[1]
      used <- c(used, cc[1])
    }
  }
  out <- out[!is.na(out)]
  extend_colors(out, lv)
}

#' Give the cell-type column the colours of the clusters it was named from
#' @param obj Seurat object. @param col Cell-type column.
#' @return The object, with the map stored.
#' @keywords internal
color_celltypes <- function(obj, col = "celltype") {
  if (!col %in% obj_meta_cols(obj)) return(obj)
  from <- obj_misc(obj, "omicone_cluster_col") %||% "seurat_clusters"
  cols <- if (from %in% obj_meta_cols(obj)) inherit_group_colors(obj, col, from) else group_colors(obj, col)
  set_group_colors(obj, col, cols)
}

#' ggplot fill / colour scale from a named colour map
#' @param colors Named colours. @param aesthetic "fill" or "colour".
#' @param ... Passed to `scale_*_manual()`.
#' @keywords internal
scale_group <- function(colors, aesthetic = "fill", ...) {
  ggplot2::scale_discrete_manual(aesthetics = aesthetic, values = colors, breaks = names(colors),
                                 limits = names(colors), drop = TRUE,
                                 na.value = style_tokens()$unassigned, ...)
}

# ---- multi-panel composition -----------------------------------------------------

# absolute size (cm) of a unit vector, on a null device if none is open
unit_cm <- function(u, dir = c("width", "height")) {
  dir <- match.arg(dir)
  if (!length(u)) return(0)
  f <- if (dir == "width") grid::convertWidth else grid::convertHeight
  if (grDevices::dev.cur() == 1) return(with_null_device(sum(f(u, "cm", valueOnly = TRUE))))
  sum(f(u, "cm", valueOnly = TRUE))
}

#' The legend of a ggplot as a grob (NULL when it has none)
#' @param p A ggplot. @param position Legend position to draw it at.
#' @keywords internal
plot_legend <- function(p, position = "bottom") {
  if (grDevices::dev.cur() == 1) return(with_null_device(plot_legend(p, position)))
  g <- ggplot2::ggplotGrob(p + ggplot2::theme(legend.position = position))
  idx <- which(grepl("guide-box", g$layout$name))
  for (i in idx) {
    gr <- g$grobs[[i]]
    if (inherits(gr, "gtable") && length(gr$grobs)) return(gr)
  }
  NULL
}

#' Arrange ggplots in a grid with their panels aligned
#'
#' Plots in the same column get the same left / right decoration width, plots
#' in the same row the same top / bottom height, so their panels line up (an
#' oncoplot's bars sit exactly over its tiles). Works with any ggplot2 >= 3.4
#' without patchwork.
#'
#' @param plots List of ggplots, grobs (placed unaligned, e.g. a nested
#'   compose_grid()) or NULL (empty cell), row-major.
#' @param nrow,ncol Grid size.
#' @param widths,heights Relative sizes of the columns / rows.
#' @param legends Optional list of ggplots whose legends are drawn in a row
#'   under the grid (the plots in `plots` should have `legend.position = "none"`).
#' @param title,subtitle Optional figure title and subtitle.
#' @return A gtable of class `omicone_grid` (draw with [draw_plot_object()]).
#' @keywords internal
compose_grid <- function(plots, nrow, ncol, widths = rep(1, ncol), heights = rep(1, nrow),
                         legends = NULL, title = NULL, subtitle = NULL) {
  stopifnot(length(plots) == nrow * ncol)
  # building grobs needs font metrics, i.e. a device; with none open R would
  # open its default one (an Rplots.pdf that is never closed)
  if (grDevices::dev.cur() == 1) {
    return(with_null_device(compose_grid(plots, nrow, ncol, widths, heights, legends, title, subtitle)))
  }
  # ggplots are aligned; a ready-made grob (e.g. a nested compose_grid()) is
  # placed as it is
  fixed <- vapply(plots, function(p) !is.null(p) && !inherits(p, "ggplot"), logical(1))
  grobs <- lapply(plots, function(p) if (is.null(p)) NULL else if (inherits(p, "ggplot")) ggplot2::ggplotGrob(p) else p)
  ext <- lapply(seq_along(grobs), function(k) {
    g <- grobs[[k]]
    if (is.null(g) || fixed[k]) return(NULL)
    lay <- g$layout[grepl("^panel", g$layout$name), , drop = FALSE]
    list(l = min(lay$l), r = max(lay$r), t = min(lay$t), b = max(lay$b))
  })
  idx <- function(i, j) (i - 1) * ncol + j
  side <- function(g, e, which) {
    switch(which,
           left = unit_cm(g$widths[seq_len(e$l - 1)], "width"),
           right = unit_cm(g$widths[setdiff(seq_along(g$widths), seq_len(e$r))], "width"),
           top = unit_cm(g$heights[seq_len(e$t - 1)], "height"),
           bottom = unit_cm(g$heights[setdiff(seq_along(g$heights), seq_len(e$b))], "height"))
  }
  for (j in seq_len(ncol)) {
    cells <- Filter(function(k) !is.null(ext[[k]]), vapply(seq_len(nrow), idx, 0, j = j))
    if (!length(cells)) next
    lt <- max(vapply(cells, function(k) side(grobs[[k]], ext[[k]], "left"), 0))
    rt <- max(vapply(cells, function(k) side(grobs[[k]], ext[[k]], "right"), 0))
    for (k in cells) {
      g <- grobs[[k]]
      nw <- length(g$widths)
      g$widths[1] <- g$widths[1] + grid::unit(lt - side(g, ext[[k]], "left"), "cm")
      g$widths[nw] <- g$widths[nw] + grid::unit(rt - side(g, ext[[k]], "right"), "cm")
      grobs[[k]] <- g
    }
  }
  for (i in seq_len(nrow)) {
    cells <- Filter(function(k) !is.null(ext[[k]]), vapply(seq_len(ncol), function(j) idx(i, j), 0))
    if (!length(cells)) next
    tp <- max(vapply(cells, function(k) side(grobs[[k]], ext[[k]], "top"), 0))
    bt <- max(vapply(cells, function(k) side(grobs[[k]], ext[[k]], "bottom"), 0))
    for (k in cells) {
      g <- grobs[[k]]
      nh <- length(g$heights)
      g$heights[1] <- g$heights[1] + grid::unit(tp - side(g, ext[[k]], "top"), "cm")
      g$heights[nh] <- g$heights[nh] + grid::unit(bt - side(g, ext[[k]], "bottom"), "cm")
      grobs[[k]] <- g
    }
  }
  tab <- gtable::gtable(widths = grid::unit(widths, "null"), heights = grid::unit(heights, "null"))
  for (i in seq_len(nrow)) for (j in seq_len(ncol)) {
    g <- grobs[[idx(i, j)]]
    if (!is.null(g)) tab <- gtable::gtable_add_grob(tab, g, t = i, l = j, clip = "off", name = paste0("cell-", i, "-", j))
  }
  leg <- Filter(Negate(is.null), lapply(legends %||% list(), plot_legend))
  if (length(leg)) {
    # each legend as wide as it is, with a gap between; the row is centred
    lw <- vapply(leg, function(g) unit_cm(g$widths, "width"), 0)
    row <- gtable::gtable(widths = grid::unit(as.vector(rbind(lw, 0.8))[-2 * length(lw)], "cm"),
                          heights = grid::unit(max(vapply(leg, function(g) unit_cm(g$heights, "height"), 0)) + 0.3, "cm"))
    for (k in seq_along(leg)) row <- gtable::gtable_add_grob(row, leg[[k]], t = 1, l = 2 * k - 1, clip = "off")
    tab <- gtable::gtable_add_rows(tab, row$heights, pos = -1)
    tab <- gtable::gtable_add_grob(tab, row, t = nrow(tab), l = 1, r = ncol(tab), clip = "off", name = "legends")
  }
  tk <- style_tokens()
  if (!is.null(subtitle)) {
    tab <- gtable::gtable_add_rows(tab, grid::unit(1.5, "lines"), pos = 0)
    tab <- gtable::gtable_add_grob(tab, grid::textGrob(subtitle, x = grid::unit(10, "pt"), hjust = 0,
                                                       gp = grid::gpar(col = tk$muted, fontsize = 10.5)),
                                   t = 1, l = 1, r = ncol(tab), name = "subtitle")
  }
  if (!is.null(title)) {
    tab <- gtable::gtable_add_rows(tab, grid::unit(2, "lines"), pos = 0)
    tab <- gtable::gtable_add_grob(tab, grid::textGrob(title, x = grid::unit(10, "pt"), hjust = 0,
                                                       gp = grid::gpar(col = tk$ink, fontsize = 13.5, fontface = "bold")),
                                   t = 1, l = 1, r = ncol(tab), name = "title")
  }
  tab <- gtable::gtable_add_padding(tab, grid::unit(c(6, 8, 6, 6), "pt"))
  tab <- gtable::gtable_add_grob(tab, grid::rectGrob(gp = grid::gpar(fill = "white", col = NA)),
                                 t = 1, l = 1, b = nrow(tab), r = ncol(tab), z = -Inf, name = "background")
  class(tab) <- c("omicone_grid", class(tab))
  tab
}
