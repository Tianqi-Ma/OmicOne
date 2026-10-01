#' Plotting layer — wraps scop's plotting functions for a unified look
#'
#' All previews use scop's plotting functions (CellDimPlot, FeatureDimPlot,
#' GroupHeatmap, DynamicHeatmap, VolcanoPlot, EnrichmentPlot, ...) with
#' `scop::palette_scp()` colours, so every figure matches the scop/SCP aesthetic.
#' Each wrapper is gated by `require_pkgs("scop")` and wrapped in tryCatch so a
#' missing package or a signature mismatch surfaces as a friendly message rather
#' than crashing the app.
#'
#' NOTE: scop signatures are verified at runtime on the user's machine; a few
#' argument names may need adjustment against the installed scop version.
#'
#' @name fct_plots
#' @keywords internal
NULL

#' Shared minimal ggplot theme (fallback when not using a scop plot)
#' @keywords internal
omicone_theme <- function() {
  ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.25, colour = "#8b98a533"),
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )
}

#' Categorical palette — scop's palette_scp when available, else a curated set
#' @param n Number of colours.
#' @param type "discrete", "continuous", or "diverging".
#' @keywords internal
sc_palette <- function(n = 8, type = "discrete") {
  if (has_pkg("scop")) {
    out <- tryCatch(scop::palette_scp(seq_len(n), n = n, type = type),
                    error = function(e) NULL)
    if (!is.null(out)) return(unname(out))
  }
  base <- c("#2f81c7", "#e4572e", "#3fb37f", "#b5179e", "#f4a261", "#4361ee",
            "#e63946", "#2a9d8f", "#9c6ade", "#ffca3a", "#577590", "#d68fb0",
            "#43aa8b", "#f9844a", "#277da1", "#f94144", "#90be6d", "#845ec2",
            "#ff9f1c", "#4d908e", "#c9184a", "#00b4d8", "#bc6c25", "#606c38")
  if (n <= length(base)) return(base[seq_len(n)])
  grDevices::colorRampPalette(base)(n)
}

#' UI slot for a (large) preview plot. scop plots are static high-res images.
#'
#' With `download = TRUE` a compact export row is added under the plot: a
#' download button, a format selector (PNG / JPEG / PDF) and — for raster
#' formats — a DPI field. The companion server side is
#' [register_figure_download()], wired with the same local `id`.
#' @param id Namespaced output id. @param height CSS height.
#' @param download Add the figure-export row under the plot.
#' @keywords internal
preview_plot_ui <- function(id, height = "100%", download = FALSE) {
  out <- shiny::plotOutput(id, height = height)
  if (!isTRUE(download)) return(out)
  shiny::tagList(
    out,
    shiny::div(
      class = "omicone-fig-dl",
      shiny::downloadButton(paste0(id, "_dl"), i18n("Download figure", "下载图片"),
                            class = "btn-sm"),
      shiny::selectInput(paste0(id, "_dlfmt"), NULL,
                         c("PNG" = "png", "JPEG" = "jpg", "PDF (vector)" = "pdf"),
                         selectize = FALSE, width = "118px"),
      shiny::conditionalPanel(
        sprintf("input['%s'] != 'pdf'", paste0(id, "_dlfmt")),
        shiny::div(class = "omicone-fig-dl-dpi",
                   shiny::tags$span("DPI", class = "omicone-fig-dl-unit"),
                   shiny::numericInput(paste0(id, "_dpi"), NULL, value = 300,
                                       min = 72, max = 1200, step = 50,
                                       width = "88px")))
    )
  )
}

#' Register a figure download handler for a preview plot
#'
#' The server half of `preview_plot_ui(..., download = TRUE)`. Replays
#' `draw_fn()` — the same closure the on-screen renderer uses — into a
#' png/jpeg/pdf device, so the exported file matches what is on screen.
#'
#' @param output,input Module server `output` / `input`.
#' @param id Local (un-namespaced) output id of the plot.
#' @param draw_fn Zero-argument function that draws the figure as a side
#'   effect (for ggplot outputs: `function() print(gg)`).
#' @param name File stem for the download.
#' @param width,height Device size in inches; may be a zero-arg function for
#'   plots whose ideal size depends on the data (e.g. an oncoplot).
#' @keywords internal
register_figure_download <- function(output, input, id, draw_fn, name,
                                     width = 10, height = 7) {
  dl  <- paste0(id, "_dl")
  fmt <- paste0(id, "_dlfmt")
  dpi <- paste0(id, "_dpi")
  output[[dl]] <- shiny::downloadHandler(
    filename = function() {
      f <- input[[fmt]] %||% "png"
      sprintf("%s_%s.%s", name, format(Sys.time(), "%Y%m%d_%H%M%S"), f)
    },
    content = function(file) {
      f <- input[[fmt]] %||% "png"
      d <- suppressWarnings(as.numeric(input[[dpi]] %||% 300))
      if (length(d) != 1 || !is.finite(d)) d <- 300
      d <- max(36, min(2400, d))
      w <- if (is.function(width)) width() else width
      h <- if (is.function(height)) height() else height
      if (identical(f, "pdf")) {
        grDevices::pdf(file, width = w, height = h, useDingbats = FALSE)
      } else if (identical(f, "jpg")) {
        grDevices::jpeg(file, width = w, height = h, units = "in",
                        res = d, quality = 95)
      } else {
        grDevices::png(file, width = w, height = h, units = "in", res = d)
      }
      on.exit(grDevices::dev.off(), add = TRUE)
      msg <- tryCatch({ draw_fn(); NULL },
                      shiny.silent.error = function(e)
                        "Nothing to export yet — run this step first.",
                      error = function(e) conditionMessage(e))
      if (!is.null(msg)) {
        p <- graphics::par(mar = c(0, 0, 0, 0)); on.exit(graphics::par(p), add = TRUE)
        graphics::plot.new()
        graphics::text(0.5, 0.5, paste0("Plot error:\n", msg), col = "#c1476b")
      }
    }
  )
}

#' Gentle text scaling: shrink as the number of labels grows
#'
#' `sqrt` falloff clamped to `[lo, hi]`: a handful of labels gets slightly
#' larger text, a crowded plot slightly smaller, without extremes. Used to
#' size maftools fonts from the number of genes / samples / labels shown.
#' @param n Number of labels on the plot.
#' @param base Size returned when `n == n_ref`. @param n_ref Reference count.
#' @param lo,hi Clamp bounds.
#' @keywords internal
adaptive_cex <- function(n, base = 1, n_ref = 20, lo = 0.55, hi = 1.25) {
  n <- suppressWarnings(as.numeric(n %||% NA))
  if (length(n) != 1 || !is.finite(n) || n <= 0) return(base)
  max(lo, min(hi, base * sqrt(n_ref / n)))
}

#' Boost base-graphics text inside a draw closure
#'
#' maftools' absolute cex defaults assume a small device; on a full-width
#' browser panel they render tiny. This wraps a draw function with a temporary
#' `par()` bump of axis / label / title text sizes. Functions that set their
#' own cex values internally keep them (explicit values win over `par`).
#' @param draw_fn Zero-argument draw function.
#' @param cex Target axis text size (labels and titles slightly larger).
#' @keywords internal
with_text_boost <- function(draw_fn, cex = 1.15) {
  force(draw_fn); force(cex)
  function() {
    op <- graphics::par(cex.axis = cex, cex.lab = cex * 1.05,
                        cex.main = cex * 1.1, cex.sub = cex)
    on.exit(graphics::par(op), add = TRUE)
    draw_fn()
  }
}

#' Render a scop/ggplot/ComplexHeatmap object to a Shiny plot output
#'
#' Accepts whatever a scop plotting function returns: a ggplot/patchwork object
#' (printed) or a ComplexHeatmap (drawn). `plot_expr` is a function returning the
#' plot object; it runs inside tryCatch so failures show as a message.
#' @keywords internal
render_scop_plot <- function(plot_expr) {
  shiny::renderPlot({
    # Build the plot object. req()/validate() (no data yet) stay silent.
    p <- tryCatch(
      plot_expr(),
      shiny.silent.error = function(e) NULL,
      error = function(e) structure(list(msg = conditionMessage(e)), class = "omicone_plot_error"))
    shiny::req(!is.null(p))

    # Draw it. Any drawing error is turned into a readable message ON the canvas
    # (and a toast) so the user never sees an opaque "[object Object]".
    show_err <- function(msg) {
      shiny::showNotification(paste("Plot error:", msg), type = "error", duration = 12)
      op <- graphics::par(mar = c(0, 0, 0, 0)); on.exit(graphics::par(op), add = TRUE)
      graphics::plot.new()
      graphics::text(0.5, 0.5, paste0("Plot error:\n", msg), col = "#c1476b", cex = 1.1)
    }
    if (inherits(p, "omicone_plot_error")) { show_err(p$msg); return(invisible()) }
    tryCatch({
      if (methods::is(p, "Heatmap") || methods::is(p, "HeatmapList")) {
        if (has_pkg("ComplexHeatmap")) ComplexHeatmap::draw(p) else print(p)
      } else {
        print(p)
      }
    }, error = function(e) show_err(conditionMessage(e)))
  })
}

# Backward-compatible alias used by older modules (renders a ggplot expr).
#' @keywords internal
render_preview_plot <- function(gg_expr, tooltip = "text") {
  render_scop_plot(gg_expr)
}

#' Render a base-graphics plot (maftools) to a Shiny plot output
#'
#' maftools draws with base graphics and returns nothing useful, so its calls
#' cannot go through [render_scop_plot()], which prints an object. `plot_expr` is
#' a function that draws as a side effect; failures become a readable message on
#' the canvas instead of a blank panel.
#'
#' @param plot_expr Function that draws a plot as a side effect.
#' @param bg Panel background, matched to the app's dark theme.
#' @keywords internal
render_base_plot <- function(plot_expr, bg = "white") {
  shiny::renderPlot({
    op <- graphics::par(bg = bg)
    on.exit(graphics::par(op), add = TRUE)
    ok <- tryCatch({ plot_expr(); TRUE },
                   shiny.silent.error = function(e) NA,
                   error = function(e) conditionMessage(e))
    if (isTRUE(ok)) return(invisible())
    if (is.na(ok)) { shiny::req(FALSE) }        # nothing to draw yet: stay blank
    shiny::showNotification(paste("Plot error:", ok), type = "error", duration = 12)
    p2 <- graphics::par(mar = c(0, 0, 0, 0)); on.exit(graphics::par(p2), add = TRUE)
    graphics::plot.new()
    graphics::text(0.5, 0.5, paste0("Plot error:\n", ok), col = "#c1476b", cex = 1.1)
  })
}

# ---- scop plotting wrappers -------------------------------------------------

#' Dimensional-reduction scatter (clusters / metadata), scop::CellDimPlot,
#' with optional mascarade cell-type outlines.
#' @keywords internal
sc_dimplot <- function(srt, group_by, reduction = NULL, mask = FALSE,
                       palette = "Paired", label = TRUE, ...) {
  if (!require_pkgs("scop", "Dimension plot")) return(NULL)
  p <- scop::CellDimPlot(srt, group.by = group_by, reduction = reduction,
                         palette = palette, label = label, ...)
  if (isTRUE(mask) && has_pkg("mascarade")) {
    p <- tryCatch(add_mascarade(p, srt, group_by, reduction),
                  error = function(e) p)
  }
  p
}

#' Feature (gene / score) on a reduction, scop::FeatureDimPlot
#' @keywords internal
sc_featureplot <- function(srt, features, reduction = NULL, ...) {
  if (!require_pkgs("scop", "Feature plot")) return(NULL)
  scop::FeatureDimPlot(srt, features = features, reduction = reduction, ...)
}

#' Grouped mean-expression heatmap, scop::GroupHeatmap (signature figure)
#' @keywords internal
sc_groupheatmap <- function(srt, features, group_by, ...) {
  if (!require_pkgs("scop", "GroupHeatmap")) return(NULL)
  scop::GroupHeatmap(srt, features = features, group.by = group_by, ...)
}

#' Dynamic (pseudotime) heatmap, scop::DynamicHeatmap
#' @keywords internal
sc_dynamicheatmap <- function(srt, lineages, ...) {
  if (!require_pkgs("scop", "DynamicHeatmap")) return(NULL)
  scop::DynamicHeatmap(srt, lineages = lineages, ...)
}

#' Composition / statistics plot, scop::CellStatPlot
#' @keywords internal
sc_cellstat <- function(srt, stat_by, group_by = NULL, plot_type = "bar", ...) {
  if (!require_pkgs("scop", "Cell statistics")) return(NULL)
  scop::CellStatPlot(srt, stat.by = stat_by, group.by = group_by,
                     plot_type = plot_type, ...)
}

#' Per-group feature distribution, scop::FeatureStatPlot
#' @keywords internal
sc_featurestat <- function(srt, stat_by, group_by, plot_type = "violin", ...) {
  if (!require_pkgs("scop", "Feature statistics")) return(NULL)
  scop::FeatureStatPlot(srt, stat.by = stat_by, group.by = group_by,
                        plot_type = plot_type, ...)
}

#' Volcano plot of DE results, scop::VolcanoPlot
#' @keywords internal
sc_volcano <- function(srt, group_by, ...) {
  if (!require_pkgs("scop", "Volcano plot")) return(NULL)
  scop::VolcanoPlot(srt, group_by = group_by, ...)
}

#' Enrichment plot (GO/KEGG/...), scop::EnrichmentPlot
#' @keywords internal
sc_enrichplot <- function(srt, group_by, plot_type = "bar", ...) {
  if (!require_pkgs("scop", "Enrichment plot")) return(NULL)
  scop::EnrichmentPlot(srt, group_by = group_by, plot_type = plot_type, ...)
}

#' GSEA running-score plot, scop::GSEAPlot
#' @keywords internal
sc_gseaplot <- function(srt, ...) {
  if (!require_pkgs("scop", "GSEA plot")) return(NULL)
  scop::GSEAPlot(srt, ...)
}

#' RNA-velocity stream/grid, scop::VelocityPlot
#' @keywords internal
sc_velocityplot <- function(srt, reduction = NULL, ...) {
  if (!require_pkgs("scop", "Velocity plot")) return(NULL)
  scop::VelocityPlot(srt, reduction = reduction, ...)
}

#' PAGA graph on an embedding, scop::PAGAPlot
#' @keywords internal
sc_pagaplot <- function(srt, ...) {
  if (!require_pkgs("scop", "PAGA plot")) return(NULL)
  scop::PAGAPlot(srt, ...)
}

#' Mascarade cell-type outlines overlaid on a scop dim plot
#'
#' Uses mascarade::generateMask() to compute polygon outlines around each group
#' on the 2D embedding and overlays them on an existing ggplot dim plot.
#' @keywords internal
add_mascarade <- function(p, srt, group_by, reduction = NULL) {
  if (!require_pkgs("mascarade", "Cell-type outlines")) return(p)
  emb <- SeuratObject::Embeddings(srt, reduction = reduction %||% SeuratObject::DefaultDimReduc(srt))
  labels <- obj_meta(srt)[[group_by]]
  mask <- mascarade::generateMask(dims = emb[, 1:2], cluster = labels)
  p + ggplot2::geom_path(
    data = mask,
    ggplot2::aes(x = .data$x, y = .data$y, group = .data$group),
    colour = "grey20", linewidth = 0.4, inherit.aes = FALSE
  )
}
