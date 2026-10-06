#' Plotting layer: renderers, figure download, palette, text scaling
#'
#' Every preview goes through [render_scop_plot()] (ggplot, compose_grid(),
#' ComplexHeatmap objects), which turns errors into a readable message on the
#' canvas. Figures are ggplot2 in the shared style (fct_style.R); base
#' graphics (maftools' own plots) are not used. [render_step_plot()] / [register_figure_download()] make the
#' same drawing downloadable. The scop plotting wrappers live in fct_scop.R.
#'
#' @name fct_plots
#' @keywords internal
NULL

# The shared theme (omicone_theme()) and the colour registry live in fct_style.R.

#' Categorical palette (curated, colour-blind-aware ordering)
#'
#' One fixed palette for every ggplot preview, so a cluster keeps its colour
#' from step to step. (scop's `palette_scp()` was used here until scop 0.9
#' moved its palettes out of the package; the fallback had silently become the
#' only path, so it is now the palette.)
#' @param n Number of colours.
#' @param type Kept for call compatibility; only "discrete" is used.
#' @keywords internal
sc_palette <- function(n = 8, type = "discrete") {
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
#' @param guide Optional bilingual `list(en =, zh =)` describing what the plot
#'   will show. Rendered as a centred placeholder inside the empty preview area
#'   before the first run, and hidden automatically (pure CSS) once the plot
#'   output has content.
#' @param caption Optional bilingual `list(en =, zh =)` figure caption shown
#'   left of the export row (always visible; the export controls keep their
#'   hover fade).
#' @param scene Explainer animation to show in the guide (a scene key from
#'   explain-sc.js / explain-wes.js). Defaults to the step key of `id`.
#' @keywords internal
preview_plot_ui <- function(id, height = "100%", download = FALSE,
                            guide = NULL, caption = NULL, scene = NULL) {
  g <- NULL
  if (!is.null(guide)) {
    # the step's looping explainer animation (explain*.js); the module id
    # ("qc-preview" -> "qc") names the scene unless one is given
    scene <- scene %||% sub("-.*$", "", id)
    g <- shiny::div(
      class = "omicone-preview-guide",
      explain_canvas(scene),
      shiny::div(class = "omicone-preview-guide-icon omicone-explain-fallback", "\U0001F4CA"),
      shiny::div(class = "omicone-preview-guide-text",
                 if (is.list(guide)) i18n(guide$en, guide$zh) else guide),
      shiny::div(class = "omicone-preview-guide-hint",
                 i18n("Set the options on the left, then click the Run button.",
                      "在左侧设置选项，然后点击运行按钮。"))
    )
  }
  # The figure box has a real size before anything is drawn: an empty plot
  # output must never collapse to 0 px (Shiny would open its PNG device at that
  # size), so the animation lies over the box instead of beside it.
  out <- shiny::div(class = "omicone-figbox", g,
                    shiny::plotOutput(id, height = if (identical(height, "100%")) NULL else height))
  g <- NULL
  cap <- NULL
  if (!is.null(caption)) {
    cap <- shiny::div(class = "omicone-figcap",
                      if (is.list(caption)) i18n(caption$en, caption$zh)
                      else caption)
  }
  if (!isTRUE(download)) return(shiny::tagList(g, out, cap))
  dl <- shiny::div(
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
  # caption stays always-visible on the left; export controls keep their
  # hover-fade on the right
  if (!is.null(cap)) {
    return(shiny::tagList(g, out,
                          shiny::div(class = "omicone-figfoot", cap, dl)))
  }
  shiny::tagList(g, out, dl)
}

#' The canvas a step's explainer animation is drawn on
#'
#' Drawn client-side by explain.js (no server work, nothing downloaded). An
#' unknown scene hides itself and the static icon next to it shows instead.
#' @param scene Scene key, normally the step key ("qc", "wes_tmb").
#' @keywords internal
explain_canvas <- function(scene) {
  shiny::tags$canvas(class = "omicone-explain", `data-scene` = scene,
                     `aria-hidden` = "true")
}

#' A step's explainer animation with a short bilingual text, for steps whose
#' main output is not a plot (report, export, tables)
#' @param scene Scene key. @param en,zh What the step produces.
#' @keywords internal
explain_scene <- function(scene, en, zh) {
  shiny::div(
    class = "omicone-explain-wrap",
    explain_canvas(scene),
    shiny::div(class = "omicone-preview-guide-text", i18n(en, zh))
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
      before <- grDevices::dev.list()
      if (identical(f, "pdf")) {
        grDevices::pdf(file, width = w, height = h, useDingbats = FALSE)
      } else if (identical(f, "jpg")) {
        grDevices::jpeg(file, width = w, height = h, units = "in",
                        res = d, quality = 95)
      } else {
        grDevices::png(file, width = w, height = h, units = "in", res = d)
      }
      dev <- grDevices::dev.cur()
      # close exactly this device (and anything the drawing left open), even
      # if the drawing switched the current device
      on.exit(close_devices_since(before), add = TRUE)
      msg <- tryCatch({ draw_fn(); NULL },
                      shiny.silent.error = function(e)
                        "Nothing to export yet — run this step first.",
                      error = function(e) conditionMessage(e))
      if (dev %in% grDevices::dev.list()) grDevices::dev.set(dev)
      if (!is.null(msg)) draw_plot_message(paste0("Plot error:\n", msg))
    }
  )
}

#' Render a scop/ggplot/ComplexHeatmap object to a Shiny plot output
#'
#' Accepts whatever a scop plotting function returns: a ggplot/patchwork object
#' (printed) or a ComplexHeatmap (drawn). `plot_expr` is a function returning the
#' plot object; it runs inside tryCatch so failures show as a message.
#' @keywords internal
render_scop_plot <- function(plot_expr, id = NULL) {
  # Never render below 50 px: a collapsed or hidden plot box would otherwise
  # get a PNG device a few pixels tall ("figure margins too large", or an
  # invalid device). Sizes are read the way renderPlot() reads them.
  sess <- shiny::getDefaultReactiveDomain()
  dim_fn <- function(which, default) {
    if (is.null(id) || is.null(sess)) return("auto")
    key <- paste0("output_", sess$ns(id), "_", which)
    function() max(50, sess$clientData[[key]] %||% default)
  }
  shiny::renderPlot({
    # Shiny's device is the current one here. Whatever the plotting code does,
    # every device opened from now on is closed when the render ends (also
    # when req() stops it), and the drawing goes to Shiny's device.
    dev <- grDevices::dev.cur()
    before <- grDevices::dev.list()
    on.exit(close_devices_since(before, keep = dev), add = TRUE)
    # Build the plot object. req()/validate() (no data yet) stay silent.
    p <- tryCatch(
      plot_expr(),
      shiny.silent.error = function(e) NULL,
      error = function(e) structure(list(msg = conditionMessage(e)), class = "omicone_plot_error"))
    shiny::req(!is.null(p))
    if (dev %in% grDevices::dev.list()) grDevices::dev.set(dev)
    if (any(grDevices::dev.size("px") < 40)) return(invisible())

    # Draw it. Any drawing error is turned into a readable message ON the canvas
    # (and a toast) so the user never sees an opaque "[object Object]".
    show_err <- function(msg) {
      shiny::showNotification(paste("Plot error:", msg), type = "error", duration = 12)
      if (dev %in% grDevices::dev.list()) grDevices::dev.set(dev)
      draw_plot_message(paste0("Plot error:\n", msg))
    }
    if (inherits(p, "omicone_plot_error")) { show_err(p$msg); return(invisible()) }
    tryCatch(draw_plot_object(p), error = function(e) show_err(conditionMessage(e)))
  }, width = dim_fn("width", 600), height = dim_fn("height", 400))
}

#' A message drawn in place of a figure (grid, so it works at any size)
#' @param msg Text. @keywords internal
draw_plot_message <- function(msg) {
  grid::grid.newpage()
  grid::grid.text(msg, gp = grid::gpar(col = style_tokens()$removed, fontsize = 12))
  invisible(NULL)
}

#' Close the graphics devices opened since `before`
#'
#' The one way the app closes devices. A bare `dev.off()` closes whichever
#' device happens to be current and then makes the *next* one current, so a
#' single misplaced call can close Shiny's device and leave ours open; leaked
#' devices pile up until R refuses to open more ("too many open devices").
#' @param before `dev.list()` taken before opening. @param keep Device to keep
#'   open and current (e.g. Shiny's plot device), or NULL.
#' @keywords internal
close_devices_since <- function(before, keep = NULL) {
  now <- grDevices::dev.list()
  extra <- setdiff(now, c(before, keep))
  for (d in rev(extra)) if (d %in% grDevices::dev.list()) grDevices::dev.off(d)
  if (!is.null(keep) && keep %in% grDevices::dev.list()) grDevices::dev.set(keep)
  invisible(NULL)
}

#' Evaluate an expression on a throw-away null device
#'
#' For code that needs a device without drawing anything visible (font metrics,
#' maftools functions that plot as a side effect, a dry run): opens `pdf(NULL)`,
#' evaluates, then closes that device and anything opened meanwhile, and makes
#' the previously current device current again.
#' @param expr Expression.
#' @keywords internal
with_null_device <- function(expr) {
  prev <- grDevices::dev.cur()
  before <- grDevices::dev.list()
  grDevices::pdf(NULL)
  on.exit({
    close_devices_since(before)
    if (prev > 1 && prev %in% grDevices::dev.list()) grDevices::dev.set(prev)
  }, add = TRUE)
  force(expr)
}

# Backward-compatible alias used by older modules (renders a ggplot expr).
#' @keywords internal
render_preview_plot <- function(gg_expr, tooltip = "text") {
  render_scop_plot(gg_expr)
}

#' Draw a plot object on the current device (ggplot / patchwork / compose_grid / Heatmap)
#' @param p A plot object.
#' @keywords internal
draw_plot_object <- function(p) {
  if (inherits(p, "gtable")) {                      # compose_grid() figures
    grid::grid.newpage()
    grid::grid.draw(p)
  } else if (methods::is(p, "Heatmap") || methods::is(p, "HeatmapList")) {
    if (has_pkg("ComplexHeatmap")) ComplexHeatmap::draw(p) else print(p)
  } else {
    print(p)
  }
  invisible(NULL)
}

#' A step's preview plot, on screen and as a download, from one closure
#'
#' The standard way a module renders a ggplot / ComplexHeatmap figure: assigns
#' `output[[id]]` with [render_scop_plot()] and wires the export row of
#' `preview_plot_ui(id, download = TRUE)` to the same `plot_fn`, so the
#' downloaded file is the figure on screen.
#'
#' @param output,input Module server `output` / `input`.
#' @param id Local output id (as given to `preview_plot_ui()`).
#' @param plot_fn Zero-argument function returning a plot object; may call
#'   `req()` while there is nothing to draw.
#' @param name File stem for downloads. @param width,height Inches (or
#'   zero-argument functions returning inches).
#' @keywords internal
render_step_plot <- function(output, input, id, plot_fn, name = id,
                             width = 10, height = 7) {
  output[[id]] <- render_scop_plot(plot_fn, id = id)
  register_figure_download(output, input, id, function() {
    dev <- grDevices::dev.cur()
    p <- plot_fn()
    shiny::req(!is.null(p))
    if (dev %in% grDevices::dev.list()) grDevices::dev.set(dev)
    draw_plot_object(p)
  }, name = name, width = width, height = height)
  invisible(NULL)
}

# The scop plotting wrappers (CellDimPlot, FeatureDimPlot, ...) and the
# mascarade outline overlay live in fct_scop.R.
