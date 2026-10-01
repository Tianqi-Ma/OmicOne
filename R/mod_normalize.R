#' Module 4: Normalization
#'
#' Put cells on a comparable scale so that differences reflect biology rather than
#' sequencing depth. Supports classic LogNormalize (default) or SCTransform.
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_normalize
NULL

#' @rdname mod_normalize
#' @keywords internal
mod_normalize_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = list(en = "Normalization", zh = "归一化"),
    what = list(
      en = "Adjust each cell's counts so that cells sequenced to different depths
            become comparable.",
      zh = "调整每个细胞的计数，使测序深度不同的细胞变得可比较。"),
    why  = list(
      en = "Raw counts depend on how deeply each cell was sequenced. Without
            normalization, deeper cells look artificially 'more expressing'; the
            comparisons downstream would reflect depth, not biology.",
      zh = "原始计数取决于每个细胞的测序深度。若不做归一化，测序更深的细胞会显得
            '表达量更高'，下游比较反映的将是深度而非生物学差异。"),
    how  = list(
      en = "<b>LogNormalize</b> scales each cell to a common total, then log
            transforms (robust default). <b>SCT</b> models counts with a
            regularized negative binomial and often needs no extra scaling.",
      zh = "<b>LogNormalize</b> 将每个细胞缩放到统一的总量后再取对数（稳健的默认方法）。
            <b>SCT</b> 用正则化负二项模型对计数建模，通常无需额外缩放。"),
    read = list(
      en = "The before/after view should show the depth effect disappearing:
            after normalization, total counts per cell no longer drive the
            spread.",
      zh = "前后对比应显示深度效应消失：归一化后，每细胞总计数不再主导数据分布。"),
    example = list(
      en = "A cell with 20,000 UMIs and one with 5,000 UMIs are put on the same
               scale so a shared marker reads similarly in both.",
      zh = "一个有 20,000 个 UMI 的细胞与一个有 5,000 个 UMI 的细胞被放到同一尺度，
               使共有的标记基因在两者中读数相近。")
  )
  controls <- shiny::tagList(
    label_with_help("Method",
                    "LogNormalize = classic log-scaled counts (recommended default). SCT = variance-stabilizing transform (SCTransform).",
                    "方法",
                    "LogNormalize = 经典的对数缩放计数（推荐默认）。SCT = 方差稳定变换（SCTransform）。"),
    shiny::selectInput(ns("method"), NULL,
                       choices = c("LogNormalize" = "LogNormalize",
                                   "SCT" = "SCT"),
                       selected = "LogNormalize"),
    shiny::conditionalPanel(
      sprintf("input['%s'] == 'LogNormalize'", ns("method")),
      label_with_help("Scale factor",
                      "Common total each cell is scaled to before log. 10,000 is the standard default.",
                      "缩放因子",
                      "取对数前每个细胞缩放到的统一总量。10,000 是标准默认值。"),
      shiny::numericInput(ns("scale_factor"), NULL, value = 1e4, min = 1, step = 1e3)
    ),
    run_button(ns("run"), "Normalize", "归一化")
  )
  step_container(title = list(en = "Normalization", zh = "归一化"),
                 subtitle = list(en = "Make cells sequenced to different depths comparable.",
                                 zh = "使测序深度不同的细胞可相互比较。"),
                 explainer = explainer, controls = controls,
                 summary = shiny::uiOutput(ns("summary")),
                 preview = preview_plot_ui(ns("preview"),
                   guide = list(en = "Before/after normalization diagnostics will be drawn here.",
                                zh = "运行后，这里将绘制归一化前后的诊断图。"),
                   caption = list(en = "Normalization diagnostics: the depth effect before vs after.",
                                  zh = "归一化诊断：深度效应的前后对比。")))
}

#' @rdname mod_normalize
#' @keywords internal
mod_normalize_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    res <- shiny::reactiveValues(done = FALSE, method = NULL, scale_factor = NULL,
                                 md = NULL)

    shiny::observeEvent(input$run, {
      shiny::req(rv$obj)
      if (!require_pkgs("Seurat", "Normalization")) return(NULL)
      method <- input$method
      sf <- input$scale_factor
      obj <- with_progress_notify({
        o <- normalize_obj(rv$obj, method = method, scale_factor = sf)
        res$md <- obj_meta(o)
        o
      }, message = "Normalizing counts...")
      if (is.null(obj)) return(NULL)
      rv$obj <- obj
      res$done <- TRUE
      res$method <- method
      res$scale_factor <- sf
      mark_done(rv, "normalize")
      log_step(log_rv, "Normalization",
               params = list(method = method,
                             scale.factor = if (method == "LogNormalize") sf else NA),
               code = if (method == "LogNormalize") {
                 sprintf("obj <- Seurat::NormalizeData(obj, normalization.method = 'LogNormalize', scale.factor = %g)", sf)
               } else {
                 "obj <- Seurat::SCTransform(obj)"
               })
      shiny::showNotification(sprintf("Normalization done (%s).", method),
                              type = "message")
    })

    output$summary <- shiny::renderUI({
      if (!isTRUE(res$done)) return(shiny::div(class = "omicone-placeholder",
                                               i18n("Pick a method and click Normalize.",
                                                    "选择一种方法并点击归一化。")))
      bslib::layout_columns(
        col_widths = c(6, 6),
        stat_tile(i18n("Method", "方法"), res$method),
        stat_tile(i18n("Scale factor", "缩放因子"),
                  if (identical(res$method, "LogNormalize"))
                    format(res$scale_factor, big.mark = ",") else "n/a")
      )
    })

    output$preview <- render_preview_plot(function() {
      md <- res$md
      shiny::req(md)
      shiny::req(!is.null(md$nCount_RNA))
      # Robust view of sequencing depth: distribution of per-cell library sizes
      df <- data.frame(cell = rownames(md),
                       nCount = md$nCount_RNA,
                       stringsAsFactors = FALSE)
      df$log_count <- log10(df$nCount + 1)
      df$text <- sprintf("cell=%s\nUMIs=%s", df$cell,
                         format(df$nCount, big.mark = ","))
      ggplot2::ggplot(df, ggplot2::aes(x = log_count)) +
        ggplot2::geom_histogram(bins = 50, fill = sc_palette(1), alpha = 0.85) +
        ggplot2::labs(x = "Library size (log10 UMIs)", y = "Cells",
                      title = "Per-cell sequencing depth / 每个细胞的测序深度") +
        omicone_theme()
    })
  })
}
