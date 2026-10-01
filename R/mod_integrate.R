#' Module: Batch integration (optional)
#'
#' Align cells from different samples/batches so that shared cell types overlap
#' instead of forming separate, technically-driven clumps. Harmony is the
#' one-click default; "none" leaves the data uncorrected.
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_integrate
NULL

#' @rdname mod_integrate
#' @keywords internal
mod_integrate_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = list(en = "Batch integration (optional)", zh = "批次整合（可选）"),
    what = list(
      en = "Correct for technical differences between samples/batches so that the
            same cell type from different samples lines up.",
      zh = "校正样本/批次之间的技术差异，使来自不同样本的相同细胞类型对齐。"),
    why  = list(
      en = "Without correction, cells often cluster by their sample of origin
            rather than by biology, which confounds downstream clustering.",
      zh = "若不校正，细胞常常按来源样本而非生物学聚在一起，从而干扰下游聚类。"),
    how  = list(
      en = "Pick the metadata column that identifies your batch/sample, then a
            method. <b>Harmony</b> is fast and a good default. Choose <b>none</b>
            if you have a single sample or want to inspect the raw structure.",
      zh = "选择标识批次/样本的元数据列，再选择方法。<b>Harmony</b> 快速且是不错的默认。
            若只有单个样本或想查看未校正的原始结构，请选择 <b>none</b>。"),
    read = list(
      en = "Compare before with after: cells of the same type should mix across
            batches while distinct types stay apart. Types collapsing into each
            other means the integration was too aggressive.",
      zh = "对比整合前后：相同类型的细胞应跨批次混合，不同类型应保持分离。类型互相塌缩说明整合强度过大。"),
    example = list(
      en = "Cells from two samples forming two separate clumps that should
               overlap — integration aligns them.",
      zh = "来自两个样本、本应重叠却形成两个独立团块的细胞——整合会将它们对齐。")
  )
  controls <- shiny::tagList(
    label_with_help("Batch / sample column",
                    "The metadata column that identifies each sample or batch.",
                    "批次/样本列",
                    "标识每个样本或批次的元数据列。"),
    shiny::uiOutput(ns("batch_ui")),
    label_with_help("Method",
                    "Harmony is fast and robust. CCA/RPCA are Seurat anchor-based. none = no correction.",
                    "方法",
                    "Harmony 快速且稳健。CCA/RPCA 是 Seurat 基于锚点的方法。none = 不做校正。"),
    shiny::selectInput(ns("method"), NULL,
                       c("none (no correction)" = "none", "Harmony" = "harmony",
                         "CCA" = "CCA", "RPCA" = "RPCA")),
    run_button(ns("run"), "Run integration", "运行整合")
  )
  step_container(title = list(en = "Batch integration", zh = "批次整合"),
                 subtitle = list(en = "Align batches so the same cell type sits together.",
                                 zh = "对齐批次，使相同细胞类型聚在一起。"),
                 explainer = explainer, controls = controls,
                 summary = shiny::uiOutput(ns("summary")),
                 preview = preview_plot_ui(ns("preview"),
                   guide = list(en = "Before/after batch-mixing will be drawn here.",
                                zh = "运行后，这里将绘制整合前后的批次混合对比。"),
                   caption = list(en = "Batch mixing before vs after integration.",
                                  zh = "整合前后的批次混合对比。")))
}

#' @rdname mod_integrate
#' @keywords internal
mod_integrate_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    res <- shiny::reactiveValues(done = FALSE, method = NULL, batch = NULL)

    # Populate the batch column selector dynamically from the object metadata.
    output$batch_ui <- shiny::renderUI({
      cols <- obj_meta_cols(rv$obj)
      if (length(cols) == 0) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("Load a dataset to choose a batch column.",
                               "加载数据集以选择批次列。")))
      }
      default <- guess_batch_col(obj_meta(rv$obj))
      if (is.null(default)) default <- cols[1]
      shiny::selectInput(session$ns("batch"), NULL, choices = cols, selected = default)
    })

    shiny::observeEvent(input$run, {
      shiny::req(rv$obj)
      shiny::req(input$batch)
      pkgs <- if (input$method == "harmony") c("Seurat", "harmony") else "Seurat"
      if (!require_pkgs(pkgs, "Integration")) return(NULL)
      method <- input$method
      batch  <- input$batch
      obj <- with_progress_notify({
        integrate_obj(rv$obj, batch = batch, method = method)
      }, message = "Integrating batches...")
      if (is.null(obj)) return(NULL)
      rv$obj <- obj
      res$done   <- TRUE
      res$method <- method
      res$batch  <- batch
      mark_done(rv, "integrate")
      log_step(log_rv, "Integration",
               params = list(method = method, batch = batch),
               code = sprintf("obj <- integrate_obj(obj, batch = '%s', method = '%s')",
                              batch, method))
      shiny::showNotification(
        if (method == "none") "Integration skipped (method = none)."
        else sprintf("Integration done using %s on '%s'.", method, batch),
        type = "message")
    })

    output$summary <- shiny::renderUI({
      if (!isTRUE(res$done)) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("Choose a batch column and method, then click Run integration.",
                               "选择批次列和方法，然后点击运行整合。")))
      }
      bslib::layout_columns(
        col_widths = c(4, 4, 4),
        stat_tile(i18n("Integrated", "已整合"),
                  if (res$method == "none") i18n("No", "否") else i18n("Yes", "是")),
        stat_tile(i18n("Method", "方法"), res$method),
        stat_tile(i18n("Batch column", "批次列"), res$batch)
      )
    })

    output$preview <- render_preview_plot(function() {
      shiny::req(res$done)
      obj <- rv$obj
      shiny::req(obj)
      batch <- res$batch
      reds <- obj_reductions(obj)
      # Prefer a 2D embedding to show batch mixing; fall back to PCA.
      reduction <- if ("umap" %in% reds) "umap" else if ("tsne" %in% reds) "tsne" else
        if ("pca" %in% reds) "pca" else NULL
      if (is.null(reduction) || is.null(batch) || !(batch %in% obj_meta_cols(obj))) {
        return(
          ggplot2::ggplot() +
            ggplot2::annotate("text", x = 0, y = 0,
                              label = "No embedding yet. / 尚无嵌入。\nRun an embedding step to preview batch mixing. / 运行嵌入步骤以预览批次混合。") +
            ggplot2::theme_void()
        )
      }
      df <- embedding_df(obj, reduction = reduction, color_by = batch)
      df$color <- factor(df$color)
      df$text <- sprintf("%s: %s\n%s: (%.2f, %.2f)",
                         batch, df$color, reduction, df$dim1, df$dim2)
      ggplot2::ggplot(df, ggplot2::aes(x = dim1, y = dim2,
                                       colour = color)) +
        ggplot2::geom_point(size = 0.5, alpha = 0.7) +
        ggplot2::scale_colour_manual(values = sc_palette(nlevels(df$color)),
                                     name = batch) +
        ggplot2::labs(x = paste0(reduction, " 1"), y = paste0(reduction, " 2"),
                      title = "Batch mixing after integration / 整合后的批次混合") +
        omicone_theme()
    })
  })
}
