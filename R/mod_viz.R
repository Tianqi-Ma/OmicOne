#' Module: Visualize
#'
#' A free exploration panel. Pick a plot type, a metadata column to group by, and
#' (for expression plots) some genes, then inspect the result and download it.
#' This module never modifies the working object -- it is read-only.
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_viz
NULL

#' @rdname mod_viz
#' @keywords internal
mod_viz_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = list(en = "Visualize", zh = "可视化"),
    what = list(
      en = "Explore your data freely: colour the embedding by any metadata column,
            or show the expression of specific genes.",
      zh = "自由探索数据：按任意元数据列为降维图着色，或展示特定基因的表达。"),
    why  = list(
      en = "A picture is the fastest way to sanity-check clustering, annotation and
            marker genes -- and to build the figures for your report.",
      zh = "作图是核查聚类、注释和标志基因是否合理的最快方式，也是为报告制作图表的手段。"),
    how  = list(
      en = "Pick a plot type. <b>Embedding by metadata</b> colours cells by a column
            on the UMAP (or, without one, the object's default reduction, named in
            the summary). <b>Violin / Dot / Feature / Heatmap</b> show expression of
            the genes you type. Nothing here changes your object.",
      zh = "选择一种图表类型。<b>按元数据着色的嵌入图</b>在 UMAP 上按某一列为细胞着色（若没有 UMAP，则使用对象的默认降维，名称显示在摘要中）。<b>小提琴图 / 点图 / 特征图 / 热图</b>展示你输入的基因的表达。此处的操作不会改动你的对象。"),
    read = list(
      en = "Pick a feature to colour cells by expression, or a grouping to
            compare populations. Use it to verify markers and annotations
            visually — the plot updates as you change the selection.",
      zh = "选择一个特征可按表达量为细胞着色，选择一个分组可比较细胞群。改选择时图会即时更新——用它直观验证标志基因与注释。"),
    example = list(
      en = "Type <code>CD3D, MS4A1, LYZ</code> and choose 'Feature plot' to see
               where T cells, B cells and monocytes sit on the UMAP.",
      zh = "输入 <code>CD3D, MS4A1, LYZ</code> 并选择“特征图”，即可查看 T 细胞、B 细胞和单核细胞在 UMAP 上的位置。")
  )
  controls <- shiny::tagList(
    label_with_help("Plot type",
                    "The embedding colours cells by metadata; the others show gene expression.",
                    "图表类型",
                    "嵌入图按元数据为细胞着色；其他类型展示基因表达。"),
    shiny::selectInput(ns("ptype"), NULL,
                       choices = c("Embedding by metadata" = "umap",
                                   "Violin plot"      = "violin",
                                   "Dot plot"         = "dotplot",
                                   "Feature plot"     = "feature",
                                   "Heatmap"          = "heatmap"),
                       selected = "umap"),
    label_with_help("Group by (metadata column)",
                    "Which metadata column to colour or split cells by (e.g. seurat_clusters, celltype).",
                    "分组依据（元数据列）",
                    "用于为细胞着色或分组的元数据列（例如 seurat_clusters、celltype）。"),
    shiny::selectInput(ns("meta_col"), NULL, choices = NULL),
    shiny::conditionalPanel(
      sprintf("input['%s'] != 'umap'", ns("ptype")),
      label_with_help("Genes",
                      "Gene names for expression plots, separated by commas, spaces or new lines.",
                      "基因",
                      "用于表达图的基因名，以逗号、空格或换行分隔。"),
      shiny::textInput(ns("genes"), NULL, placeholder = "e.g. CD3D, MS4A1, LYZ")
    ),
    shiny::checkboxInput(ns("mask"),
                         i18n("Outline cell types (mascarade)",
                              "细胞类型轮廓 (mascarade)"),
                         value = FALSE),
    label_with_help("Download format", "File type for the downloaded figure.",
                    "下载格式", "下载图片的文件类型。"),
    shiny::radioButtons(ns("fmt"), NULL, c("PNG" = "png", "PDF" = "pdf"), inline = TRUE),
    shiny::downloadButton(ns("download"),
                          i18n("Download plot", "下载图片"), class = "w-100")
  )
  step_container(
    title     = list(en = "Visualize", zh = "可视化"),
    subtitle  = list(en = "Explore any gene or grouping on the embedding.",
                     zh = "在嵌入图上自由查看任何基因或分组。"),
    explainer = explainer,
    controls  = controls,
    summary   = shiny::uiOutput(ns("summary")),
    preview   = shiny::uiOutput(ns("plot_slot"))
  )
}

#' @rdname mod_viz
#' @keywords internal
mod_viz_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    shiny::observe({
      obj <- rv$obj
      cols <- categorical_cols(obj_meta(obj))
      def <- default_group_col(cols, obj_misc(obj, "omicone_cluster_col"))
      shiny::updateSelectInput(session, "meta_col", choices = cols,
                               selected = keep_selected(shiny::isolate(input$meta_col), cols, def))
    })

    # The embedding drawn: UMAP when present, else the object's default
    # reduction (which can be a PCA, so its name is shown in the summary).
    reduction <- shiny::reactive({
      obj <- rv$obj
      if (has_reduction(obj, "umap")) return("umap")
      tryCatch(SeuratObject::DefaultDimReduc(obj), error = function(e) NULL)
    })

    # Is this plot type worth making interactive (ggplotly)?
    #
    # Only a plain ggplot can go through ggplotly: a multi-gene FeaturePlot comes
    # back as a patchwork, and the mascarade/scop dim plot is a composed object --
    # ggplotly errors on both, so those stay static.
    is_interactive <- shiny::reactive({
      if (!has_pkg("plotly") || !isTRUE(input$ptype %in% c("umap", "feature"))) {
        return(FALSE)
      }
      if (identical(input$ptype, "umap") && isTRUE(input$mask)) return(FALSE)
      if (identical(input$ptype, "feature") && length(parse_genes(input$genes)) != 1) {
        return(FALSE)
      }
      # Before anything is plotted (no data yet) fall back to the static output.
      p <- tryCatch(current_plot(), error = function(e) NULL)
      isTRUE(inherits(p, "ggplot") && !inherits(p, "patchwork"))
    })

    # Build the current plot as a ggplot object (or NULL on failure).
    current_plot <- shiny::reactive({
      obj <- rv$obj
      shiny::req(obj)
      if (!require_pkgs("Seurat", "Visualization")) return(NULL)
      genes <- parse_genes(input$genes)
      red <- reduction()
      tryCatch({
        switch(input$ptype,
          umap = {
            shiny::validate(shiny::need(length(red) && !is.na(red),
                                        "No embedding found. Run an embedding first."))
            # Outlined (mascarade) view: outlines follow the chosen column.
            if (isTRUE(input$mask) && has_pkg("scop") && !is.null(input$meta_col)) {
              p <- sc_dimplot(obj, group_by = input$meta_col, reduction = red,
                              mask = TRUE)
              if (!is.null(p)) return(p)
            }
            Seurat::DimPlot(obj, reduction = red, group.by = input$meta_col) +
              omicone_theme()
          },
          feature = {
            shiny::validate(shiny::need(length(genes) > 0, "Enter at least one gene."))
            shiny::validate(shiny::need(length(red) && !is.na(red),
                                        "No embedding found. Run an embedding first."))
            Seurat::FeaturePlot(obj, features = genes, reduction = red) &
              omicone_theme()
          },
          violin = {
            shiny::validate(shiny::need(length(genes) > 0, "Enter at least one gene."))
            Seurat::VlnPlot(obj, features = genes, group.by = input$meta_col) &
              omicone_theme()
          },
          dotplot = {
            shiny::validate(shiny::need(length(genes) > 0, "Enter at least one gene."))
            Seurat::DotPlot(obj, features = genes, group.by = input$meta_col) +
              omicone_theme()
          },
          heatmap = {
            shiny::validate(shiny::need(length(genes) > 0, "Enter at least one gene."))
            Seurat::DoHeatmap(obj, features = genes, group.by = input$meta_col)
          })
      },
      # req()/validate() are not errors: let them through so the output stays
      # blank (or shows the validate message) instead of raising a red toast.
      shiny.silent.error = function(e) stop(e),
      error = function(e) {
        shiny::showNotification(
          i18n(paste("Plot error:", conditionMessage(e)),
               paste("作图出错：", conditionMessage(e))),
          type = "error", duration = 10)
        NULL
      })
    })

    # Choose the correct output widget for the current plot type.
    output$plot_slot <- shiny::renderUI({
      if (is.null(rv$obj)) {
        return(explain_scene("viz",
                             "Load and process data, then pick a plot type on the left.",
                             "加载并处理数据后，在左侧选择图表类型。"))
      }
      if (is_interactive()) {
        plotly::plotlyOutput(ns("iplot"), height = "480px")
      } else {
        shiny::plotOutput(ns("splot"), height = "480px")
      }
    })

    # Only define the interactive output when plotly is installed; otherwise
    # referencing plotly:: at setup would error on machines without it.
    if (has_pkg("plotly")) {
      output$iplot <- plotly::renderPlotly({
        gg <- current_plot()
        shiny::req(gg)
        # "all" (not "text"): these are Seurat/ggplot layers with no `text`
        # aesthetic, so asking for "text" produced empty hover boxes.
        plotly::ggplotly(gg, tooltip = "all") |>
          plotly::config(displayModeBar = FALSE)
      })
    }

    output$splot <- shiny::renderPlot({
      gg <- current_plot()
      shiny::req(gg)
      gg
    })

    output$summary <- shiny::renderUI({
      obj <- rv$obj
      if (is.null(obj)) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("Load and process data first, then explore it here.",
                               "请先加载并处理数据，然后在此探索。")))
      }
      dims <- obj_dims(obj)
      genes <- parse_genes(input$genes)
      shiny::tagList(
        stat_tile(i18n("Cells", "细胞数"), format(dims$cells, big.mark = ",")),
        stat_tile(i18n("Embedding", "嵌入"), reduction() %||% i18n("none", "无")),
        stat_tile(i18n("Genes requested", "请求的基因数"), length(genes))
      )
    })

    output$download <- shiny::downloadHandler(
      filename = function() paste0("omicone_", input$ptype, "_",
                                   format(Sys.time(), "%Y%m%d_%H%M%S"), ".", input$fmt),
      content = function(file) {
        gg <- current_plot()
        shiny::req(gg)
        ggplot2::ggsave(file, plot = gg, width = 8, height = 6, dpi = 300,
                        device = input$fmt)
        mark_done(rv, "viz")
      }
    )
  })
}
