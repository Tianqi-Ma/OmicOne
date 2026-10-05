#' Module: Non-linear embedding (visualization)
#'
#' Compute a 2D embedding (UMAP or t-SNE) for visualizing the cellular
#' structure. UMAP is the default. The embedding is for display only; clustering
#' and statistics use the linear reduction (PCA/Harmony), not the 2D map.
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_embed
NULL

#' @rdname mod_embed
#' @keywords internal
mod_embed_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = list(en = "Non-linear embedding (visualization)",
                 zh = "非线性降维（可视化）"),
    what = list(
      en = "Project cells into 2D so you can see clusters and structure.",
      zh = "将细胞投影到二维，以便直观查看簇和结构。"),
    why  = list(
      en = "High-dimensional data is hard to inspect; a 2D map reveals groups,
            gradients, and rare populations at a glance. Distances between
            islands on a UMAP are not quantitative.",
      zh = "高维数据难以直接查看；二维图能一眼揭示细胞群、连续梯度和稀有细胞群体。UMAP 上岛屿之间的距离不具定量意义。"),
    how  = list(
      en = "UMAP is the default. Fewer neighbors / smaller min-dist emphasize
            local structure; larger values emphasize global layout. Base the
            embedding on the same reduction you clustered on (selected by
            default after integration).",
      zh = "默认使用 UMAP。更少的邻居数 / 更小的 min-dist 强调局部结构，更大的取值强调全局布局。降维应基于你聚类时所用的同一个线性降维结果（整合后默认选中）。"),
    read = list(
      en = "Each point is a cell; nearby cells have similar expression. Islands
            are populations; thin bridges between them may be doublets or
            genuine transition states. Colours follow the grouping chosen under
            Colour by.",
      zh = "每个点是一个细胞；相邻细胞的表达相似。岛屿是细胞群；岛屿之间细窄的桥可能是双细胞或真实的过渡状态。颜色按“着色依据”中选择的分组显示。"),
    example = list(
      en = "Distinct cell types appear as visually separated islands on the
               UMAP; colour by sample to check that batches overlap after
               integration.",
      zh = "不同的细胞类型会在 UMAP 上呈现为彼此分离的“岛屿”；按样本着色可检查整合后各批次是否重叠。")
  )
  controls <- shiny::tagList(
    label_with_help("Method",
                    "UMAP is the default. t-SNE emphasizes local structure.",
                    "方法",
                    "默认使用 UMAP。t-SNE 强调局部结构。"),
    shiny::selectInput(ns("method"), NULL, c("UMAP" = "umap", "t-SNE" = "tsne")),
    label_with_help("Reduction", "Which reduction to embed (PCA, Harmony or Seurat's integrated.dr).",
                    "线性降维", "用于嵌入的线性降维结果（PCA、Harmony 或 Seurat 的 integrated.dr）。"),
    shiny::uiOutput(ns("reduction_ui")),
    label_with_help("Dimensions", "Number of leading dimensions to use (capped at those available).",
                    "维度数", "使用的前若干个维度的数量（不超过已有维度）。"),
    shiny::numericInput(ns("dims"), NULL, value = 30, min = 2, max = 100),
    shiny::conditionalPanel(
      sprintf("input['%s'] == 'umap'", ns("method")),
      label_with_help("n_neighbors", "Balances local vs global structure (UMAP).",
                      "n_neighbors", "在局部结构与全局结构之间取得平衡（UMAP）。"),
      shiny::numericInput(ns("n_neighbors"), NULL, value = 30, min = 2, max = 200),
      label_with_help("min_dist", "Minimum spacing between points (UMAP).",
                      "min_dist", "点与点之间的最小间距（UMAP）。"),
      shiny::numericInput(ns("min_dist"), NULL, value = 0.3, min = 0, max = 1, step = 0.05)
    ),
    shiny::conditionalPanel(
      sprintf("input['%s'] == 'tsne'", ns("method")),
      label_with_help("perplexity", "Effective number of neighbors (t-SNE); capped below (cells - 1) / 3.",
                      "perplexity", "有效邻居数（t-SNE）；上限为 (细胞数 - 1) / 3。"),
      shiny::numericInput(ns("perplexity"), NULL, value = 30, min = 5, max = 100)
    ),
    label_with_help("Colour by", "Metadata column used to colour the map (display only; changing it does not re-run anything).",
                    "着色依据", "用于给图着色的元数据列（仅影响显示，修改不会重新计算）。"),
    shiny::uiOutput(ns("color_ui")),
    shiny::checkboxInput(ns("mask"),
                         i18n("Outline groups (mascarade)",
                              "分组轮廓 (mascarade)"),
                         value = FALSE),
    run_button(ns("run"), "Run embedding", "运行降维")
  )
  step_container(id = id, title = list(en = "Embedding (UMAP / t-SNE)", zh = "降维可视化"),
                 subtitle = list(en = "A 2D view of the cell neighbourhood graph.",
                                 zh = "细胞邻域图的二维视图。"),
                 explainer = explainer, controls = controls,
                 summary = shiny::uiOutput(ns("summary")),
                 preview = preview_plot_ui(ns("preview"), download = TRUE,
                   guide = list(en = "The embedding will be drawn here.",
                                zh = "运行后，这里将绘制降维嵌入图。"),
                   caption = list(en = "One point = one cell, coloured by the column chosen under Colour by.",
                                  zh = "每个点为一个细胞，按“着色依据”中选择的列着色。")))
}

#' @rdname mod_embed
#' @keywords internal
mod_embed_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    res <- step_results(rv, "sc", done = FALSE, method = NULL, reduction = NULL,
                        params = NULL)
    last_default <- NULL

    output$reduction_ui <- shiny::renderUI({
      choices <- graph_reduction_choices(rv$obj)
      if (length(choices) == 0) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("No reductions yet — run PCA first.",
                               "还没有降维结果 —— 请先运行 PCA。")))
      }
      def <- default_graph_reduction(rv$obj)
      sel <- reduction_selection(shiny::isolate(input$reduction), choices, def, last_default)
      last_default <<- def
      shiny::selectInput(session$ns("reduction"), NULL, choices = choices, selected = sel)
    })

    output$color_ui <- shiny::renderUI({
      md <- obj_meta(rv$obj)
      cols <- sc_group_cols(md, min_levels = 1)
      if (!length(cols)) {
        return(shiny::div(class = "omicone-note",
                          i18n("No categorical metadata yet.", "暂无分类元数据。")))
      }
      active <- obj_misc(rv$obj, "omicone_cluster_col")
      default <- intersect(c("celltype", active, "seurat_clusters"), cols)
      default <- if (length(default)) default[1] else cols[1]
      shiny::selectInput(session$ns("color_by"), NULL, choices = cols,
                         selected = keep_selected(shiny::isolate(input$color_by), cols, default))
    })

    shiny::observeEvent(input$run, {
      shiny::req(rv$obj)
      shiny::req(input$reduction)
      if (!require_pkgs("Seurat", "Embedding")) return(NULL)
      method    <- input$method
      reduction <- input$reduction
      dims      <- int_input(input$dims, 2, 500)
      n_neighbors <- int_input(input$n_neighbors, 2, 1000)
      min_dist    <- num_input(input$min_dist, 0, 1)
      perplexity  <- num_input(input$perplexity, 1, 1000)
      if (is.na(dims) || (method == "umap" && (is.na(n_neighbors) || is.na(min_dist))) ||
          (method == "tsne" && is.na(perplexity))) {
        shiny::showNotification(i18n("Enter valid numbers for the embedding parameters.",
                                     "请为降维参数输入有效数值。"), type = "error")
        return(NULL)
      }
      if (!reduction %in% graph_reduction_choices(rv$obj)) return(NULL)
      pr <- embed_params(ncol(rv$obj),
                         ncol(SeuratObject::Embeddings(rv$obj, reduction = reduction)),
                         dims, if (is.na(n_neighbors)) 30L else n_neighbors,
                         if (is.na(perplexity)) 30 else perplexity)
      obj <- with_progress_notify({
        embed_obj(rv$obj, method = method, reduction = reduction, dims = pr$dims,
                  n_neighbors = pr$n_neighbors, min_dist = min_dist,
                  perplexity = pr$perplexity)
      }, message = sprintf("Computing %s...", toupper(method)))
      if (is.null(obj)) return(NULL)
      rv$obj <- obj
      res$done      <- TRUE
      res$method    <- method
      res$reduction <- reduction
      res$params <- if (method == "umap") {
        sprintf("%d dims, n_neighbors = %d, min_dist = %s", pr$dims, pr$n_neighbors, min_dist)
      } else {
        sprintf("%d dims, perplexity = %s", pr$dims, pr$perplexity)
      }
      mark_done(rv, "embed")
      log_step(log_rv, "Embedding",
               params = c(list(method = method, reduction = reduction, dims = pr$dims),
                          if (method == "umap") list(n_neighbors = pr$n_neighbors,
                                                     min_dist = min_dist, seed = 42)
                          else list(perplexity = pr$perplexity, seed = 1)),
               code = embed_log_code(method, reduction, pr, min_dist))
      shiny::showNotification(i18n(sprintf("%s embedding done.", toupper(method)),
                                   sprintf("%s 降维完成。", toupper(method))),
                              type = "message")
    })

    output$summary <- shiny::renderUI({
      if (!isTRUE(res$done)) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("Set parameters and click Run embedding.",
                               "设置参数后点击“运行降维”。")))
      }
      shiny::tagList(
        stat_tile(i18n("Method", "方法"), toupper(res$method)),
        stat_tile(i18n("Reduction", "线性降维"), res$reduction),
        stat_tile(i18n("Parameters", "参数"), res$params)
      )
    })

    render_step_plot(output, input, "preview", function() {
      shiny::req(res$done)
      obj <- rv$obj
      shiny::req(obj, has_reduction(obj, res$method))
      color_by <- input$color_by
      if (!is.null(color_by) && !color_by %in% obj_meta_cols(obj)) color_by <- NULL
      if (isTRUE(input$mask) && has_pkg("scop") && !is.null(color_by)) {
        p <- tryCatch(sc_dimplot(obj, group_by = color_by, reduction = res$method,
                                 mask = TRUE),
                      error = function(e) NULL)
        if (!is.null(p)) return(p)
      }
      embed_plot(obj, res$method, color_by)
    }, name = "embed")
  })
}
