#' Module: Non-linear embedding (visualization)
#'
#' Compute a 2D embedding (UMAP, t-SNE, or PaCMAP) for visualizing the cellular
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
            gradients, and rare populations at a glance.",
      zh = "高维数据难以直接查看；二维图能一眼揭示细胞群、连续梯度和稀有细胞群体。"),
    how  = list(
      en = "UMAP is the default. Fewer neighbors / smaller min-dist emphasize
            local structure; larger values emphasize global layout. Base the
            embedding on the same reduction you clustered on.",
      zh = "默认使用 UMAP。更少的邻居数 / 更小的 min-dist 强调局部结构，更大的取值强调全局布局。降维应基于你聚类时所用的同一个线性降维结果。"),
    read = list(
      en = "Each point is a cell; distance approximates expression similarity.
            Islands are populations; thin bridges between them may be doublets
            or genuine transition states. Colours follow the current grouping.",
      zh = "每个点是一个细胞；距离近似表达相似性。岛屿是细胞群；岛屿之间细窄的桥可能是双细胞或真实的过渡状态。颜色跟随当前分组。"),
    example = list(
      en = "Distinct cell types appear as visually separated islands on the
               UMAP.<br><b>Note:</b> PaCMAP requires a Python backend and may not
               be available in every install.",
      zh = "不同的细胞类型会在 UMAP 上呈现为彼此分离的“岛屿”。<br><b>注意：</b>PaCMAP 需要 Python 后端，并非每个安装环境都可用。")
  )
  controls <- shiny::tagList(
    label_with_help("Method",
                    "UMAP is the default. t-SNE emphasizes local structure. PaCMAP needs Python.",
                    "方法",
                    "默认使用 UMAP。t-SNE 强调局部结构。PaCMAP 需要 Python。"),
    shiny::selectInput(ns("method"), NULL,
                       c("UMAP" = "umap", "t-SNE" = "tsne", "PaCMAP" = "pacmap")),
    label_with_help("Reduction", "Which reduction to embed (PCA or Harmony).",
                    "线性降维", "用于嵌入的线性降维结果（PCA 或 Harmony）。"),
    shiny::uiOutput(ns("reduction_ui")),
    label_with_help("Dimensions", "Number of leading dimensions to use.",
                    "维度数", "使用的前若干个维度的数量。"),
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
      label_with_help("perplexity", "Effective number of neighbors (t-SNE).",
                      "perplexity", "有效邻居数（t-SNE）。"),
      shiny::numericInput(ns("perplexity"), NULL, value = 30, min = 5, max = 100)
    ),
    shiny::checkboxInput(ns("mask"),
                         i18n("Outline cell types (mascarade)",
                              "细胞类型轮廓 (mascarade)"),
                         value = FALSE),
    run_button(ns("run"), "Run embedding", "运行降维")
  )
  step_container(title = list(en = "Embedding (UMAP / t-SNE)", zh = "降维可视化"),
                 subtitle = list(en = "A 2D view of the cell neighbourhood graph.",
                                 zh = "细胞邻域图的二维视图。"),
                 explainer = explainer, controls = controls,
                 summary = shiny::uiOutput(ns("summary")),
                 preview = preview_plot_ui(ns("preview"),
                   guide = list(en = "The embedding will be drawn here.",
                                zh = "运行后，这里将绘制降维嵌入图。"),
                   caption = list(en = "2D embedding coloured by the active grouping.",
                                  zh = "按当前分组着色的二维嵌入图。")))
}

#' @rdname mod_embed
#' @keywords internal
mod_embed_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    res <- shiny::reactiveValues(done = FALSE, method = NULL, reduction = NULL,
                                 params = NULL)

    # Offer the linear reductions the embedding can be built on.
    output$reduction_ui <- shiny::renderUI({
      reds <- obj_reductions(rv$obj)
      choices <- intersect(c("pca", "harmony"), reds)
      if (length(choices) == 0) choices <- reds
      if (length(choices) == 0) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("No reductions yet — run PCA first.",
                               "还没有降维结果 —— 请先运行 PCA。")))
      }
      selected <- if ("harmony" %in% choices) "harmony" else choices[1]
      shiny::selectInput(session$ns("reduction"), NULL,
                         choices = choices, selected = selected)
    })

    shiny::observeEvent(input$run, {
      shiny::req(rv$obj)
      shiny::req(input$reduction)
      if (!require_pkgs("Seurat", "Embedding")) return(NULL)
      method    <- input$method
      reduction <- input$reduction
      dims      <- input$dims
      n_neighbors <- input$n_neighbors
      min_dist    <- input$min_dist
      perplexity  <- input$perplexity
      obj <- with_progress_notify({
        embed_obj(rv$obj, method = method, reduction = reduction, dims = dims,
                  n_neighbors = n_neighbors, min_dist = min_dist,
                  perplexity = perplexity)
      }, message = sprintf("Computing %s...", toupper(method)))
      if (is.null(obj)) return(NULL)
      rv$obj <- obj
      res$done      <- TRUE
      res$method    <- method
      res$reduction <- reduction
      res$params <- if (method == "umap") {
        sprintf("n_neighbors=%s, min_dist=%s", n_neighbors, min_dist)
      } else if (method == "tsne") {
        sprintf("perplexity=%s", perplexity)
      } else {
        sprintf("dims=%s", dims)
      }
      mark_done(rv, "embed")
      log_step(log_rv, "Embedding",
               params = list(method = method, reduction = reduction, dims = dims,
                             n_neighbors = n_neighbors, min_dist = min_dist,
                             perplexity = perplexity),
               code = sprintf(
                 "obj <- embed_obj(obj, method = '%s', reduction = '%s', dims = %d)",
                 method, reduction, dims))
      shiny::showNotification(sprintf("%s embedding done.", toupper(method)),
                              type = "message")
    })

    output$summary <- shiny::renderUI({
      if (!isTRUE(res$done)) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("Set parameters and click Run embedding.",
                               "设置参数后点击“运行降维”。")))
      }
      bslib::layout_columns(
        col_widths = c(4, 4, 4),
        stat_tile(i18n("Method", "方法"), toupper(res$method)),
        stat_tile(i18n("Reduction", "线性降维"), res$reduction),
        stat_tile(i18n("Parameters", "参数"), res$params)
      )
    })

    output$preview <- render_preview_plot(function() {
      shiny::req(res$done)
      obj <- rv$obj
      shiny::req(obj)
      # The embedding is stored under the method name (e.g. "umap", "tsne").
      reduction <- res$method
      if (!has_reduction(obj, reduction)) {
        return(
          ggplot2::ggplot() +
            ggplot2::annotate("text", x = 0, y = 0,
                              label = paste0("Embedding '", reduction, "' not available.")) +
            ggplot2::theme_void()
        )
      }
      # Colour by a cluster column if one exists.
      cols <- obj_meta_cols(obj)
      color_by <- intersect(c("seurat_clusters"), cols)
      color_by <- if (length(color_by)) color_by[1] else NULL
      # Outlined (mascarade) view: use scop's dim plot with mask overlays.
      if (isTRUE(input$mask) && has_pkg("scop") && !is.null(color_by)) {
        p <- tryCatch(sc_dimplot(obj, group_by = color_by, reduction = reduction,
                                 mask = TRUE),
                      error = function(e) NULL)
        if (!is.null(p)) return(p)
      }
      df <- embedding_df(obj, reduction = reduction, color_by = color_by)
      if (!is.null(color_by)) {
        df$color <- factor(df$color)
        df$text <- sprintf("cluster: %s\n%s: (%.2f, %.2f)",
                           df$color, reduction, df$dim1, df$dim2)
        p <- ggplot2::ggplot(df, ggplot2::aes(x = dim1, y = dim2,
                                              colour = color)) +
          ggplot2::geom_point(size = 0.5, alpha = 0.7) +
          ggplot2::scale_colour_manual(values = sc_palette(nlevels(df$color)),
                                       name = "cluster")
      } else {
        df$text <- sprintf("%s: (%.2f, %.2f)", reduction, df$dim1, df$dim2)
        p <- ggplot2::ggplot(df, ggplot2::aes(x = dim1, y = dim2)) +
          ggplot2::geom_point(size = 0.5, alpha = 0.7, colour = "#3b6ea5")
      }
      p +
        ggplot2::labs(x = paste0(toupper(reduction), " 1"),
                      y = paste0(toupper(reduction), " 2"),
                      title = sprintf("%s embedding", toupper(res$method))) +
        omicone_theme()
    })
  })
}
