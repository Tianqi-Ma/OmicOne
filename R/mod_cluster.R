#' Module: Clustering
#'
#' Group cells into clusters (candidate cell populations) with a
#' neighbor-graph + community-detection approach, at one or more resolutions.
#' Leiden (algorithm 4) is the modern default; Louvain (algorithm 1) is classic.
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_cluster
NULL

#' @rdname mod_cluster
#' @keywords internal
mod_cluster_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = list(en = "Clustering", zh = "聚类"),
    what = list(
      en = "Partition cells into clusters that likely correspond to distinct
            cell types or states.",
      zh = "将细胞划分为可能对应不同细胞类型或状态的簇。"),
    why  = list(
      en = "Clusters are the units you annotate and compare. Good clustering
            separates real populations without over-splitting noise.",
      zh = "簇是你用来注释和比较的单位。好的聚类能分开真实的细胞群，而不会把噪声过度切分。"),
    how  = list(
      en = "Higher <b>resolution</b> = more, smaller clusters. Try several
            resolutions and compare. Use the reduction you want to cluster on
            (PCA, or an integrated embedding like Harmony).",
      zh = "<b>分辨率</b>越高 = 簇越多、越小。尝试多个分辨率并比较。使用你想在其上聚类的
            降维（PCA，或像 Harmony 这样的整合嵌入）。"),
    read = list(
      en = "Each colour is one cluster. Resolution controls granularity — more
            clusters split subtypes, fewer merge them. If one cluster swallows
            most cells, the resolution is too low.",
      zh = "每种颜色是一个簇。分辨率控制粒度——簇多则拆分亚型，簇少则合并。若一个簇吞掉大多数细胞，说明分辨率太低。"),
    example = list(
      en = "At resolution 0.2 you may get 6 broad clusters; at 1.0 they split
               into finer subtypes.<br><b>Note:</b> Leiden (algorithm 4) needs the
               Python <code>leidenalg</code> (or <code>leidenbase</code>) backend;
               if unavailable, switch to Louvain.",
      zh = "在分辨率 0.2 时你可能得到 6 个大簇；在 1.0 时它们会分裂为更细的亚型。
               <br><b>注意：</b>Leiden（算法 4）需要 Python 的 <code>leidenalg</code>
               （或 <code>leidenbase</code>）后端；若不可用，请切换到 Louvain。")
  )
  controls <- shiny::tagList(
    label_with_help("Method",
                    "Leiden (algorithm 4) is the modern default; Louvain (algorithm 1) is classic.",
                    "方法",
                    "Leiden（算法 4）是现代默认方法；Louvain（算法 1）是经典方法。"),
    shiny::selectInput(ns("method"), NULL,
                       c("Leiden (algorithm 4)" = "leiden",
                         "Louvain (algorithm 1)" = "louvain")),
    label_with_help("Reduction",
                    "Which dimensional reduction to build the neighbor graph on.",
                    "降维",
                    "在哪个降维结果上构建近邻图。"),
    shiny::uiOutput(ns("reduction_ui")),
    label_with_help("Dimensions", "Number of leading dimensions to use (e.g. PCs).",
                    "维度", "使用的前若干个维度的数量（例如主成分）。"),
    shiny::numericInput(ns("dims"), NULL, value = 30, min = 2, max = 100),
    label_with_help("Neighbors (k)", "Neighbors used to build the graph.",
                    "近邻数 (k)", "用于构建图的近邻数量。"),
    shiny::numericInput(ns("neighbors"), NULL, value = 20, min = 2, max = 100),
    label_with_help("Resolutions",
                    "Comma-separated list; each is clustered. The last one is used for the summary/preview.",
                    "分辨率",
                    "逗号分隔的列表；每个都会进行聚类。最后一个用于摘要/预览。"),
    shiny::textInput(ns("resolutions"), NULL, value = "0.2,0.5,0.8,1.0"),
    run_button(ns("run"), "Run clustering", "运行聚类")
  )
  step_container(title = list(en = "Clustering", zh = "聚类"),
                 subtitle = list(en = "Group cells by expression similarity.",
                                 zh = "按表达相似性对细胞分群。"),
                 explainer = explainer, controls = controls,
                 summary = shiny::uiOutput(ns("summary")),
                 preview = preview_plot_ui(ns("preview"),
                   guide = list(en = "The clustering result will be drawn here.",
                                zh = "运行后，这里将绘制聚类结果。"),
                   caption = list(en = "Cells coloured by cluster at the chosen resolution.",
                                  zh = "按所选分辨率对细胞按簇着色。")))
}

#' @rdname mod_cluster
#' @keywords internal
mod_cluster_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    res <- shiny::reactiveValues(done = FALSE, col = NULL, res_used = NULL,
                                 n_clusters = NA_integer_)

    # Offer reductions suitable for clustering (pca / harmony if present).
    output$reduction_ui <- shiny::renderUI({
      reds <- obj_reductions(rv$obj)
      choices <- intersect(c("pca", "harmony"), reds)
      if (length(choices) == 0) choices <- reds
      if (length(choices) == 0) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("No reductions yet — run PCA first.",
                               "尚无降维结果——请先运行 PCA。")))
      }
      selected <- if ("harmony" %in% choices) "harmony" else choices[1]
      shiny::selectInput(session$ns("reduction"), NULL,
                         choices = choices, selected = selected)
    })

    shiny::observeEvent(input$run, {
      shiny::req(rv$obj)
      shiny::req(input$reduction)
      if (!require_pkgs("Seurat", "Clustering")) return(NULL)
      # Parse the comma-separated resolutions into a numeric vector.
      resolutions <- suppressWarnings(as.numeric(
        trimws(strsplit(input$resolutions, ",", fixed = TRUE)[[1]])))
      resolutions <- resolutions[!is.na(resolutions)]
      if (length(resolutions) == 0) {
        shiny::showNotification("Enter at least one valid resolution.",
                                type = "error")
        return(NULL)
      }
      algorithm <- if (input$method == "leiden") 4 else 1
      reduction <- input$reduction
      dims <- input$dims
      obj <- with_progress_notify({
        cluster_obj(rv$obj, reduction = reduction, dims = dims,
                    resolutions = resolutions, algorithm = algorithm)
      }, message = "Building graph and clustering...")
      if (is.null(obj)) return(NULL)
      rv$obj <- obj
      # Seurat stores the last FindClusters result in `seurat_clusters`.
      last_res <- resolutions[length(resolutions)]
      clusters <- obj_meta(obj)$seurat_clusters
      res$done       <- TRUE
      res$col        <- "seurat_clusters"
      res$res_used   <- last_res
      res$n_clusters <- if (is.null(clusters)) NA_integer_ else nlevels(factor(clusters))
      mark_done(rv, "cluster")
      log_step(log_rv, "Clustering",
               params = list(method = input$method, algorithm = algorithm,
                             reduction = reduction, dims = dims,
                             resolutions = resolutions),
               code = sprintf(
                 "obj <- cluster_obj(obj, reduction = '%s', dims = %d, resolutions = c(%s), algorithm = %d)",
                 reduction, dims, paste(resolutions, collapse = ", "), algorithm))
      shiny::showNotification(
        sprintf("Clustering done: %d clusters at resolution %s.",
                res$n_clusters, format(last_res)),
        type = "message")
    })

    output$summary <- shiny::renderUI({
      if (!isTRUE(res$done)) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("Set parameters and click Run clustering.",
                               "设置参数并点击运行聚类。")))
      }
      bslib::layout_columns(
        col_widths = c(6, 6),
        stat_tile(i18n("Clusters", "簇数"), format(res$n_clusters)),
        stat_tile(i18n("Resolution", "分辨率"), format(res$res_used))
      )
    })

    output$preview <- render_preview_plot(function() {
      shiny::req(res$done)
      md <- obj_meta(rv$obj)
      shiny::req(res$col %in% colnames(md))
      cl <- factor(md[[res$col]])

      # If a 2D embedding already exists, show clusters on the map (most useful).
      red <- if (has_reduction(rv$obj, "umap")) "umap"
             else if (has_reduction(rv$obj, "tsne")) "tsne" else NULL
      if (!is.null(red)) {
        df <- embedding_df(rv$obj, red, color_by = res$col)
        df$text <- sprintf("Cluster %s", df$color)
        return(
          ggplot2::ggplot(df, ggplot2::aes(dim1, dim2, colour = factor(color))) +
            ggplot2::geom_point(size = 0.6, alpha = 0.75) +
            ggplot2::scale_colour_manual(values = sc_palette(nlevels(cl)), name = "Cluster / 簇") +
            ggplot2::labs(x = paste0(toupper(red), " 1"), y = paste0(toupper(red), " 2"),
                          title = sprintf("Clusters on %s (resolution %s) / %s 上的聚类（分辨率 %s）",
                                          toupper(red), format(res$res_used),
                                          toupper(red), format(res$res_used))) +
            omicone_theme()
        )
      }

      # No embedding yet: show cluster sizes and point users to the Embed step.
      counts <- as.data.frame(table(cluster = cl), stringsAsFactors = FALSE)
      names(counts) <- c("cluster", "n")
      counts$cluster <- factor(counts$cluster, levels = levels(cl))
      counts$text <- sprintf("Cluster %s\n%s cells",
                             counts$cluster, format(counts$n, big.mark = ","))
      ggplot2::ggplot(counts, ggplot2::aes(x = cluster, y = n,
                                           fill = cluster)) +
        ggplot2::geom_col() +
        ggplot2::scale_fill_manual(values = sc_palette(nlevels(counts$cluster)),
                                   guide = "none") +
        ggplot2::labs(x = "Cluster", y = "Cells",
                      title = sprintf("Cluster sizes (resolution %s) / 簇大小（分辨率 %s）",
                                      format(res$res_used), format(res$res_used)),
                      caption = "No 2D map yet - run step 8 (Embed) to see clusters on a UMAP. / 尚无二维图——运行第 8 步（嵌入）以在 UMAP 上查看聚类。") +
        omicone_theme()
    })
  })
}
