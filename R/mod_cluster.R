#' Module: Clustering
#'
#' Group cells into clusters (candidate cell populations) with a
#' neighbor-graph + community-detection approach, at one or more resolutions.
#' Leiden (algorithm 4) is the modern default; Louvain (algorithm 1) is classic.
#' One resolution is the active clustering (Idents / `seurat_clusters`), chosen
#' explicitly by the user.
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
      en = "Higher <b>resolution</b> = more, smaller clusters; several can be run
            at once, and the <b>active resolution</b> chosen after the run is
            what markers, annotation and later steps use. Cluster on the
            reduction the Integrate step produced (Harmony / integrated.dr) when
            you have several batches. <b>k</b> is the number of neighbours in
            the graph (smaller = finer local structure).",
      zh = "<b>分辨率</b>越高 = 簇越多、越小；可一次运行多个分辨率，运行后选定的<b>当前分辨率</b>就是标志基因、注释及后续步骤使用的聚类。多批次数据请在整合步骤产生的降维（Harmony / integrated.dr）上聚类。<b>k</b> 为图中的近邻数（越小越能体现局部细节）。"),
    read = list(
      en = "Each colour is one cluster of the active resolution, on the UMAP if
            one exists, otherwise as cluster sizes. If one cluster swallows most
            cells, the resolution is too low; many tiny clusters suggest it is
            too high.",
      zh = "每种颜色是当前分辨率下的一个簇；若已有 UMAP 则画在 UMAP 上，否则显示各簇大小。若一个簇吞掉大多数细胞，说明分辨率太低；大量很小的簇则提示分辨率过高。"),
    example = list(
      en = "At resolution 0.2 you may get 6 broad clusters; at 1.0 they split
               into finer subtypes.<br><b>Note:</b> Leiden needs the R package
               <code>leidenbase</code> with Seurat 5.2 or later (Seurat versions that
               offer it fall back to igraph); otherwise choose Louvain.",
      zh = "在分辨率 0.2 时你可能得到 6 个大簇；在 1.0 时它们会分裂为更细的亚型。<br><b>注意：</b>Leiden 需要 Seurat 5.2 及以上版本配合 R 包 <code>leidenbase</code>（支持的 Seurat 版本可回退到 igraph）；否则请选择 Louvain。")
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
                    "Which dimensional reduction to build the neighbor graph on. Defaults to the space the Integrate step produced.",
                    "降维",
                    "在哪个降维结果上构建近邻图。默认使用整合步骤产生的空间。"),
    shiny::uiOutput(ns("reduction_ui")),
    label_with_help("Dimensions", "Number of leading dimensions to use (capped at those available).",
                    "维度", "使用的前若干个维度的数量（不超过已有维度）。"),
    shiny::numericInput(ns("dims"), NULL, value = 30, min = 2, max = 100),
    label_with_help("Neighbors (k)", "Neighbors used to build the shared-nearest-neighbour graph (k.param, Seurat default 20).",
                    "近邻数 (k)", "构建共享近邻图所用的近邻数（k.param，Seurat 默认 20）。"),
    shiny::numericInput(ns("neighbors"), NULL, value = 20, min = 2, max = 100),
    label_with_help("Resolutions",
                    "Comma-separated list; each is clustered and kept as a column. The first becomes the active clustering; switch it after the run.",
                    "分辨率",
                    "逗号分隔的列表；每个都会聚类并保存为一列。第一个成为当前聚类；运行后可切换。"),
    shiny::textInput(ns("resolutions"), NULL, value = "0.5,0.2,0.8,1.0"),
    run_button(ns("run"), "Run clustering", "运行聚类"),
    shiny::uiOutput(ns("active_ui"))
  )
  step_container(title = list(en = "Clustering", zh = "聚类"),
                 subtitle = list(en = "Group cells by expression similarity, then pick the active resolution.",
                                 zh = "按表达相似性对细胞分群，再选定当前分辨率。"),
                 explainer = explainer, controls = controls,
                 summary = shiny::uiOutput(ns("summary")),
                 preview = shiny::tagList(
                   shiny::uiOutput(ns("insight")),
                   preview_plot_ui(ns("preview"), download = TRUE,
                     guide = list(en = "The active clustering will be drawn here.",
                                  zh = "运行后，这里将绘制当前聚类结果。"),
                     caption = list(en = "Cells coloured by cluster at the active resolution (or cluster sizes before an embedding exists).",
                                    zh = "按当前分辨率的簇为细胞着色（尚无嵌入图时显示各簇大小）。"))))
}

#' @rdname mod_cluster
#' @keywords internal
mod_cluster_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    res <- step_results(rv, "sc", done = FALSE, cols = NULL, active = NULL,
                        n_clusters = NA_integer_, params = NULL)
    last_default <- NULL

    output$reduction_ui <- shiny::renderUI({
      choices <- graph_reduction_choices(rv$obj)
      if (length(choices) == 0) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("No reductions yet — run PCA first.",
                               "尚无降维结果——请先运行 PCA。")))
      }
      def <- default_graph_reduction(rv$obj)
      sel <- reduction_selection(shiny::isolate(input$reduction), choices, def, last_default)
      last_default <<- def
      shiny::selectInput(session$ns("reduction"), NULL, choices = choices, selected = sel)
    })

    # The active resolution: written to Idents / seurat_clusters on change.
    output$active_ui <- shiny::renderUI({
      cols <- res$cols
      if (is.null(cols)) return(NULL)
      shiny::tagList(
        label_with_help("Active resolution",
                        "The clustering used by markers, annotation and every later step (written to Idents and seurat_clusters).",
                        "当前分辨率",
                        "标志基因、注释及后续所有步骤使用的聚类（写入 Idents 和 seurat_clusters）。"),
        shiny::selectInput(session$ns("active"), NULL,
                           choices = stats::setNames(unname(cols), paste("resolution", names(cols))),
                           selected = res$active)
      )
    })

    write_log <- function(p, active) {
      log_step(log_rv, "Clustering",
               params = list(method = p$method, algorithm = p$algorithm,
                             reduction = p$reduction, dims = p$dims, k = p$k,
                             resolutions = p$resolutions, seed = p$seed,
                             active = active),
               code = cluster_log_code(p$reduction, p$dims, p$resolutions, p$algorithm,
                                       p$k, p$seed, p$assay, active, p$leiden_igraph))
    }

    shiny::observeEvent(input$run, {
      shiny::req(rv$obj)
      shiny::req(input$reduction)
      if (!require_pkgs("Seurat", "Clustering")) return(NULL)
      resolutions <- suppressWarnings(as.numeric(
        trimws(strsplit(input$resolutions %||% "", ",", fixed = TRUE)[[1]])))
      resolutions <- unique(resolutions[is.finite(resolutions) & resolutions > 0])
      dims <- int_input(input$dims, 2, 500)
      k <- int_input(input$neighbors, 2, 500)
      if (length(resolutions) == 0 || is.na(dims) || is.na(k)) {
        shiny::showNotification(i18n("Enter whole numbers for dimensions and k, and at least one positive resolution.",
                                     "请为维度和 k 输入整数，并至少输入一个正的分辨率。"),
                                type = "error")
        return(NULL)
      }
      method <- input$method
      algorithm <- if (method == "leiden") 4 else 1
      reduction <- input$reduction
      if (!reduction %in% graph_reduction_choices(rv$obj)) return(NULL)
      seed <- 0
      dims <- min(dims, ncol(SeuratObject::Embeddings(rv$obj, reduction = reduction)))
      leiden_igraph <- algorithm == 4 && identical(leiden_backend(), "igraph")
      obj <- with_progress_notify({
        cluster_obj(rv$obj, reduction = reduction, dims = dims,
                    resolutions = resolutions, algorithm = algorithm, k = k, seed = seed)
      }, message = "Building graph and clustering...")
      if (is.null(obj)) return(NULL)
      rv$obj <- obj
      cols <- obj_misc(obj, "omicone_cluster_cols")
      active <- obj_misc(obj, "omicone_cluster_col")
      res$params <- list(method = method, algorithm = algorithm, reduction = reduction,
                         dims = dims, k = k, resolutions = resolutions, seed = seed,
                         assay = obj_default_assay(obj) %||% "RNA",
                         leiden_igraph = leiden_igraph)
      res$done <- TRUE
      res$cols <- cols
      res$active <- active
      res$n_clusters <- nlevels(factor(obj_meta(obj)[[active]]))
      mark_done(rv, "cluster")
      write_log(res$params, active)
      shiny::showNotification(
        i18n(sprintf("Clustering done: %d clusters at resolution %s (active).",
                     res$n_clusters, names(cols)[cols == active]),
             sprintf("聚类完成：当前分辨率 %s 下有 %d 个簇。",
                     names(cols)[cols == active], res$n_clusters)),
        type = "message")
    })

    shiny::observeEvent(input$active, {
      col <- input$active
      shiny::req(rv$obj, res$cols, col)
      if (identical(col, obj_misc(rv$obj, "omicone_cluster_col"))) return(NULL)
      if (!col %in% obj_meta_cols(rv$obj)) return(NULL)
      obj <- with_progress_notify(cluster_set_active(rv$obj, col),
                                  message = "Switching the active clustering...")
      if (is.null(obj)) return(NULL)
      rv$obj <- obj
      res$active <- col
      res$n_clusters <- nlevels(factor(obj_meta(obj)[[col]]))
      mark_done(rv, "cluster")
      write_log(res$params, col)
    }, ignoreInit = TRUE)

    output$summary <- shiny::renderUI({
      if (!isTRUE(res$done)) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("Set parameters and click Run clustering.",
                               "设置参数并点击运行聚类。")))
      }
      cols <- res$cols
      shiny::tagList(
        stat_tile(i18n("Clusters", "簇数"), format(res$n_clusters)),
        stat_tile(i18n("Active resolution", "当前分辨率"), names(cols)[cols == res$active]),
        stat_tile(i18n("Graph", "近邻图"),
                  sprintf("%s, %d dims, k = %d", res$params$reduction, res$params$dims,
                          res$params$k))
      )
    })

    output$insight <- shiny::renderUI({
      if (!isTRUE(res$done)) return(NULL)
      md <- obj_meta(rv$obj)
      cols <- res$cols
      cols <- cols[cols %in% names(md)]
      if (!length(cols)) return(NULL)
      n_each <- vapply(cols, function(cl) nlevels(factor(md[[cl]])), integer(1))
      tab <- table(md[[res$active]])
      small <- sum(tab < 10)
      insight_bar(
        sprintf("Clusters per resolution: %s. Active: resolution %s, %d clusters (largest %.0f%% of cells%s).",
                paste(sprintf("%s -> %d", names(cols), n_each), collapse = ", "),
                names(cols)[cols == res$active], length(tab), 100 * max(tab) / sum(tab),
                if (small) sprintf("; %d with fewer than 10 cells", small) else ""),
        sprintf("各分辨率的簇数：%s。当前：分辨率 %s，%d 个簇（最大簇占 %.0f%% 细胞%s）。",
                paste(sprintf("%s -> %d", names(cols), n_each), collapse = "，"),
                names(cols)[cols == res$active], length(tab), 100 * max(tab) / sum(tab),
                if (small) sprintf("；%d 个簇少于 10 个细胞", small) else ""))
    })

    render_step_plot(output, input, "preview", function() {
      shiny::req(res$done, res$active)
      shiny::req(res$active %in% obj_meta_cols(rv$obj))
      cols <- res$cols
      cluster_plot(rv$obj, res$active, names(cols)[cols == res$active])
    }, name = "cluster")
  })
}
