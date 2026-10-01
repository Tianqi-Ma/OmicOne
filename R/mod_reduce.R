#' Module 5: Feature selection & PCA
#'
#' Pick highly variable genes (HVGs), scale them, and run PCA to compress the data
#' into a handful of informative components used by every later step.
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_reduce
NULL

#' @rdname mod_reduce
#' @keywords internal
mod_reduce_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = list(en = "Feature selection & PCA", zh = "特征选择与 PCA"),
    what = list(
      en = "Select the most informative genes (highly variable features), then run
            <b>PCA</b> to summarise them as a few principal components.",
      zh = "选出信息量最大的基因（高变基因），再运行 <b>PCA</b> 将它们概括为少数几个主成分。"),
    why  = list(
      en = "Most genes vary little between cells and just add noise. Focusing on
            highly variable genes and compressing them with PCA makes clustering
            and embedding faster and cleaner.",
      zh = "大多数基因在细胞间变化很小，只会增加噪声。聚焦高变基因并用 PCA 压缩它们，
            能让聚类和嵌入更快、更干净。"),
    how  = list(
      en = "Choose how variable genes are ranked, how many to keep, and how many
            principal components to compute. The elbow plot on the right shows how
            much variation each component captures &mdash; keep components before the
            curve flattens.",
      zh = "选择高变基因的排序方式、保留多少个，以及要计算多少个主成分。右侧的肘部图显示
            每个主成分捕获了多少变异 &mdash; 保留曲线变平之前的主成分。"),
    read = list(
      en = "Left: the most variable genes — these drive everything downstream.
            Right: the elbow plot; keep the PCs before the curve flattens —
            anything past the elbow is mostly noise.",
      zh = "左图：变异最大的基因——它们驱动下游所有分析。右图：肘部图；保留曲线变平之前的主成分——拐点之后基本是噪声。"),
    example = list(
      en = "From 20,000 genes you keep ~2,000 variable ones, then summarise them
               as 50 PCs; the first ~20 usually carry the real structure.",
      zh = "从 20,000 个基因中保留约 2,000 个高变基因，再概括为 50 个主成分；
               通常前约 20 个承载了真正的结构。")
  )
  controls <- shiny::tagList(
    label_with_help("HVG method",
                    "How variable genes are ranked. vst = variance-stabilizing (recommended); mvp/dispersion = mean-variance/dispersion based.",
                    "高变基因方法",
                    "高变基因的排序方式。vst = 方差稳定（推荐）；mvp/dispersion = 基于均值-方差/离散度。"),
    shiny::selectInput(ns("hvg_method"), NULL,
                       choices = c("vst" = "vst", "mvp" = "mvp",
                                   "dispersion" = "dispersion"),
                       selected = "vst"),
    label_with_help("Number of variable genes",
                    "How many highly variable genes to keep. 2,000 is a common default.",
                    "高变基因数量",
                    "保留多少个高变基因。2,000 是常用默认值。"),
    shiny::sliderInput(ns("n_hvg"), NULL, min = 500, max = 5000, value = 2000,
                       step = 100),
    label_with_help("Number of principal components",
                    "How many PCs to compute. 50 is a common default; you rarely use them all downstream.",
                    "主成分数量",
                    "计算多少个主成分。50 是常用默认值；下游很少会全部用到。"),
    shiny::numericInput(ns("npcs"), NULL, value = 50, min = 2, max = 200, step = 1),
    run_button(ns("run"), "Select features & run PCA", "选择特征并运行 PCA")
  )
  step_container(title = list(en = "Feature selection & PCA", zh = "特征选择与 PCA"),
                 subtitle = list(en = "Highly variable genes, then PCA compression.",
                                 zh = "高变基因筛选，再做 PCA 压缩。"),
                 explainer = explainer, controls = controls,
                 summary = shiny::uiOutput(ns("summary")),
                 preview = preview_plot_ui(ns("preview"),
                   guide = list(en = "Variable features and the PCA elbow plot will be drawn here.",
                                zh = "运行后，这里将绘制高变基因图与 PCA 肘部图。"),
                   caption = list(en = "Highly variable genes (left) and variance captured per PC (right).",
                                  zh = "高变基因（左）与每个主成分捕获的方差（右）。")))
}

#' @rdname mod_reduce
#' @keywords internal
mod_reduce_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    res <- shiny::reactiveValues(done = FALSE, n_hvg = NA, npcs = NA, stdev = NULL)

    shiny::observeEvent(input$run, {
      shiny::req(rv$obj)
      if (!require_pkgs("Seurat", "Feature selection & PCA")) return(NULL)
      n_hvg <- input$n_hvg
      npcs <- input$npcs
      hvg_method <- input$hvg_method
      obj <- with_progress_notify({
        reduce_obj(rv$obj, n_hvg = n_hvg, npcs = npcs, hvg_method = hvg_method)
      }, message = "Selecting features and running PCA...")
      if (is.null(obj)) return(NULL)
      rv$obj <- obj
      res$done <- TRUE
      res$n_hvg <- length(tryCatch(Seurat::VariableFeatures(obj),
                                   error = function(e) character(0)))
      res$stdev <- tryCatch(obj@reductions$pca@stdev, error = function(e) NULL)
      res$npcs <- length(res$stdev)
      mark_done(rv, "reduce")
      log_step(log_rv, "Feature selection & PCA",
               params = list(hvg_method = hvg_method, n_hvg = n_hvg, npcs = npcs),
               code = sprintf(paste0(
                 "obj <- Seurat::FindVariableFeatures(obj, selection.method = '%s', nfeatures = %d)\n",
                 "obj <- Seurat::ScaleData(obj)\n",
                 "obj <- Seurat::RunPCA(obj, npcs = %d)"),
                 hvg_method, n_hvg, npcs))
      shiny::showNotification(sprintf("PCA done: %d HVGs, %d PCs.",
                                      res$n_hvg, res$npcs), type = "message")
    })

    output$summary <- shiny::renderUI({
      if (!isTRUE(res$done)) return(shiny::div(class = "omicone-placeholder",
                                               i18n("Set parameters and click Select features & run PCA.",
                                                    "设置参数并点击选择特征并运行 PCA。")))
      bslib::layout_columns(
        col_widths = c(6, 6),
        stat_tile(i18n("Variable genes", "高变基因"), format(res$n_hvg, big.mark = ",")),
        stat_tile(i18n("Principal components", "主成分"), format(res$npcs, big.mark = ","))
      )
    })

    output$preview <- render_preview_plot(function() {
      stdev <- res$stdev
      shiny::req(stdev)
      df <- data.frame(PC = seq_along(stdev), stdev = stdev)
      df$text <- sprintf("PC%d\nstdev=%.3f", df$PC, df$stdev)
      ggplot2::ggplot(df, ggplot2::aes(x = PC, y = stdev)) +
        ggplot2::geom_line(colour = "#7d8b8f", linewidth = 0.4) +
        ggplot2::geom_point(colour = sc_palette(1), size = 1.4) +
        ggplot2::labs(x = "Principal component", y = "Standard deviation",
                      title = "PCA elbow plot / PCA 肘部图") +
        omicone_theme()
    })
  })
}
