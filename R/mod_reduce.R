#' Module 5: Feature selection & PCA
#'
#' Pick highly variable genes (HVGs), scale them, and run PCA to compress the data
#' into a handful of informative components used by every later step. After
#' SCTransform, SCT's own variable genes and Pearson residuals are used as is.
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
      en = "Choose how variable genes are ranked and how many to keep (<b>mvp</b>
            selects by mean/dispersion cut-offs, so the number is not used), then
            how many PCs to compute. After <b>SCT</b> normalisation the HVG
            controls do not apply: SCT's own variable genes and residuals are
            used. Re-running PCA removes Harmony / integrated spaces and 2-D maps
            built on the old PCA.",
      zh = "选择高变基因的排序方式和保留数量（<b>mvp</b> 按均值/离散度阈值筛选，不使用数量参数），再选择要计算的主成分数。若用 <b>SCT</b> 归一化，高变基因控件不适用：直接使用 SCT 自身的高变基因和残差。重新运行 PCA 会删除基于旧 PCA 的 Harmony / 整合空间和二维图。"),
    read = list(
      en = "Left: every gene's mean (x) against its variability statistic (y);
            red = selected variable genes, the top 10 labelled. Right: % of the
            total variance each PC explains; keep the PCs before the curve
            flattens.",
      zh = "左图：每个基因的均值（横轴）与其变异统计量（纵轴）；红色 = 入选的高变基因，并标注前 10 个。右图：每个主成分解释的总方差百分比；保留曲线变平之前的主成分。"),
    example = list(
      en = "From 20,000 genes you keep ~2,000 variable ones, then summarise them
               as 50 PCs; the first ~20 usually carry the real structure.",
      zh = "从 20,000 个基因中保留约 2,000 个高变基因，再概括为 50 个主成分；
               通常前约 20 个承载了真正的结构。")
  )
  controls <- shiny::tagList(
    shiny::uiOutput(ns("hvg_ui")),
    label_with_help("Number of principal components",
                    "How many PCs to compute. 50 is a common default; you rarely use them all downstream.",
                    "主成分数量",
                    "计算多少个主成分。50 是常用默认值；下游很少会全部用到。"),
    shiny::numericInput(ns("npcs"), NULL, value = 50, min = 2, max = 200, step = 1),
    run_button(ns("run"), "Select features & run PCA", "选择特征并运行 PCA")
  )
  step_container(id = id, title = list(en = "Feature selection & PCA", zh = "特征选择与 PCA"),
                 subtitle = list(en = "Highly variable genes, then PCA compression.",
                                 zh = "高变基因筛选，再做 PCA 压缩。"),
                 explainer = explainer, controls = controls,
                 summary = shiny::uiOutput(ns("summary")),
                 preview = shiny::tagList(
                   shiny::uiOutput(ns("insight")),
                   preview_plot_ui(ns("preview"), download = TRUE,
                     guide = list(en = "The variable-gene plot and the variance explained per PC will be drawn here.",
                                  zh = "运行后，这里将绘制高变基因图与每个主成分解释的方差。"),
                     caption = list(en = "Left: one point = one gene (red = variable). Right: % of total variance per PC.",
                                    zh = "左：每个点为一个基因（红色 = 高变基因）。右：每个主成分解释的总方差百分比。"))))
}

#' @rdname mod_reduce
#' @keywords internal
mod_reduce_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    res <- step_results(rv, "sc", done = FALSE, n_hvg = NA, npcs = NA, hvg = NULL,
                        pv = NULL, method = NULL)

    # The HVG controls only apply to a log-normalised RNA assay.
    output$hvg_ui <- shiny::renderUI({
      ns <- session$ns
      if (identical(obj_default_assay(rv$obj), "SCT")) {
        return(shiny::div(class = "omicone-note",
                          i18n("SCT is the active assay: its variable genes and Pearson residuals are used as they are, so the HVG method and number do not apply.",
                               "当前 assay 为 SCT：直接使用其高变基因和 Pearson 残差，高变基因方法与数量不适用。")))
      }
      meth <- c("vst" = "vst", "mvp" = "mvp", "dispersion" = "dispersion")
      shiny::tagList(
        label_with_help("HVG method",
                        "How variable genes are ranked. vst = variance-stabilizing (recommended); mvp = mean/dispersion cut-offs (gene number set by the cut-offs); dispersion = top genes by dispersion.",
                        "高变基因方法",
                        "高变基因的排序方式。vst = 方差稳定（推荐）；mvp = 均值/离散度阈值（基因数由阈值决定）；dispersion = 按离散度取前若干个基因。"),
        shiny::selectInput(ns("hvg_method"), NULL, choices = meth,
                           selected = keep_selected(shiny::isolate(input$hvg_method), meth, "vst")),
        shiny::conditionalPanel(
          sprintf("input['%s'] != 'mvp'", ns("hvg_method")),
          label_with_help("Number of variable genes",
                          "How many highly variable genes to keep. 2,000 is a common default.",
                          "高变基因数量",
                          "保留多少个高变基因。2,000 是常用默认值。"),
          shiny::sliderInput(ns("n_hvg"), NULL, min = 500, max = 5000,
                             value = shiny::isolate(input$n_hvg) %||% 2000, step = 100)
        ),
        shiny::conditionalPanel(
          sprintf("input['%s'] == 'mvp'", ns("hvg_method")),
          shiny::div(class = "omicone-note",
                     i18n("mvp keeps genes with mean 0.1-8 and scaled dispersion above 1 (Seurat defaults); the number of genes follows from these cut-offs.",
                          "mvp 保留均值在 0.1-8 且标准化离散度大于 1 的基因（Seurat 默认值）；基因数量由这些阈值决定。")))
      )
    })

    shiny::observeEvent(input$run, {
      shiny::req(rv$obj)
      if (!require_pkgs("Seurat", "Feature selection & PCA")) return(NULL)
      sct <- identical(obj_default_assay(rv$obj), "SCT")
      n_hvg <- int_input(input$n_hvg %||% 2000, 50)
      npcs <- int_input(input$npcs, 2, 500)
      hvg_method <- if (sct) "SCT" else input$hvg_method %||% "vst"
      if (is.na(npcs) || is.na(n_hvg)) {
        shiny::showNotification(i18n("Enter whole numbers for the gene and PC counts.",
                                     "请为基因数和主成分数输入整数。"), type = "error")
        return(NULL)
      }
      dropped <- obj_reductions(rv$obj)
      out <- with_progress_notify({
        o <- reduce_obj(rv$obj, n_hvg = n_hvg, npcs = npcs,
                        hvg_method = if (sct) "vst" else hvg_method)
        list(obj = o, hvg = tryCatch(hvg_plot_data(o, hvg_method), error = function(e) NULL),
             pv = pca_variance(o))
      }, message = "Selecting features and running PCA...")
      if (is.null(out)) return(NULL)
      obj <- out$obj
      rv$obj <- obj
      res$done <- TRUE
      res$method <- hvg_method
      res$n_hvg <- length(SeuratObject::VariableFeatures(obj))
      res$hvg <- out$hvg
      res$pv <- out$pv
      res$npcs <- nrow(out$pv)
      mark_done(rv, "reduce")
      log_step(log_rv, "Feature selection & PCA",
               params = list(hvg_method = hvg_method,
                             n_hvg = if (hvg_method %in% c("vst", "dispersion")) n_hvg else NULL,
                             npcs = res$npcs, seed = 42),
               code = reduce_log_code(sct, hvg_method, n_hvg, res$npcs, 42, dropped))
      shiny::showNotification(i18n(sprintf("PCA done: %s HVGs, %d PCs.",
                                           format(res$n_hvg, big.mark = ","), res$npcs),
                                   sprintf("PCA 完成：%s 个高变基因，%d 个主成分。",
                                           format(res$n_hvg, big.mark = ","), res$npcs)),
                              type = "message")
    })

    output$summary <- shiny::renderUI({
      if (!isTRUE(res$done)) return(shiny::div(class = "omicone-placeholder",
                                               i18n("Set parameters and click Select features & run PCA.",
                                                    "设置参数并点击选择特征并运行 PCA。")))
      shiny::tagList(
        stat_tile(i18n("HVG method", "高变基因方法"), res$method),
        stat_tile(i18n("Variable genes", "高变基因"), format(res$n_hvg, big.mark = ",")),
        stat_tile(i18n("Principal components", "主成分"), format(res$npcs, big.mark = ","))
      )
    })

    output$insight <- shiny::renderUI({
      pv <- res$pv
      if (is.null(pv)) return(NULL)
      k <- min(10, nrow(pv))
      total <- identical(pv$basis[1], "total")
      insight_bar(
        sprintf("%d PCs from %s variable genes; the first %d explain %.1f%% of the %s variance. Choose the dims for clustering where the curve flattens, not from this number alone.",
                nrow(pv), format(res$n_hvg, big.mark = ","), k, pv$cum[k],
                if (total) "total scaled" else "computed PCs'"),
        sprintf("基于 %s 个高变基因得到 %d 个主成分；前 %d 个解释了%s %.1f%% 的方差。聚类维度应在曲线变平处选择，而不只看这个数字。",
                format(res$n_hvg, big.mark = ","), nrow(pv), k,
                if (total) "缩放后总" else "已计算主成分", pv$cum[k]))
    })

    render_step_plot(output, input, "preview", function() {
      shiny::req(res$pv)
      reduce_plot(res$hvg, res$pv)
    }, name = "reduce")
  })
}
