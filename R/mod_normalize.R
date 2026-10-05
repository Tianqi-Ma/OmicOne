#' Module 4: Normalization
#'
#' Put cells on a comparable scale so that differences reflect biology rather than
#' sequencing depth. Supports classic LogNormalize (default) or SCTransform,
#' always computed from the RNA counts.
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
            regularized negative binomial; its Pearson residuals replace
            scaling, and it can regress out the mitochondrial % or cell-cycle
            scores. Both always start from the RNA counts; switching back to
            LogNormalize removes the SCT assay.",
      zh = "<b>LogNormalize</b> 将每个细胞缩放到统一的总量后再取对数（稳健的默认方法）。<b>SCT</b> 用正则化负二项模型对计数建模，其 Pearson 残差替代缩放，并可回归掉线粒体比例或细胞周期评分。两者都从 RNA 计数出发；切回 LogNormalize 会删除 SCT assay。"),
    read = list(
      en = "Two panels, before and after: each point is a cell, x = library
            size (log10 UMIs), y = mean count of the 500 most-expressed genes
            (raw counts before, normalised counts after; log axis). Before
            normalisation the cloud climbs with depth (Spearman rho near 1);
            after, the trend should be much weaker.",
      zh = "前后两个面板：每个点是一个细胞，横轴 = 文库大小（log10 UMI），纵轴 = 表达最高的 500 个基因的平均计数（归一化前为原始计数，归一化后为归一化计数；对数轴）。归一化前点云随深度上升（Spearman rho 接近 1）；归一化后趋势应明显减弱。"),
    example = list(
      en = "A cell with 20,000 UMIs and one with 5,000 UMIs are put on the same
               scale. On the bundled PBMC 3k demo, LogNormalize takes rho from 0.99
               to 0.10.",
      zh = "一个有 20,000 个 UMI 的细胞与一个有 5,000 个 UMI 的细胞被放到同一尺度。在内置的 PBMC 3k 演示数据上，LogNormalize 将 rho 从 0.99 降到 0.10。")
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
    shiny::conditionalPanel(
      sprintf("input['%s'] == 'SCT'", ns("method")),
      label_with_help("Regress out (optional)",
                      "Metadata covariates removed from the SCT residuals. Offered only when present: percent.mt (after QC), S.Score / G2M.Score (after cell-cycle scoring).",
                      "回归掉的变量（可选）",
                      "从 SCT 残差中去除的元数据协变量。仅在存在时提供：percent.mt（质控后）、S.Score / G2M.Score（细胞周期评分后）。"),
      shiny::uiOutput(ns("regress_ui"))
    ),
    run_button(ns("run"), "Normalize", "归一化")
  )
  step_container(title = list(en = "Normalization", zh = "归一化"),
                 subtitle = list(en = "Make cells sequenced to different depths comparable.",
                                 zh = "使测序深度不同的细胞可相互比较。"),
                 explainer = explainer, controls = controls,
                 summary = shiny::uiOutput(ns("summary")),
                 preview = shiny::tagList(
                   shiny::uiOutput(ns("insight")),
                   preview_plot_ui(ns("preview"), download = TRUE,
                     guide = list(en = "Library size vs the mean count of the top 500 genes, before and after normalization, will be drawn here.",
                                  zh = "运行后，这里将绘制归一化前后文库大小与前 500 个高表达基因平均计数的关系。"),
                     caption = list(en = "One point = one cell; rho = Spearman correlation with library size (lower after = less depth effect).",
                                    zh = "每个点为一个细胞；rho = 与文库大小的 Spearman 相关（归一化后越低说明深度效应越小）。"))))
}

#' @rdname mod_normalize
#' @keywords internal
mod_normalize_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    res <- step_results(rv, "sc", done = FALSE, method = NULL, scale_factor = NULL,
                        vars = NULL, diag = NULL)

    output$regress_ui <- shiny::renderUI({
      cols <- intersect(c("percent.mt", "percent.ribo", "S.Score", "G2M.Score"),
                        obj_meta_cols(rv$obj))
      if (!length(cols)) {
        return(shiny::div(class = "omicone-note",
                          i18n("No covariates available yet (run QC or cell-cycle scoring first).",
                               "暂无可用协变量（请先运行质控或细胞周期评分）。")))
      }
      shiny::checkboxGroupInput(session$ns("vars"), NULL, choices = cols,
                                selected = intersect(shiny::isolate(input$vars), cols))
    })

    shiny::observeEvent(input$run, {
      shiny::req(rv$obj)
      if (!require_pkgs("Seurat", "Normalization")) return(NULL)
      method <- input$method
      sf <- num_input(input$scale_factor, 1)
      vars <- if (method == "SCT") intersect(input$vars, obj_meta_cols(rv$obj)) else NULL
      if (method == "LogNormalize" && is.na(sf)) {
        shiny::showNotification(i18n("Enter a positive scale factor.", "请输入正的缩放因子。"),
                                type = "error")
        return(NULL)
      }
      had_sct <- "SCT" %in% obj_assays(rv$obj)
      out <- with_progress_notify({
        o <- normalize_obj(rv$obj, method = method, scale_factor = sf,
                           vars_to_regress = vars)
        list(obj = o, diag = tryCatch(normalize_diag_data(o), error = function(e) NULL))
      }, message = "Normalizing counts...")
      if (is.null(out)) return(NULL)
      rv$obj <- out$obj
      res$done <- TRUE
      res$method <- method
      res$scale_factor <- sf
      res$vars <- vars
      res$diag <- out$diag
      mark_done(rv, "normalize")
      log_step(log_rv, "Normalization",
               params = list(method = method,
                             scale.factor = if (method == "LogNormalize") sf else NULL,
                             vars.to.regress = vars),
               code = normalize_log_code(method, sf, vars,
                                         dropped_sct = method == "LogNormalize" && had_sct))
      shiny::showNotification(i18n(sprintf("Normalization done (%s).", method),
                                   sprintf("归一化完成（%s）。", method)),
                              type = "message")
    })

    output$summary <- shiny::renderUI({
      if (!isTRUE(res$done)) return(shiny::div(class = "omicone-placeholder",
                                               i18n("Pick a method and click Normalize.",
                                                    "选择一种方法并点击归一化。")))
      shiny::tagList(
        stat_tile(i18n("Method", "方法"), res$method),
        if (identical(res$method, "LogNormalize")) {
          stat_tile(i18n("Scale factor", "缩放因子"), format(res$scale_factor, big.mark = ","))
        } else {
          stat_tile(i18n("Regressed", "回归变量"),
                    if (length(res$vars)) paste(res$vars, collapse = ", ") else i18n("none", "无"))
        },
        stat_tile(i18n("Active assay", "当前 assay"),
                  if (identical(res$method, "SCT")) "SCT" else "RNA")
      )
    })

    output$insight <- shiny::renderUI({
      dd <- res$diag
      if (is.null(dd)) return(NULL)
      rho <- attr(dd, "rho")
      insight_bar(
        sprintf("Correlation of the top genes' mean count with library size: Spearman rho %.2f before, %.2f after %s (%s cells shown). Residual correlation can also reflect real differences in cell size or type.",
                rho[["before"]], rho[["after"]], res$method,
                format(nrow(dd) / 2, big.mark = ",")),
        sprintf("头部基因平均计数与文库大小的相关性：%s 前 Spearman rho 为 %.2f，后为 %.2f（展示 %s 个细胞）。剩余相关也可能反映细胞大小或类型的真实差异。",
                res$method, rho[["before"]], rho[["after"]],
                format(nrow(dd) / 2, big.mark = ",")))
    })

    render_step_plot(output, input, "preview", function() {
      dd <- res$diag
      shiny::req(dd)
      normalize_diag_plot(dd, res$method)
    }, name = "normalize")
  })
}
