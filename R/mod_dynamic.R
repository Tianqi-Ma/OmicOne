#' Module: Dynamic features (pseudotime-varying genes)
#'
#' Identify genes whose expression changes along one or more trajectory
#' lineages (pseudotime columns written by the Trajectory step), then show
#' them as a DynamicHeatmap. Wraps scop::RunDynamicFeatures() via sc_dynamic().
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_dynamic
NULL

#' @rdname mod_dynamic
#' @keywords internal
mod_dynamic_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = list(en = "Dynamic features", zh = "动态特征"),
    what = list(
      en = "Find genes whose expression rises or falls along a lineage's
            pseudotime, and order them into a smoothed dynamic heatmap.",
      zh = "找出沿某条谱系拟时序表达上升或下降的基因，并将它们排列成平滑的动态热图。"),
    why  = list(
      en = "Cluster markers describe discrete states; dynamic features describe the
            programme a cell runs as it differentiates. The pseudotime was itself
            inferred from these expression data, so testing genes against it is
            double dipping: p-values come out too small (Neufeld et al. 2022).",
      zh = "簇标志基因描述离散状态；动态特征描述细胞分化过程中运行的程序。拟时序本身就是由这些表达数据推断出来的，再用它检验基因属于 double dipping：p 值会偏小（Neufeld et al. 2022）。"),
    how  = list(
      en = "<b>Run the Trajectory step first</b>; its pseudotime columns are listed
            here. Each gene is fitted with a GAM along pseudotime. Raise
            <b>candidate features</b> to scan more variable genes (slower).",
      zh = "<b>请先运行轨迹步骤</b>；其拟时序列会列在这里。每个基因沿拟时序用 GAM 拟合。增大<b>候选特征数</b>可扫描更多高变基因（更慢）。"),
    read = list(
      en = "Rows are genes ordered by where they peak, columns are cells in
            pseudotime order; colour is the smoothed, scaled expression. Diagonal
            bands are waves of activation. Use the ranking, not the p-values.",
      zh = "行为按峰值位置排序的基因，列为按拟时序排列的细胞；颜色为平滑并标准化后的表达。对角带即一波波的激活。请使用排序，而不是 p 值。"),
    example = list(
      en = "Along a stem-to-mature lineage, stemness genes fade early while
            maturation markers switch on late, drawing one diagonal band.",
      zh = "沿干细胞到成熟细胞的谱系，干性基因较早减弱，成熟标志基因较晚开启，形成一条对角带。")
  )
  controls <- shiny::tagList(
    label_with_help("Lineages (pseudotime columns)",
                    "Pseudotime columns produced by the Trajectory step.",
                    label_zh = "谱系（拟时序列）",
                    tip_zh = "由轨迹步骤生成的拟时序列。"),
    shiny::selectInput(ns("lineages"), NULL, choices = NULL, multiple = TRUE),
    label_with_help("Candidate features",
                    "Number of highly variable genes scanned for dynamic behaviour. Higher = more thorough but slower.",
                    label_zh = "候选特征数",
                    tip_zh = "扫描多少个高变基因以寻找动态行为。越高越全面但越慢。"),
    shiny::numericInput(ns("n_candidates"), NULL, value = 1000, min = 50, max = 5000, step = 50),
    run_button(ns("run"), "Detect dynamic features", "检测动态特征")
  )
  step_container(id = id, title = list(en = "Dynamic features", zh = "动态特征"),
                 subtitle = list(en = "Genes that switch on or off along a trajectory.",
                                 zh = "沿轨迹开启或关闭的基因。"),
                 explainer = explainer, controls = controls,
                 summary = shiny::uiOutput(ns("summary")),
                 preview = shiny::tagList(
                   shiny::uiOutput(ns("insight")),
                   preview_plot_ui(ns("preview"), download = TRUE,
                     guide = list(en = "Dynamic gene programmes will be drawn here.",
                                  zh = "运行后，这里将绘制动态基因程序。"),
                     caption = list(en = "Genes (rows) passing R² > 0.2, deviance explained > 0.2 and adjusted p < 0.05, ordered by peak; cells (columns) in pseudotime order.",
                                    zh = "通过 R² > 0.2、偏差解释度 > 0.2 且校正 p < 0.05 的基因（行），按峰值排序；细胞（列）按拟时序排列。"))))
}

#' @rdname mod_dynamic
#' @keywords internal
mod_dynamic_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    res <- step_results(rv, "sc", done = FALSE, lineages = NULL, n = NA, counts = NULL)

    available <- shiny::reactive({
      obj <- rv$obj
      rec <- obj_misc(obj, "omicone_trajectory")$pt_cols
      dynamic_lineage_cols(obj_meta(obj), recorded = rec)
    })

    shiny::observe({
      cols <- available()
      cur <- shiny::isolate(input$lineages)
      sel <- if (length(cur) && all(cur %in% cols)) cur else cols
      shiny::updateSelectInput(session, "lineages", choices = cols, selected = sel)
    })

    shiny::observeEvent(input$run, {
      shiny::req(rv$obj)
      lineages <- intersect(input$lineages %||% available(), available())
      n_candidates <- int_input(input$n_candidates, lo = 50, hi = 5000)
      shiny::req(!is.na(n_candidates))
      if (!length(lineages)) {
        shiny::showNotification(
          i18n("No pseudotime column found. Run the Trajectory step first.",
               "未找到拟时序列。请先运行轨迹步骤。"),
          type = "warning", duration = 10)
        return(NULL)
      }
      if (!require_pkgs("scop", "Dynamic features")) return(NULL)
      obj <- with_progress_notify({
        sc_dynamic(rv$obj, lineages = lineages, n_candidates = n_candidates)
      }, message = "Fitting genes along pseudotime...")
      if (is.null(obj)) return(NULL)
      rv$obj <- obj
      res$done     <- TRUE
      res$lineages <- lineages
      res$n        <- n_candidates
      res$counts   <- dynamic_counts(obj, lineages)
      mark_done(rv, "dynamic")
      log_step(log_rv, "Dynamic features",
               params = list(lineages = lineages, n_candidates = n_candidates),
               code = dynamic_log_code(lineages, n_candidates))
      shiny::showNotification(
        i18n("Dynamic feature detection done. See the dynamic heatmap.",
             "动态特征检测完成。请查看动态热图。"),
        type = "message")
    })

    output$summary <- shiny::renderUI({
      if (!isTRUE(res$done)) {
        msg <- if (length(available())) {
          i18n("Choose the lineages, then click <b>Detect dynamic features</b>.",
               "选择谱系，然后点击<b>检测动态特征</b>。")
        } else {
          i18n("No pseudotime yet: run the Trajectory step first.",
               "尚无拟时序：请先运行轨迹步骤。")
        }
        return(shiny::div(class = "omicone-placeholder", msg))
      }
      shiny::tagList(
        stat_tile(i18n("Lineages", "谱系"), length(res$lineages)),
        stat_tile(i18n("Candidates", "候选数"), format(res$n, big.mark = ",")),
        stat_tile(i18n("Dynamic genes", "动态基因"),
                  format(sum(res$counts), big.mark = ","))
      )
    })

    output$insight <- shiny::renderUI({
      if (!isTRUE(res$done)) return(NULL)
      per <- paste(sprintf("%s: %d", names(res$counts), res$counts), collapse = "; ")
      insight_bar(
        sprintf("Dynamic genes per lineage (of %s candidates): %s. The pseudotime was inferred from the same data, so these p-values are anti-conservative (double dipping); rank genes, do not report them as significant.",
                format(res$n, big.mark = ","), per),
        sprintf("各谱系的动态基因数（候选 %s 个）：%s。拟时序由同一数据推断而来，这些 p 值偏乐观（double dipping）；请用于排序，不要作为显著性结果报告。",
                format(res$n, big.mark = ","), per))
    })

    render_step_plot(output, input, "preview", function() {
      shiny::req(res$done, sum(res$counts) > 0)
      # Lineages without a passing gene would make DynamicHeatmap() fail.
      sc_dynamicheatmap(rv$obj, lineages = names(res$counts)[res$counts > 0])
    }, name = "dynamic")
  })
}
