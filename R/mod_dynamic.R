#' Module: Dynamic features (pseudotime-varying genes)
#'
#' Identify genes whose expression changes significantly along one or more
#' trajectory lineages (pseudotime), then visualise them as a DynamicHeatmap.
#' Requires a trajectory to have been computed first (lineages / pseudotime in
#' the object metadata).
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
      en = "Find genes whose expression rises or falls along a trajectory
            (pseudotime), and order them into a smooth dynamic heatmap.",
      zh = "找出沿轨迹（拟时序）表达上升或下降的基因，并将它们排列成平滑的动态热图。"),
    why  = list(
      en = "Cluster markers describe discrete states; dynamic features describe
            the continuous programme a cell runs as it differentiates, revealing
            waves of activation along the lineage.",
      zh = "簇标志基因描述离散状态；动态特征则描述细胞分化过程中运行的连续程序，
            揭示沿谱系的激活波次。"),
    how  = list(
      en = "<b>Run a trajectory step first</b> so the object carries lineages and
            pseudotime. Leave the lineage box empty to use all detected lineages,
            or list specific lineage names (comma-separated). Increase
            <b>candidate features</b> to scan more genes (slower).",
      zh = "<b>请先运行轨迹步骤</b>，使对象携带谱系与拟时序信息。谱系框留空则使用所有检测到的谱系，
            或填写特定谱系名称（逗号分隔）。增大<b>候选特征数</b>可扫描更多基因（更慢）。"),
    read = list(
      en = "Heatmap rows are genes ordered by where they peak along the
            trajectory; columns are cells in pseudotime. Blocks of genes
            switching on or off together are the dynamic programmes.",
      zh = "热图的行为按轨迹上峰值位置排序的基因，列为按拟时序排列的细胞。成块同时开启或关闭的基因即为动态程序。"),
    example = list(
      en = "Along a stem-to-mature lineage, stemness genes fade early while
               maturation markers switch on later — the heatmap shows this as a
               diagonal band.",
      zh = "沿干细胞到成熟细胞的谱系，干性基因较早减弱，而成熟标志基因较晚开启——
               热图会将其显示为一条对角带。")
  )
  controls <- shiny::tagList(
    label_with_help("Lineages (optional)",
                    "Comma-separated lineage names (e.g. Lineage1, Lineage2). Leave empty to use all detected lineages.",
                    label_zh = "谱系（可选）",
                    tip_zh = "逗号分隔的谱系名称（如 Lineage1, Lineage2）。留空则使用所有检测到的谱系。"),
    shiny::textInput(ns("lineages"), NULL, value = "",
                     placeholder = "Lineage1, Lineage2"),
    label_with_help("Candidate features",
                    "How many candidate genes to scan for dynamic behaviour. Higher = more thorough but slower.",
                    label_zh = "候选特征数",
                    tip_zh = "扫描多少候选基因以寻找动态行为。越高越彻底但越慢。"),
    shiny::numericInput(ns("n_candidates"), NULL, value = 200, min = 10, step = 10),
    run_button(ns("run"), "Detect dynamic features", "检测动态特征")
  )
  step_container(title = list(en = "Dynamic features", zh = "动态特征"),
                 subtitle = list(en = "Genes that switch on or off along a trajectory.",
                                 zh = "沿轨迹开启或关闭的基因。"),
                 explainer = explainer, controls = controls,
                 summary = shiny::uiOutput(ns("summary")),
                 preview = preview_plot_ui(ns("preview"),
                   guide = list(en = "Dynamic gene programmes will be drawn here.",
                                zh = "运行后，这里将绘制动态基因程序。"),
                   caption = list(en = "Gene expression dynamics along the trajectory.",
                                  zh = "基因表达沿轨迹的动态变化。")))
}

#' @rdname mod_dynamic
#' @keywords internal
mod_dynamic_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    res <- shiny::reactiveValues(done = FALSE, lineages = NULL, n = NA)

    # Parse the comma-separated lineage text into a character vector (or NULL).
    parse_lineages <- function(txt) {
      if (is.null(txt)) return(NULL)
      parts <- trimws(strsplit(txt, ",", fixed = TRUE)[[1]])
      parts <- parts[nzchar(parts)]
      if (length(parts) == 0) NULL else parts
    }

    shiny::observeEvent(input$run, {
      shiny::req(rv$obj)
      if (!require_pkgs("scop", "Dynamic features")) return(NULL)
      lineages <- parse_lineages(input$lineages)
      n_candidates <- as.integer(input$n_candidates)
      obj <- with_progress_notify({
        sc_dynamic(rv$obj, lineages = lineages, n_candidates = n_candidates)
      }, message = "Detecting dynamic features along pseudotime...")
      if (is.null(obj)) return(NULL)
      rv$obj <- obj
      res$done     <- TRUE
      res$lineages <- lineages
      res$n        <- n_candidates
      mark_done(rv, "dynamic")
      log_step(log_rv, "Dynamic features",
               params = list(lineages = lineages %||% "all",
                             n_candidates = n_candidates),
               code = sprintf(
                 "obj <- scop::RunDynamicFeatures(obj, lineages = %s, n_candidates = %d)",
                 if (is.null(lineages)) "NULL" else
                   paste0("c(", paste(sprintf('"%s"', lineages), collapse = ", "), ")"),
                 n_candidates))
      shiny::showNotification(
        i18n("Dynamic feature detection done. See the dynamic heatmap.",
             "动态特征检测完成。请查看动态热图。"),
        type = "message")
    })

    output$summary <- shiny::renderUI({
      if (!isTRUE(res$done)) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("Run a trajectory step first, then detect dynamic features.",
                               "请先运行轨迹步骤，然后检测动态特征。")))
      }
      lin_lab <- if (is.null(res$lineages)) i18n("all", "全部") else
        paste(res$lineages, collapse = ", ")
      bslib::layout_columns(
        col_widths = c(6, 6),
        stat_tile(i18n("Lineages", "谱系"), lin_lab),
        stat_tile(i18n("Candidates", "候选数"), format(res$n, big.mark = ","))
      )
    })

    output$preview <- render_scop_plot(function() {
      shiny::req(res$done)
      obj <- rv$obj
      shiny::req(obj)
      if (!require_pkgs("scop", "DynamicHeatmap")) return(NULL)
      # scop::DynamicHeatmap draws the pseudotime-ordered feature heatmap.
      tryCatch(
        sc_dynamicheatmap(obj, lineages = res$lineages),
        error = function(e) {
          shiny::showNotification(
            i18n(paste("Dynamic heatmap unavailable:", conditionMessage(e)),
                 paste("动态热图不可用：", conditionMessage(e))),
            type = "error", duration = 10)
          NULL
        })
    })
  })
}
