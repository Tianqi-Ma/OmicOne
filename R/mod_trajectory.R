#' Module: Trajectory / pseudotime
#'
#' Order cells along inferred developmental / differentiation trajectories and
#' assign each cell a pseudotime. Wraps scop's trajectory engines
#' (Slingshot / Monocle2 / Monocle3 / PAGA / Palantir / WOT) via sc_trajectory().
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_trajectory
NULL

#' @rdname mod_trajectory
#' @keywords internal
mod_trajectory_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = list(en = "Trajectory / pseudotime", zh = "轨迹 / 拟时序"),
    what = list(
      en = "Infer how cells transition between states and assign each cell a
            <i>pseudotime</i> along that path.",
      zh = "推断细胞在不同状态间如何转变，并沿该路径为每个细胞赋予一个<i>拟时序</i>值。"),
    why  = list(
      en = "Development and differentiation are continuous, but a snapshot mixes
            all stages together. Ordering cells by pseudotime reconstructs the
            process and reveals genes that change along it.",
      zh = "发育与分化是连续的，但单次快照会把所有阶段混在一起。按拟时序排序细胞可以重建这一过程，并揭示沿途变化的基因。"),
    how  = list(
      en = "<b>Slingshot</b> is a robust R-only default. <b>Monocle2/3</b> offer
            alternative graph models. Methods marked <b>*</b>
            (PAGA / Palantir / WOT) run in Python and need a scop conda
            environment (scop::PrepareEnv). Optionally name a start cluster to
            root the trajectory.",
      zh = "<b>Slingshot</b> 是稳健的纯 R 默认方法。<b>Monocle2/3</b> 提供另一类图模型。标有 <b>*</b> 的方法（PAGA / Palantir / WOT）在 Python 中运行，需要 scop 的 conda 环境（scop::PrepareEnv）。可选地指定起始簇，以确定轨迹的根。"),
    read = list(
      en = "The tree or curve is the inferred trajectory; colour is pseudotime
            measured from the chosen root. Branches are fate decisions — check
            that the root sits in the stem-like cluster, otherwise the arrow of
            time runs backwards.",
      zh = "树或曲线是推断出的轨迹；颜色为从所选起点计算的拟时序。分叉即命运决定——请确认起点位于干性的簇，否则时间方向是反的。"),
    example = list(
      en = "Starting from hematopoietic stem cells, pseudotime should increase
               smoothly toward the mature myeloid and lymphoid tips.",
      zh = "从造血干细胞出发，拟时序应平滑地增大，直至成熟的髓系和淋巴系终端。")
  )
  controls <- shiny::tagList(
    label_with_help("Method",
                    "Slingshot/Monocle are R; PAGA/Palantir/WOT (*) run in Python via a scop conda env.",
                    label_zh = "方法",
                    tip_zh = "Slingshot/Monocle 为 R 实现；PAGA/Palantir/WOT（*）通过 scop 的 conda 环境在 Python 中运行。"),
    shiny::selectInput(ns("method"), NULL,
                       choices = c("Slingshot" = "slingshot",
                                   "Monocle2"  = "monocle2",
                                   "Monocle3"  = "monocle3",
                                   "PAGA *"    = "paga",
                                   "Palantir *" = "palantir",
                                   "WOT *"     = "wot"),
                       selected = "slingshot"),
    shiny::div(class = "omicone-note",
               i18n("* PAGA / Palantir / WOT need a Python conda environment (run scop::PrepareEnv() once).",
                    "* PAGA / Palantir / WOT 需要 Python conda 环境（首次使用请运行 scop::PrepareEnv()）。")),
    label_with_help("Group by (metadata column)",
                    "Cell grouping used to build the trajectory graph (e.g. seurat_clusters, celltype).",
                    label_zh = "分组依据（元数据列）",
                    tip_zh = "用于构建轨迹图的细胞分组（例如 seurat_clusters、celltype）。"),
    shiny::selectInput(ns("group_by"), NULL, choices = NULL),
    label_with_help("Start cluster (optional)",
                    "Name of the group to root the trajectory at (leave blank to auto-detect).",
                    label_zh = "起始簇（可选）",
                    tip_zh = "作为轨迹根节点的分组名称（留空则自动判定）。"),
    shiny::textInput(ns("start"), NULL, placeholder = "e.g. HSC / 0"),
    run_button(ns("run"), "Run trajectory", "运行轨迹分析")
  )
  step_container(
    title     = list(en = "Trajectory / pseudotime", zh = "轨迹 / 拟时序"),
    subtitle  = list(en = "Order cells along a differentiation path.",
                     zh = "沿分化路径为细胞排序。"),
    explainer = explainer,
    controls  = controls,
    summary   = shiny::uiOutput(ns("summary")),
    preview   = preview_plot_ui(ns("preview"),
      guide = list(en = "The inferred trajectory and pseudotime will be drawn here.",
                   zh = "运行后，这里将绘制推断的轨迹与拟时序。"),
      caption = list(en = "Trajectory coloured by pseudotime from the chosen root.",
                     zh = "按拟时序着色的轨迹（起点为所选根节点）。"))
  )
}

#' @rdname mod_trajectory
#' @keywords internal
mod_trajectory_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    res <- shiny::reactiveValues(done = FALSE, method = NULL, group_by = NULL,
                                 pt_col = NULL)

    # Keep the group-by selector in sync with the current object.
    shiny::observe({
      obj <- rv$obj
      cols <- if (is.null(obj)) character(0) else obj_meta_cols(obj)
      sel <- if ("celltype" %in% cols) "celltype"
             else if ("seurat_clusters" %in% cols) "seurat_clusters"
             else if (length(cols)) cols[1] else NULL
      shiny::updateSelectInput(session, "group_by", choices = cols, selected = sel)
    })

    # Best-guess a pseudotime metadata column that a trajectory run just added.
    guess_pt_col <- function(before, after) {
      new_cols <- setdiff(after, before)
      pat <- "pseudotime|Pseudotime|Lineage|dpt|palantir|Palantir|latent_time"
      hit <- new_cols[grepl(pat, new_cols)]
      if (length(hit)) return(hit[1])
      hit <- after[grepl(pat, after)]
      if (length(hit)) hit[1] else NULL
    }

    shiny::observeEvent(input$run, {
      shiny::req(rv$obj)
      if (!require_pkgs("scop", "Trajectory / pseudotime")) return(NULL)
      # Python-backed methods additionally require a conda environment.
      if (input$method %in% c("paga", "palantir", "wot")) {
        shiny::showNotification(
          i18n("This method runs in Python. If it fails, run scop::PrepareEnv() to set up the conda environment.",
               "该方法在 Python 中运行。如果失败，请运行 scop::PrepareEnv() 配置 conda 环境。"),
          type = "warning", duration = 8)
      }
      group_by <- input$group_by
      shiny::req(group_by)
      start <- trimws(input$start %||% "")
      before <- obj_meta_cols(rv$obj)
      obj <- with_progress_notify({
        args <- list(rv$obj, method = input$method, group_by = group_by)
        if (nzchar(start)) args$start <- start
        do.call(sc_trajectory, args)
      }, message = sprintf("Running %s trajectory...", input$method))
      if (is.null(obj)) return(NULL)
      rv$obj <- obj
      res$done     <- TRUE
      res$method   <- input$method
      res$group_by <- group_by
      res$pt_col   <- guess_pt_col(before, obj_meta_cols(obj))
      mark_done(rv, "trajectory")
      log_step(log_rv, "Trajectory",
               params = list(method = input$method, group_by = group_by,
                             start = if (nzchar(start)) start else NULL),
               code = sprintf(
                 'obj <- sc_trajectory(obj, method="%s", group_by="%s"%s)',
                 input$method, group_by,
                 if (nzchar(start)) sprintf(', start="%s"', start) else ""))
      shiny::showNotification(
        i18n(sprintf("Trajectory (%s) finished on '%s'.", input$method, group_by),
             sprintf("轨迹分析（%s）已完成，分组：%s。", input$method, group_by)),
        type = "message")
    })

    output$summary <- shiny::renderUI({
      if (!isTRUE(res$done)) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("Pick a method and group, then click <b>Run trajectory</b>.",
                               "选择方法与分组，然后点击<b>运行轨迹分析</b>。")))
      }
      bslib::layout_columns(
        col_widths = c(4, 4, 4),
        stat_tile(i18n("Method", "方法"), res$method),
        stat_tile(i18n("Group by", "分组依据"), res$group_by),
        stat_tile(i18n("Pseudotime", "拟时序"),
                  res$pt_col %||% i18n("n/a", "无"))
      )
    })

    output$preview <- render_scop_plot(function() {
      shiny::req(res$done)
      # Prefer colouring the embedding by the pseudotime column; fall back to a
      # dim plot coloured by the grouping if no pseudotime column was found.
      pt <- res$pt_col
      p <- NULL
      if (!is.null(pt)) {
        p <- tryCatch(sc_featureplot(rv$obj, features = pt),
                      error = function(e) NULL)
      }
      if (is.null(p)) {
        p <- tryCatch(sc_dimplot(rv$obj, group_by = res$group_by),
                      error = function(e) NULL)
      }
      p
    })
  })
}
