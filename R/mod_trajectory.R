#' Module: Trajectory / pseudotime
#'
#' Order cells along an inferred differentiation path and give each cell a
#' pseudotime, with one of scop's engines (Slingshot / Monocle2 / Monocle3 /
#' PAGA / Palantir / WOT) via sc_trajectory(). Each method's own rooting
#' argument is filled from the start group the user picks.
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
      en = "Infer how cells move between states and give each cell a
            <i>pseudotime</i>: its position along that path from a chosen root.",
      zh = "推断细胞在状态之间如何转变，并为每个细胞赋予<i>拟时序</i>：即它从所选根出发沿路径所处的位置。"),
    why  = list(
      en = "Differentiation is continuous but a snapshot mixes all stages; ordering
            the cells reconstructs the process. Fitting on a 2-D UMAP would follow
            the drawing's distortions, so Slingshot is fitted in the PCA (or
            integrated) space.",
      zh = "分化是连续的，但单次快照会把各阶段混在一起；为细胞排序可以重建这一过程。在二维 UMAP 上拟合会跟随作图的形变，因此 Slingshot 在 PCA（或整合后）空间中拟合。"),
    how  = list(
      en = "<b>Slingshot</b> (R) is the default and uses the first 20 dimensions of
            the space shown. Pick the <b>start group</b> from the grouping's levels:
            Monocle3, PAGA and Palantir require it; Slingshot without one picks a
            root itself. PAGA / Palantir / WOT (<b>*</b>) run in Python
            (scop::PrepareEnv()); WOT needs a numeric sampling-time column.",
      zh = "<b>Slingshot</b>（R）为默认方法，使用所示空间的前 20 维。请从分组的取值中选择<b>起始分组</b>：Monocle3、PAGA 和 Palantir 必须指定；Slingshot 未指定时会自行选根。PAGA / Palantir / WOT（<b>*</b>）在 Python 中运行（scop::PrepareEnv()）；WOT 需要数值型的采样时间列。"),
    read = list(
      en = "Slingshot: cells coloured by group with one curve per lineage, arrows
            toward later pseudotime. Other methods: cells coloured by pseudotime
            (Monocle adds its principal graph). Check the root sits in the
            stem-like group, otherwise time runs backwards.",
      zh = "Slingshot：细胞按分组着色，每条谱系一条曲线，箭头指向更晚的拟时序。其他方法：细胞按拟时序着色（Monocle 会叠加其主图）。请确认根位于干性分组，否则时间方向是反的。"),
    example = list(
      en = "Rooted at haematopoietic stem cells, pseudotime should rise smoothly
            toward the myeloid and lymphoid tips, giving two lineages.",
      zh = "以造血干细胞为根时，拟时序应平滑地增大至髓系和淋巴系终端，得到两条谱系。")
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
                    "Cell grouping the trajectory is built on (e.g. seurat_clusters, celltype).",
                    label_zh = "分组依据（元数据列）",
                    tip_zh = "构建轨迹所依据的细胞分组（例如 seurat_clusters、celltype）。"),
    shiny::selectInput(ns("group_by"), NULL, choices = NULL),
    label_with_help("Start group (root)",
                    "Group the pseudotime starts from. Required for Monocle3, PAGA and Palantir.",
                    label_zh = "起始分组（根）",
                    tip_zh = "拟时序的起点分组。Monocle3、PAGA 和 Palantir 必须指定。"),
    shiny::selectInput(ns("start"), NULL, choices = c("(none)" = "")),
    shiny::conditionalPanel(
      sprintf("['slingshot','paga','palantir'].indexOf(input['%s']) >= 0", ns("method")),
      label_with_help("Fitting space",
                      "Linear space the trajectory is fitted in (PCA or the integrated space).",
                      label_zh = "拟合空间",
                      tip_zh = "拟合轨迹所用的线性空间（PCA 或整合后的空间）。"),
      shiny::selectInput(ns("reduction"), NULL, choices = NULL)),
    shiny::conditionalPanel(
      sprintf("input['%s'] == 'slingshot'", ns("method")),
      label_with_help("Dimensions (Slingshot)",
                      "Number of leading dimensions of the fitting space Slingshot uses.",
                      label_zh = "维数（Slingshot）",
                      tip_zh = "Slingshot 使用拟合空间的前多少维。"),
      shiny::numericInput(ns("dims"), NULL, value = 20, min = 2, max = 100, step = 1)),
    shiny::conditionalPanel(
      sprintf("input['%s'] == 'wot'", ns("method")),
      label_with_help("Sampling time column (WOT)",
                      "Numeric metadata column with the real collection time of each cell (e.g. day).",
                      label_zh = "采样时间列（WOT）",
                      tip_zh = "数值型元数据列，记录每个细胞的真实采样时间（例如天数）。"),
      shiny::selectInput(ns("time_col"), NULL, choices = NULL)),
    run_button(ns("run"), "Run trajectory", "运行轨迹分析")
  )
  step_container(id = id, 
    title     = list(en = "Trajectory / pseudotime", zh = "轨迹 / 拟时序"),
    subtitle  = list(en = "Order cells along a differentiation path.",
                     zh = "沿分化路径为细胞排序。"),
    explainer = explainer,
    controls  = controls,
    summary   = shiny::uiOutput(ns("summary")),
    preview   = shiny::tagList(
      shiny::uiOutput(ns("insight")),
      preview_plot_ui(ns("preview"), download = TRUE,
        guide = list(en = "The inferred trajectory and pseudotime will be drawn here.",
                     zh = "运行后，这里将绘制推断的轨迹与拟时序。"),
        caption = list(
          en = "Slingshot: cells coloured by group, one curve per lineage. Other methods: cells coloured by pseudotime (Monocle adds its graph). The line above names what this run drew.",
          zh = "Slingshot：细胞按分组着色，每条谱系一条曲线。其他方法：细胞按拟时序着色（Monocle 叠加其主图）。图上方的说明给出本次绘制的内容。")))
  )
}

#' @rdname mod_trajectory
#' @keywords internal
mod_trajectory_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    res <- step_results(rv, "sc", done = FALSE, method = NULL, group_by = NULL,
                        start = NULL, reduction = NULL, nd = NA, pt_cols = character(0),
                        fate_cols = character(0), n_na = NA, n_cells = NA)

    shiny::observe({
      obj <- rv$obj
      cols <- categorical_cols(obj_meta(obj))
      def <- default_group_col(cols, obj_misc(obj, "omicone_cluster_col"))
      shiny::updateSelectInput(session, "group_by", choices = cols,
                               selected = keep_selected(shiny::isolate(input$group_by),
                                                        cols, def))
      reds <- traj_reduction_choices(obj)
      shiny::updateSelectInput(session, "reduction", choices = reds,
                               selected = keep_selected(shiny::isolate(input$reduction),
                                                        reds, traj_default_reduction(obj)))
      tcols <- traj_time_cols(obj_meta(obj))
      shiny::updateSelectInput(session, "time_col", choices = tcols,
                               selected = keep_selected(shiny::isolate(input$time_col), tcols))
    })

    # Start-group choices follow the grouping.
    shiny::observe({
      lv <- group_levels(obj_meta(rv$obj), input$group_by)
      choices <- c("(none)" = "", stats::setNames(lv, lv))
      shiny::updateSelectInput(session, "start", choices = choices,
                               selected = keep_selected(shiny::isolate(input$start),
                                                        unname(choices), ""))
    })

    shiny::observeEvent(input$run, {
      shiny::req(rv$obj)
      method <- input$method
      group_by <- input$group_by
      start <- input$start %||% ""
      reduction <- input$reduction %||% traj_default_reduction(rv$obj)
      dims <- int_input(input$dims, lo = 2, hi = 100)
      time_col <- input$time_col
      shiny::req(method, group_by, !is.na(dims))
      if (!require_pkgs("scop", "Trajectory / pseudotime")) return(NULL)
      md <- obj_meta(rv$obj)
      bad <- traj_check_start(method, start, group_levels(md, group_by),
                              reductions = obj_reductions(rv$obj), time_col = time_col,
                              time_values = if (length(time_col) && time_col %in% names(md))
                                md[[time_col]])
      if (!is.null(bad)) {
        shiny::showNotification(i18n(bad$en, bad$zh), type = "warning", duration = 10)
        return(NULL)
      }
      if (method %in% c("paga", "palantir", "wot")) {
        shiny::showNotification(
          i18n("This method runs in Python. If it fails, run scop::PrepareEnv() to set up the conda environment.",
               "该方法在 Python 中运行。如果失败，请运行 scop::PrepareEnv() 配置 conda 环境。"),
          type = "warning", duration = 8)
      }
      out <- with_progress_notify({
        sc_trajectory(rv$obj, method = method, group_by = group_by, start = start,
                      reduction = reduction, dims = dims, time_col = time_col)
      }, message = sprintf("Running %s trajectory...", method))
      if (is.null(out)) return(NULL)
      res$done      <- TRUE
      res$method    <- method
      res$group_by  <- group_by
      res$start     <- if (nzchar(start)) start else NULL
      res$reduction <- reduction
      res$nd        <- out$nd
      res$pt_cols   <- out$pt_cols
      res$fate_cols <- out$fate_cols
      res$n_cells   <- ncol(out$obj)
      res$n_na      <- if (length(out$pt_cols)) {
        sum(rowSums(!is.na(obj_meta(out$obj)[, out$pt_cols, drop = FALSE])) == 0)
      } else {
        NA
      }
      if (!length(out$pt_cols) && !length(out$fate_cols)) {
        shiny::showNotification(
          i18n(sprintf("%s produced no pseudotime column; nothing was stored.", method),
               sprintf("%s 未产生拟时序列；未保存任何结果。", method)),
          type = "warning", duration = 10)
        return(NULL)
      }
      rv$obj <- out$obj
      mark_done(rv, "trajectory")
      log_step(log_rv, "Trajectory",
               params = list(method = method, group_by = group_by,
                             start = if (nzchar(start)) start else NULL,
                             reduction = if (method %in% c("slingshot", "paga", "palantir"))
                               reduction else NULL,
                             dims = if (method == "slingshot") out$nd else NULL,
                             time_col = if (method == "wot") time_col else NULL),
               code = trajectory_log_code(method, group_by, start = start,
                                          reduction = reduction, nd = out$nd,
                                          time_col = time_col))
      shiny::showNotification(
        i18n(sprintf("Trajectory (%s) finished on '%s'.", method, group_by),
             sprintf("轨迹分析（%s）已完成，分组：%s。", method, group_by)),
        type = "message")
    })

    output$summary <- shiny::renderUI({
      if (!isTRUE(res$done)) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("Pick a method, a grouping and a start group, then click <b>Run trajectory</b>.",
                               "选择方法、分组与起始分组，然后点击<b>运行轨迹分析</b>。")))
      }
      n_pt <- length(res$pt_cols)
      shiny::tagList(
        stat_tile(i18n("Method", "方法"), res$method),
        stat_tile(i18n("Group by", "分组依据"), res$group_by),
        stat_tile(i18n("Start", "起点"), res$start %||% i18n("auto", "自动")),
        stat_tile(i18n("Pseudotime columns", "拟时序列"),
                  if (n_pt) n_pt else i18n("none", "无"))
      )
    })

    output$insight <- shiny::renderUI({
      if (!isTRUE(res$done)) return(NULL)
      if (!length(res$pt_cols) && !length(res$fate_cols)) {
        return(insight_bar(
          sprintf("<b>%s</b> produced no pseudotime in this run, so nothing is shown. Check the R console for scop's messages.",
                  res$method),
          sprintf("<b>%s</b> 本次运行未产生拟时序，因此不显示任何结果。请查看 R 控制台中 scop 的提示信息。",
                  res$method)))
      }
      root_en <- if (is.null(res$start)) {
        "none chosen, so the method picked one and the direction may be reversed"
      } else {
        sprintf("<b>%s</b>", res$start)
      }
      root_zh <- if (is.null(res$start)) {
        "未指定，由方法自行选根，方向可能相反"
      } else {
        sprintf("<b>%s</b>", res$start)
      }
      drawn <- switch(
        res$method,
        slingshot = c(sprintf("%d Slingshot lineage(s), fitted on the first %d dims of %s; curves are smoothed onto the map for display.",
                              length(res$pt_cols), res$nd, res$reduction),
                      sprintf("%d 条 Slingshot 谱系，在 %s 的前 %d 维上拟合；曲线为显示而平滑映射到二维图上。",
                              length(res$pt_cols), res$reduction, res$nd)),
        monocle2 = c("Monocle2 pseudotime on its DDRTree embedding, with the principal tree.",
                     "Monocle2 拟时序，显示在其 DDRTree 嵌入上并叠加主树。"),
        monocle3 = c("Monocle3 pseudotime, principal graph learnt on the UMAP.",
                     "Monocle3 拟时序，主图在 UMAP 上学习。"),
        paga = c("Diffusion pseudotime (PAGA/DPT) as colour; no curve is drawn.",
                 "扩散拟时序（PAGA/DPT）以颜色显示；不绘制曲线。"),
        palantir = c("Palantir pseudotime as colour; no curve is drawn.",
                     "Palantir 拟时序以颜色显示；不绘制曲线。"),
        wot = c("WOT transport-map trajectory score from the earliest time point; this is not a pseudotime.",
                "WOT 自最早时间点起的运输映射轨迹分数；它不是拟时序。"))
      na_en <- if (!is.na(res$n_na) && res$n_na > 0) {
        sprintf(" %s of %s cells have no pseudotime (not on any lineage / unreachable from the root).",
                format(res$n_na, big.mark = ","), format(res$n_cells, big.mark = ","))
      } else ""
      na_zh <- if (!is.na(res$n_na) && res$n_na > 0) {
        sprintf("%s 个细胞中有 %s 个没有拟时序（不在任何谱系上或无法从根到达）。",
                format(res$n_cells, big.mark = ","), format(res$n_na, big.mark = ","))
      } else ""
      insight_bar(
        sprintf("%s Root: %s.%s Pseudotime is an ordering, not clock time; treat branches as hypotheses to check with markers.",
                drawn[1], root_en, na_en),
        sprintf("%s根：%s。%s拟时序是一种排序而非真实时间；分支应作为假设，用标志基因加以验证。",
                drawn[2], root_zh, na_zh))
    })

    render_step_plot(output, input, "preview", function() {
      shiny::req(res$done, length(res$pt_cols) > 0 || length(res$fate_cols) > 0)
      sc_trajectory_plot(rv$obj, method = res$method, group_by = res$group_by,
                         pt_cols = res$pt_cols, fate_cols = res$fate_cols)
    }, name = "trajectory")
  })
}
