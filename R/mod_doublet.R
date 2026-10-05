#' Module 3: Doublet removal (去双细胞)
#'
#' Detect droplets that captured two cells (doublets) with scDblFinder, run
#' within each sample, then either flag them or drop them from the object.
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_doublet
NULL

#' @rdname mod_doublet
#' @keywords internal
mod_doublet_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = list(en = "Doublet removal", zh = "去除双细胞"),
    what = list(
      en = "Detect and remove <b>doublets</b>: droplets that accidentally captured
            two cells instead of one.",
      zh = "检测并去除<b>双细胞</b>：意外捕获了两个细胞而非一个的液滴。"),
    why  = list(
      en = "A doublet masquerades as a fake 'intermediate' cell type. Left in,
            doublets create spurious clusters and confuse downstream
            annotation.",
      zh = "双细胞会伪装成虚假的“中间”细胞类型。若保留，它们会产生假的细胞群并干扰下游注释。"),
    how  = list(
      en = "<b>scDblFinder</b> simulates artificial doublets and scores every cell
            by how much it resembles them. Doublets only form within one
            sample (one channel), so choose the sample column: each sample is
            then scored on its own. Choose whether to flag doublets (keep all
            cells, labelled) or remove them.",
      zh = "<b>scDblFinder</b> 会模拟人工双细胞，并按每个细胞与它们的相似程度打分。双细胞只会在同一样本（同一通道）内形成，因此请选择样本列，让每个样本单独打分。可选择仅标记（保留所有细胞并加标签）或去除。"),
    read = list(
      en = "Histogram of doublet scores, one panel per sample: blue = called
            singlet, red = called doublet; a dashed line marks a custom
            threshold. scDblFinder's prior is about 0.8% doublets per 1,000
            cells recovered in a sample (8% at 10,000 cells).",
      zh = "双细胞分数直方图，每个样本一个面板：蓝色 = 判为单细胞，红色 = 判为双细胞；虚线为自定义阈值。scDblFinder 的先验约为每个样本每回收 1,000 个细胞 0.8% 双细胞（10,000 个细胞时约 8%）。"),
    example = list(
      en = "Two samples of 5,000 cells each: run per sample, the prior is ~4%
            each; pooled as 10,000 cells it would be ~8%, and cross-sample
            doublets that cannot exist would be simulated.",
      zh = "两个各 5,000 个细胞的样本：按样本运行时每个样本的先验约为 4%；若合并为 10,000 个细胞运行则约为 8%，且会模拟出物理上不存在的跨样本双细胞。")
  )
  controls <- shiny::tagList(
    shiny::div(class = "omicone-note",
               i18n("Method: scDblFinder (Germain et al. 2021), seed 42.",
                    "方法：scDblFinder（Germain 等，2021），随机种子 42。")),
    label_with_help("Sample column",
                    "Doublets are detected within each capture. Pick the column naming the 10x channel (lane) of each cell; if several donors were pooled in one channel, use the channel, not the donor. Choose single sample only if all cells come from one channel.",
                    label_zh = "样本列",
                    tip_zh = "双细胞在每次捕获内检测。请选择标识每个细胞所属 10x 通道（lane）的列；若多个供体混合在同一通道上机，应选通道而非供体。只有当所有细胞来自同一通道时才选单样本。"),
    shiny::uiOutput(ns("sample_ui")),
    label_with_help("Action",
                    "Flag = keep every cell but label it doublet/singlet. Remove = drop cells classed as doublets.",
                    label_zh = "操作",
                    tip_zh = "标记 = 保留所有细胞，但标注为双细胞/单细胞。去除 = 丢弃被判为双细胞的细胞。"),
    shiny::radioButtons(ns("action"), NULL,
                        c("Flag only" = "flag", "Remove doublets" = "remove"),
                        selected = "remove", inline = TRUE),
    label_with_help("Custom score threshold (optional)",
                    "Leave blank to use scDblFinder's own call. Set a value (0-1) to class cells with doublet_score above it as doublets.",
                    label_zh = "自定义评分阈值（可选）",
                    tip_zh = "留空则使用 scDblFinder 自身的判定。设定数值（0-1）后，doublet_score 高于该值的细胞将被判为双细胞。"),
    shiny::numericInput(ns("threshold"), NULL, value = NA, min = 0, max = 1, step = 0.05),
    run_button(ns("run"), "Detect doublets", "检测双细胞")
  )
  step_container(title = list(en = "Doublet removal", zh = "去除双细胞"),
                 subtitle = list(en = "Score and remove droplets that captured two cells, per sample.",
                                 zh = "按样本为捕获了两个细胞的液滴打分并去除。"),
                 explainer = explainer, controls = controls,
                 summary = shiny::uiOutput(ns("summary")),
                 preview = shiny::tagList(
                   shiny::uiOutput(ns("insight")),
                   preview_plot_ui(ns("preview"), download = TRUE,
                     guide = list(en = "Doublet-score histograms per sample will be drawn here.",
                                  zh = "运行后，这里将绘制每个样本的双细胞分数直方图。"),
                     caption = list(en = "Cells per doublet-score bin; blue = singlet call, red = doublet call.",
                                    zh = "每个分数区间的细胞数；蓝色 = 判为单细胞，红色 = 判为双细胞。"))))
}

#' @rdname mod_doublet
#' @keywords internal
mod_doublet_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    res <- step_results(rv, "sc", before = NA, after = NA, n_doublet = NA,
                        md = NULL, action = NULL, samples = NULL, threshold = NULL)

    output$sample_ui <- shiny::renderUI({
      md <- obj_meta(rv$obj)
      cols <- sc_group_cols(md)
      choices <- c(stats::setNames(".single", "(single sample)"), stats::setNames(cols, cols))
      default <- guess_batch_col(md[, cols, drop = FALSE]) %||% ".single"
      shiny::selectInput(session$ns("samples"), NULL, choices = choices,
                         selected = keep_selected(shiny::isolate(input$samples), choices, default))
    })

    shiny::observeEvent(input$run, {
      shiny::req(rv$obj)
      if (!require_pkgs(c("Seurat", "scDblFinder", "SingleCellExperiment"),
                        "Doublet removal")) return(NULL)
      thr <- num_input(input$threshold, 0, 1)
      if (is.na(thr)) thr <- NULL
      action <- input$action
      samples <- input$samples %||% ".single"
      samples <- if (samples %in% obj_meta_cols(rv$obj)) samples else NULL
      seed <- 42
      out <- with_progress_notify({
        # re-runs start from the same input, so removals never compound
        o <- run_doublets(step_input(rv, "doublet", rv$obj), samples = samples,
                          seed = seed, threshold = thr)
        md <- obj_meta(o)
        is_doub <- md$doublet_class == "doublet"
        is_doub[is.na(is_doub)] <- FALSE
        list(obj = if (action == "remove") o[, !is_doub] else o, md = md,
             n_doublet = sum(is_doub))
      }, message = "Detecting doublets...")
      if (is.null(out)) return(NULL)
      keep_cols <- intersect(c("doublet_score", "doublet_class", samples), names(out$md))
      res$before <- nrow(out$md)
      res$after <- ncol(out$obj)
      res$md <- out$md[, keep_cols, drop = FALSE]
      res$n_doublet <- out$n_doublet
      res$action <- action
      res$samples <- samples
      res$threshold <- thr
      rv$obj <- out$obj
      mark_done(rv, "doublet")
      log_step(log_rv, "Doublet removal",
               params = list(method = "scDblFinder", samples = samples, seed = seed,
                             action = action, threshold = thr),
               code = doublet_log_code(samples, seed, thr, action))
      shiny::showNotification(
        i18n(sprintf("Doublet detection done: %s doublets (%s).",
                     format(res$n_doublet, big.mark = ","),
                     if (action == "remove") "removed" else "flagged"),
             sprintf("双细胞检测完成：%s 个双细胞（%s）。",
                     format(res$n_doublet, big.mark = ","),
                     if (action == "remove") "已去除" else "已标记")),
        type = "message")
    })

    output$summary <- shiny::renderUI({
      if (is.na(res$before)) return(shiny::div(class = "omicone-placeholder",
                                               i18n("Pick the sample column and click Detect doublets.",
                                                    "选择样本列后点击“检测双细胞”。")))
      pct <- if (res$before > 0) 100 * res$n_doublet / res$before else 0
      shiny::tagList(
        stat_tile(i18n("Cells in", "输入细胞数"), format(res$before, big.mark = ",")),
        stat_tile(i18n("Doublets", "双细胞"), sprintf("%s (%.1f%%)",
                                                      format(res$n_doublet, big.mark = ","), pct)),
        stat_tile(if (identical(res$action, "remove")) i18n("Kept", "已保留")
                  else i18n("Cells out", "输出细胞数"),
                  format(res$after, big.mark = ",")),
        stat_tile(i18n("Run per", "运行单位"),
                  if (is.null(res$samples)) i18n("single sample", "单样本") else res$samples)
      )
    })

    output$insight <- shiny::renderUI({
      md <- res$md
      if (is.null(md)) return(NULL)
      grp <- if (is.null(res$samples)) rep("all", nrow(md)) else as.character(md[[res$samples]])
      n <- tapply(rep(1, nrow(md)), grp, sum)
      obs <- tapply(md$doublet_class == "doublet", grp, mean, na.rm = TRUE)
      expct <- doublet_expected_rate(n)
      cav_en <- if (!is.null(res$threshold)) {
        sprintf(" Calls use your score cut-off %s, not scDblFinder's own threshold.", res$threshold)
      } else ""
      cav_zh <- if (!is.null(res$threshold)) {
        sprintf("判定使用你设定的分数阈值 %s，而非 scDblFinder 自身的阈值。", res$threshold)
      } else ""
      insight_bar(
        sprintf("%s of %s cells called doublets (%.1f%%); per sample %.1f%%-%.1f%% observed vs %.1f%%-%.1f%% expected from scDblFinder's prior (0.8%% per 1,000 cells in the sample).%s",
                format(res$n_doublet, big.mark = ","), format(res$before, big.mark = ","),
                100 * res$n_doublet / max(1, res$before), 100 * min(obs), 100 * max(obs),
                100 * min(expct), 100 * max(expct), cav_en),
        sprintf("%s 个细胞中有 %s 个被判为双细胞（%.1f%%）；各样本实际 %.1f%%-%.1f%%，按 scDblFinder 先验（样本内每 1,000 个细胞 0.8%%）预期 %.1f%%-%.1f%%。%s",
                format(res$before, big.mark = ","), format(res$n_doublet, big.mark = ","),
                100 * res$n_doublet / max(1, res$before), 100 * min(obs), 100 * max(obs),
                100 * min(expct), 100 * max(expct), cav_zh))
    })

    render_step_plot(output, input, "preview", function() {
      md <- res$md
      shiny::req(md, !is.null(md$doublet_score))
      doublet_plot(doublet_plot_data(md, res$samples), res$threshold)
    }, name = "doublet")
  })
}
