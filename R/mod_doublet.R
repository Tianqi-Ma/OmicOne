#' Module 3: Doublet removal (去双细胞)
#'
#' Detect droplets that captured two cells (doublets) with scDblFinder (default)
#' or DoubletFinder, then either flag them or drop them from the object.
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_doublet
NULL

#' @rdname mod_doublet
#' @keywords internal
mod_doublet_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = list(en = "Doublet removal (去双细胞)", zh = "去除双细胞"),
    what = list(
      en = "Detect and remove <b>doublets</b>: droplets that accidentally captured
            two cells instead of one.",
      zh = "检测并去除<b>双细胞</b>：意外捕获了两个细胞而非一个的液滴。"),
    why  = list(
      en = "A doublet = one droplet captured two cells, masquerading as a fake
            'intermediate' cell type. Left in, they create spurious clusters and
            confuse downstream annotation; remove them.",
      zh = "双细胞 = 一个液滴捕获了两个细胞，伪装成虚假的“中间”细胞类型。若保留，
            它们会产生假的细胞群并干扰下游注释，因此应予去除。"),
    how  = list(
      en = "The recommended <b>scDblFinder</b> method scores every cell for how
            doublet-like it is. Choose whether to just flag doublets (keep all
            cells but label them) or remove them (drop the flagged cells).",
      zh = "推荐的 <b>scDblFinder</b> 方法会为每个细胞打分，评估其像双细胞的程度。
            可选择仅标记双细胞（保留所有细胞但加上标签）或去除它们（丢弃被标记的细胞）。"),
    read = list(
      en = "Every cell gets a doublet score — high means its profile looks like
            two cells averaged. The summary says how many were flagged; expect
            roughly 1% per 1,000 loaded cells.",
      zh = "每个细胞得到一个双细胞分数——越高表示其表达谱越像两个细胞的平均。汇总会给出被标记的数量；经验上约为每 1,000 个上机细胞 1%。"),
    example = list(
      en = "Two cells of different types share a droplet and look like a novel
               'hybrid' population &mdash; scDblFinder flags them so you can drop them.",
      zh = "两个不同类型的细胞共享一个液滴，看起来像一个新的“杂合”细胞群——
               scDblFinder 会将它们标记出来，以便您丢弃。")
  )
  controls <- shiny::tagList(
    label_with_help("Detection method",
                    "scDblFinder = fast one-click default (recommended). DoubletFinder = classic method needing a per-dataset pK sweep.",
                    label_zh = "检测方法",
                    tip_zh = "scDblFinder = 快速的一键默认方法（推荐）。DoubletFinder = 经典方法，需针对每个数据集进行 pK 扫描。"),
    shiny::selectInput(ns("method"), NULL,
                       choices = c("scDblFinder" = "scDblFinder",
                                   "DoubletFinder" = "DoubletFinder"),
                       selected = "scDblFinder"),
    label_with_help("Action",
                    "Flag = keep every cell but label it doublet/singlet. Remove = drop cells classed as doublets.",
                    label_zh = "操作",
                    tip_zh = "标记 = 保留所有细胞，但标注为双细胞/单细胞。去除 = 丢弃被判为双细胞的细胞。"),
    shiny::radioButtons(ns("action"), NULL,
                        c("Flag only" = "flag", "Remove doublets" = "remove"),
                        selected = "remove", inline = TRUE),
    label_with_help("Custom score threshold (optional)",
                    "Leave blank to use the method's own call. Set a value to class cells with doublet_score above it as doublets.",
                    label_zh = "自定义评分阈值（可选）",
                    tip_zh = "留空则使用方法自身的判定。设定数值后，doublet_score 高于该值的细胞将被判为双细胞。"),
    shiny::numericInput(ns("threshold"), NULL, value = NA, min = 0, max = 1, step = 0.05),
    run_button(ns("run"), "Detect doublets", "检测双细胞")
  )
  step_container(title = list(en = "Doublet removal", zh = "去除双细胞"),
                 subtitle = list(en = "Score and remove droplets that captured two cells.",
                                 zh = "为捕获了两个细胞的液滴打分并去除。"),
                 explainer = explainer, controls = controls,
                 summary = shiny::uiOutput(ns("summary")),
                 preview = preview_plot_ui(ns("preview"),
                   guide = list(en = "Doublet scores and the flagged/kept breakdown will be drawn here.",
                                zh = "运行后，这里将绘制双细胞分数与标记/保留构成。"),
                   caption = list(en = "Doublet score distribution and the flagging result.",
                                  zh = "双细胞分数分布与标记结果。")))
}

#' @rdname mod_doublet
#' @keywords internal
mod_doublet_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    res <- shiny::reactiveValues(before = NA, after = NA, n_doublet = NA,
                                 md = NULL, action = NULL)

    shiny::observeEvent(input$run, {
      shiny::req(rv$obj)
      if (!require_pkgs(c("Seurat", "scDblFinder"), "Doublet removal")) return(NULL)
      thr <- if (is.null(input$threshold) || is.na(input$threshold)) NULL else input$threshold
      action <- input$action
      obj <- with_progress_notify({
        o <- run_doublets(rv$obj, method = input$method)
        md <- obj_meta(o)
        # Optional custom score threshold overrides the method's own call
        if (!is.null(thr) && !is.null(md$doublet_score)) {
          cls <- ifelse(md$doublet_score > thr, "doublet", "singlet")
          o[["doublet_class"]] <- cls
        }
        res$before <- ncol(o)
        res$md <- obj_meta(o)
        is_doub <- res$md$doublet_class == "doublet"
        res$n_doublet <- sum(is_doub, na.rm = TRUE)
        if (action == "remove") o[, !is_doub] else o
      }, message = "Detecting doublets...")
      if (is.null(obj)) return(NULL)
      res$after <- ncol(obj)
      res$action <- action
      rv$obj <- obj
      mark_done(rv, "doublet")
      log_step(log_rv, "Doublet removal",
               params = list(method = input$method, action = action,
                             threshold = thr),
               code = if (action == "remove") {
                 "obj <- run_doublets(obj); obj <- obj[, obj$doublet_class == 'singlet']"
               } else {
                 "obj <- run_doublets(obj)  # flag only"
               })
      shiny::showNotification(
        i18n(sprintf("Doublet detection done: %d doublets (%s).",
                     res$n_doublet, if (action == "remove") "removed" else "flagged"),
             sprintf("双细胞检测完成：%d 个双细胞（%s）。",
                     res$n_doublet, if (action == "remove") "已去除" else "已标记")),
        type = "message")
    })

    output$summary <- shiny::renderUI({
      if (is.na(res$before)) return(shiny::div(class = "omicone-placeholder",
                                               i18n("Pick a method and click Detect doublets.",
                                                    "选择一种方法后点击“检测双细胞”。")))
      pct <- if (res$before > 0) 100 * res$n_doublet / res$before else 0
      bslib::layout_columns(
        col_widths = c(4, 4, 4),
        stat_tile(i18n("Cells in", "输入细胞数"), format(res$before, big.mark = ",")),
        stat_tile(i18n("Doublets", "双细胞"), sprintf("%s (%.1f%%)",
                                      format(res$n_doublet, big.mark = ","), pct)),
        stat_tile(if (identical(res$action, "remove")) i18n("Kept", "已保留") else i18n("Cells out", "输出细胞数"),
                  format(res$after, big.mark = ","))
      )
    })

    output$preview <- render_preview_plot(function() {
      md <- res$md
      shiny::req(md)
      shiny::req(!is.null(md$doublet_score))
      md$cellid <- rownames(md)
      md$doublet_class <- ifelse(is.na(md$doublet_class), "singlet", md$doublet_class)
      md$text <- sprintf("score=%.3f\nclass=%s", md$doublet_score, md$doublet_class)
      ggplot2::ggplot(md, ggplot2::aes(x = doublet_score, fill = doublet_class)) +
        ggplot2::geom_histogram(bins = 50, alpha = 0.75, position = "identity") +
        ggplot2::scale_fill_manual(values = c(singlet = "#3b6ea5", doublet = "#c1476b"),
                                   name = NULL) +
        ggplot2::labs(x = "Doublet score", y = "Cells",
                      title = "Doublet score distribution") +
        omicone_theme()
    })
  })
}
