#' Module 2: Quality control
#'
#' Compute per-cell QC metrics and filter low-quality cells using either
#' adaptive MAD-based thresholds (recommended) or manual cutoffs.
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_qc
NULL

#' @rdname mod_qc
#' @keywords internal
mod_qc_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = list(en = "Quality control", zh = "质量控制"),
    what = list(
      en = "Remove low-quality cells: empty droplets, dying cells, and debris.",
      zh = "去除低质量细胞：空液滴、濒死细胞和碎片。"),
    why  = list(
      en = "Dying cells leak cytoplasmic RNA and show high mitochondrial content;
            empty droplets have very few genes. Keeping them adds noise.",
      zh = "濒死细胞会泄漏胞质 RNA 并表现出高线粒体含量；空液滴则检测到的基因极少。保留它们会引入噪声。"),
    how  = list(
      en = "The recommended <b>MAD</b> method flags cells that are statistical
            outliers for their own dataset (no guessing fixed numbers). Increase
            the MAD multiplier to keep more cells; decrease to be stricter.",
      zh = "推荐的 <b>MAD</b> 方法会标记出相对于自身数据集的统计离群细胞（无需猜测固定数值）。
            增大 MAD 倍数可保留更多细胞；减小则更严格。"),
    read = list(
      en = "Violins and scatter show per-cell gene count, UMI count and
            mitochondrial fraction; cells outside the accepted band get flagged.
            Watch the mitochondrial tail — a long one means many dying cells.",
      zh = "小提琴图与散点图展示每细胞的基因数、UMI 数与线粒体比例；落在接受区间外的细胞会被标记。留意线粒体比例的长尾——长尾意味着大量濒死细胞。"),
    example = list(
      en = "A cell with 40% mitochondrial reads is likely dying and gets flagged;
               a healthy cell (~5%) is kept.",
      zh = "线粒体读段占 40% 的细胞很可能正在濒死，会被标记；健康细胞（约 5%）则被保留。")
  )
  controls <- shiny::tagList(
    label_with_help("Species", "Sets gene-name patterns for mitochondrial/ribosomal/hemoglobin genes.",
                    label_zh = "物种", tip_zh = "设定线粒体/核糖体/血红蛋白基因的基因名匹配模式。"),
    shiny::selectInput(ns("species"), NULL, c("Human" = "human", "Mouse" = "mouse")),
    label_with_help("Threshold method",
                    "MAD = adaptive, data-driven (recommended). Manual = you set fixed cutoffs.",
                    label_zh = "阈值方法",
                    tip_zh = "MAD = 自适应、数据驱动（推荐）。手动 = 由您设定固定的阈值。"),
    shiny::radioButtons(ns("method"), NULL,
                        c("MAD (adaptive)" = "mad", "Manual" = "manual"), inline = TRUE),
    shiny::conditionalPanel(
      sprintf("input['%s'] == 'mad'", ns("method")),
      label_with_help("MAD multiplier (library size / genes)",
                      "Higher = more permissive. 5 is a common default.",
                      label_zh = "MAD 倍数（文库大小 / 基因数）",
                      tip_zh = "越大越宽松。5 是常用的默认值。"),
      shiny::sliderInput(ns("nmad_lib"), NULL, min = 2, max = 8, value = 5, step = 0.5),
      label_with_help("MAD multiplier (mito %)", "Upper-tail only. 3 is common.",
                      label_zh = "MAD 倍数（线粒体 %）", tip_zh = "仅针对上尾。3 是常用值。"),
      shiny::sliderInput(ns("nmad_mt"), NULL, min = 2, max = 8, value = 3, step = 0.5)
    ),
    shiny::conditionalPanel(
      sprintf("input['%s'] == 'manual'", ns("method")),
      shiny::numericInput(ns("min_genes"), i18n("Min genes/cell", "每个细胞最少基因数"), 200, min = 0),
      shiny::numericInput(ns("max_genes"), i18n("Max genes/cell", "每个细胞最多基因数"), 6000, min = 0),
      shiny::numericInput(ns("max_mt"), i18n("Max mito %", "最大线粒体 %"), 15, min = 0, max = 100)
    ),
    run_button(ns("run"), "Compute & filter", "计算并过滤")
  )
  step_container(title = list(en = "Quality control", zh = "质量控制"),
                 subtitle = list(en = "Flag and remove low-quality cells before they add noise.",
                                 zh = "在低质量细胞引入噪声之前将其标记并去除。"),
                 explainer = explainer, controls = controls,
                 summary = shiny::uiOutput(ns("summary")),
                 preview = preview_plot_ui(ns("preview"),
                   guide = list(en = "QC violins and scatter plots will be drawn here.",
                                zh = "运行后，这里将绘制 QC 小提琴图与散点图。"),
                   caption = list(en = "Per-cell QC metrics; flagged cells are highlighted.",
                                  zh = "每细胞 QC 指标；被标记的细胞高亮显示。")))
}

#' @rdname mod_qc
#' @keywords internal
mod_qc_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    res <- shiny::reactiveValues(before = NA, after = NA, keep = NULL, md = NULL)

    # Pre-select the species from gene-name casing when data arrives, so mouse
    # datasets do not silently get human mito/ribo patterns (and 0% readings).
    shiny::observeEvent(rv$obj, {
      shiny::req(rv$obj)
      shiny::updateSelectInput(session, "species", selected = guess_species(rv$obj))
    }, ignoreNULL = TRUE)

    shiny::observeEvent(input$run, {
      shiny::req(rv$obj)
      if (!require_pkgs("Seurat", "QC")) return(NULL)
      obj <- with_progress_notify({
        o <- qc_add_metrics(rv$obj, species = input$species)
        keep <- if (input$method == "mad") {
          qc_mad_keep(o, input$nmad_lib, input$nmad_mt)
        } else {
          qc_manual_keep(o, input$min_genes, input$max_genes, input$max_mt)
        }
        res$before <- ncol(o)
        res$md <- obj_meta(o)
        res$md$keep <- keep
        res$keep <- keep
        o[, keep]
      }, message = "Computing QC and filtering...")
      if (is.null(obj)) return(NULL)
      res$after <- ncol(obj)
      rv$obj <- obj
      mark_done(rv, "qc")
      log_step(log_rv, "QC",
               params = list(method = input$method, species = input$species,
                             nmad_lib = input$nmad_lib, nmad_mt = input$nmad_mt),
               code = "obj <- subset(obj, cells = keep_cells)  # MAD/manual QC")
      shiny::showNotification(i18n(sprintf("QC done: kept %d of %d cells.",
                                           res$after, res$before),
                                   sprintf("质量控制完成：在 %2$d 个细胞中保留了 %1$d 个。",
                                           res$after, res$before)),
                              type = "message")
    })

    output$summary <- shiny::renderUI({
      if (is.na(res$before)) return(shiny::div(class = "omicone-placeholder",
                                               i18n("Set thresholds and click Compute & filter.",
                                                    "设定阈值后点击“计算并过滤”。")))
      removed <- res$before - res$after
      bslib::layout_columns(
        col_widths = c(4, 4, 4),
        stat_tile(i18n("Before", "过滤前"), format(res$before, big.mark = ",")),
        stat_tile(i18n("Kept", "已保留"), format(res$after, big.mark = ",")),
        stat_tile(i18n("Removed", "已去除"), format(removed, big.mark = ","))
      )
    })

    output$preview <- render_preview_plot(function() {
      md <- res$md
      shiny::req(md)
      md$cellid <- rownames(md)
      md$text <- sprintf("nCount=%s\nnGenes=%s\nmito=%.1f%%\n%s",
                         format(md$nCount_RNA, big.mark = ","),
                         format(md$nFeature_RNA, big.mark = ","),
                         ifelse(is.null(md$percent.mt), NA, md$percent.mt),
                         ifelse(md$keep, "kept", "flagged"))
      ggplot2::ggplot(md, ggplot2::aes(x = nCount_RNA, y = percent.mt,
                                       colour = keep)) +
        ggplot2::geom_point(size = 0.5, alpha = 0.6) +
        ggplot2::scale_x_log10() +
        ggplot2::scale_colour_manual(values = c(`TRUE` = "#3b6ea5", `FALSE` = "#c1476b"),
                                     labels = c(`TRUE` = "kept", `FALSE` = "flagged"),
                                     name = NULL) +
        ggplot2::labs(x = "UMI count (log10)", y = "Mitochondrial %",
                      title = "QC: cells kept vs flagged") +
        omicone_theme()
    })
  })
}
