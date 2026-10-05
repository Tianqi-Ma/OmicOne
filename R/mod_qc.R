#' Module 2: Quality control
#'
#' Compute per-cell QC metrics and filter low-quality cells using either
#' adaptive MAD-based thresholds (recommended; per sample when a sample column
#' exists) or manual cutoffs.
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
      en = "Remove low-quality cells: low-UMI droplets and debris, dying cells, and
            unusually large libraries.",
      zh = "去除低质量细胞：低 UMI 的液滴与碎片、濒死细胞，以及文库异常大的细胞。"),
    why  = list(
      en = "Dying cells leak cytoplasmic RNA and show high mitochondrial content;
            debris and barely captured droplets have very few genes. Keeping
            them adds noise. This step does not test for empty droplets
            (emptyDrops); it assumes Cell Ranger's filtered matrix.",
      zh = "濒死细胞会泄漏胞质 RNA 并表现出高线粒体含量；碎片和捕获很差的液滴检测到的基因极少。保留它们会引入噪声。本步骤不做空液滴检验（emptyDrops），默认输入已是 Cell Ranger 过滤后的矩阵。"),
    how  = list(
      en = "The recommended <b>MAD</b> method flags cells that are outliers for
            their own sample: more than <i>n</i> median absolute deviations
            from the median of log UMIs or log genes, or above the median
            mitochondrial % by the mito multiplier. Pick the sample column so
            each sample gets its own thresholds; increase a multiplier to keep
            more cells.",
      zh = "推荐的 <b>MAD</b> 方法标记相对于其所在样本的离群细胞：log UMI 或 log 基因数偏离中位数超过 <i>n</i> 个绝对中位差，或线粒体比例高出中位数超过设定倍数。选择样本列可让每个样本使用各自的阈值；增大倍数可保留更多细胞。"),
    read = list(
      en = "Three panels: genes per cell, UMIs per cell (both log10) and
            mitochondrial %. Each point is a cell (blue kept, red flagged);
            dashed lines are the thresholds actually applied, one set per
            sample. A long red tail in the mitochondrial panel means many
            dying cells.",
      zh = "三个面板：每细胞基因数、每细胞 UMI（均为 log10）和线粒体比例。每个点是一个细胞（蓝色保留，红色标记）；虚线是实际使用的阈值，每个样本一组。线粒体面板中红色长尾意味着大量濒死细胞。"),
    example = list(
      en = "In a sample whose median mitochondrial fraction is 4% with a MAD of
            1.5%, a 3-MAD cut flags cells above 8.5%; a cell at 40% is flagged.",
      zh = "某样本线粒体比例中位数为 4%、绝对中位差为 1.5%，按 3 倍 MAD 截断会标记高于 8.5% 的细胞；比例为 40% 的细胞会被标记。")
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
      label_with_help("Thresholds per sample",
                      "Compute medians and MADs within each sample (recommended when samples differ in depth). Only categorical columns without missing values are offered.",
                      label_zh = "按样本计算阈值",
                      tip_zh = "在每个样本内分别计算中位数与绝对中位差（样本测序深度不同时推荐）。仅列出无缺失值的分类列。"),
      shiny::uiOutput(ns("batch_ui")),
      label_with_help("MAD multiplier (library size / genes)",
                      "Higher = more permissive. 5 is a common default.",
                      label_zh = "MAD 倍数（文库大小 / 基因数）",
                      tip_zh = "越大越宽松。5 是常用的默认值。"),
      shiny::sliderInput(ns("nmad_lib"), NULL, min = 2, max = 8, value = 5, step = 0.5),
      label_with_help("MAD multiplier (mito %)", "Upper-tail only. 3 is common.",
                      label_zh = "MAD 倍数（线粒体 %）", tip_zh = "仅针对上尾。3 是常用值。"),
      shiny::sliderInput(ns("nmad_mt"), NULL, min = 2, max = 8, value = 3, step = 0.5),
      label_with_help("Hard mito cap % (optional)",
                      "Also drop cells above this mitochondrial %, whatever the MAD says. Useful when most cells read 0% (MAD = 0, so no adaptive mito cut is possible). Leave blank for none.",
                      label_zh = "线粒体硬上限 %（可选）",
                      tip_zh = "无论 MAD 结果如何，线粒体比例高于该值的细胞也会被去除。当大多数细胞为 0%（MAD = 0，无法自适应截断）时有用。留空表示不设上限。"),
      shiny::numericInput(ns("max_mt_cap"), NULL, value = NA, min = 0, max = 100, step = 1)
    ),
    shiny::conditionalPanel(
      sprintf("input['%s'] == 'manual'", ns("method")),
      shiny::numericInput(ns("min_genes"), i18n("Min genes/cell", "每个细胞最少基因数"), 200, min = 0),
      shiny::numericInput(ns("max_genes"), i18n("Max genes/cell", "每个细胞最多基因数"), 6000, min = 0),
      shiny::numericInput(ns("max_mt"), i18n("Max mito %", "最大线粒体 %"), 15, min = 0, max = 100)
    ),
    run_button(ns("run"), "Compute & filter", "计算并过滤")
  )
  step_container(id = id, title = list(en = "Quality control", zh = "质量控制"),
                 subtitle = list(en = "Flag and remove low-quality cells before they add noise.",
                                 zh = "在低质量细胞引入噪声之前将其标记并去除。"),
                 explainer = explainer, controls = controls,
                 summary = shiny::uiOutput(ns("summary")),
                 preview = shiny::tagList(
                   shiny::uiOutput(ns("insight")),
                   preview_plot_ui(ns("preview"), download = TRUE,
                     guide = list(en = "Violins of genes, UMIs and mitochondrial % per cell, with the thresholds, will be drawn here.",
                                  zh = "运行后，这里将绘制每细胞基因数、UMI 与线粒体比例的小提琴图及阈值线。"),
                     caption = list(en = "One point = one cell before filtering (blue kept, red flagged); dashed lines = thresholds applied.",
                                    zh = "每个点为过滤前的一个细胞（蓝色保留，红色标记）；虚线为实际使用的阈值。"))))
}

#' @rdname mod_qc
#' @keywords internal
mod_qc_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    res <- step_results(rv, "sc", before = NA, after = NA, keep = NULL, md = NULL,
                        thr = NULL, batch = NULL, method = NULL, id_type = NULL)

    # Pre-select the species from gene names when a new dataset arrives (not
    # after every step, which would overwrite the user's choice), so mouse
    # datasets do not silently get human mito/ribo patterns (and 0% readings).
    shiny::observeEvent(rv$epoch_sc, {
      obj <- shiny::isolate(rv$obj)
      shiny::req(obj)
      shiny::updateSelectInput(session, "species", selected = guess_species(obj))
    }, ignoreNULL = TRUE)

    output$batch_ui <- shiny::renderUI({
      md <- obj_meta(rv$obj)
      cols <- sc_group_cols(md)
      choices <- c(stats::setNames(".all", "(all cells together)"), stats::setNames(cols, cols))
      default <- guess_batch_col(md[, cols, drop = FALSE]) %||% ".all"
      shiny::selectInput(session$ns("batch"), NULL, choices = choices,
                         selected = keep_selected(shiny::isolate(input$batch), choices, default))
    })

    shiny::observeEvent(input$run, {
      shiny::req(rv$obj)
      if (!require_pkgs("Seurat", "QC")) return(NULL)
      species <- input$species
      method <- input$method
      nmad_lib <- input$nmad_lib
      nmad_mt <- input$nmad_mt
      batch <- input$batch %||% ".all"
      batch <- if (method == "mad" && batch %in% obj_meta_cols(rv$obj)) batch else NULL
      cap <- num_input(input$max_mt_cap, 0, 100)
      min_genes <- num_input(input$min_genes, 0)
      max_genes <- num_input(input$max_genes, 0)
      max_mt <- num_input(input$max_mt, 0, 100)
      if (method == "manual" && (is.na(min_genes) || is.na(max_genes) || is.na(max_mt))) {
        shiny::showNotification(i18n("Enter numeric manual thresholds.",
                                     "请输入数值型的手动阈值。"), type = "error")
        return(NULL)
      }
      out <- with_progress_notify({
        # always filter the cells QC first received, never its own output
        o <- qc_add_metrics(step_input(rv, "qc", rv$obj), species = species)
        md <- obj_meta(o)
        keep <- if (method == "mad") {
          qc_mad_keep_from_meta(md, nmad_lib, nmad_mt, batch = batch, max_mt = cap)
        } else {
          qc_manual_keep_from_meta(md, min_genes, max_genes, max_mt)
        }
        list(obj = o[, as.vector(keep)], md = md, keep = keep,
             id_type = feature_id_type(rownames(o)))
      }, message = "Computing QC and filtering...")
      if (is.null(out)) return(NULL)
      thr <- attr(out$keep, "thresholds")
      res$before <- nrow(out$md)
      res$after <- ncol(out$obj)
      res$md <- out$md
      res$keep <- as.vector(out$keep)
      res$thr <- thr
      res$batch <- batch
      res$method <- method
      res$id_type <- out$id_type
      rv$obj <- out$obj
      mark_done(rv, "qc")
      rule <- if (method == "mad") {
        sprintf("MAD thresholds%s: %s MADs on log1p(nCount_RNA) / log1p(nFeature_RNA), %s MADs above the median percent.mt%s",
                if (is.null(batch)) "" else paste0(" per '", batch, "'"), nmad_lib, nmad_mt,
                if (is.na(cap)) "" else sprintf(", hard cap %s%%", cap))
      } else {
        "Manual thresholds"
      }
      log_step(log_rv, "QC",
               params = list(method = method, species = species, batch = batch,
                             nmad_lib = if (method == "mad") nmad_lib else NULL,
                             nmad_mt = if (method == "mad") nmad_mt else NULL,
                             max_mt_cap = if (method == "mad" && !is.na(cap)) cap else NULL,
                             min_genes = if (method == "manual") min_genes else NULL,
                             max_genes = if (method == "manual") max_genes else NULL,
                             max_mt = if (method == "manual") max_mt else NULL),
               code = qc_log_code(species, thr, batch, rule))
      if (identical(out$id_type, "ensembl")) {
        shiny::showNotification(
          i18n("Feature names look like Ensembl IDs: mitochondrial / ribosomal / hemoglobin percentages read 0%. Convert them to gene symbols for a meaningful mito filter.",
               "基因名看起来是 Ensembl ID：线粒体 / 核糖体 / 血红蛋白比例均为 0%。请先转换为基因符号，线粒体过滤才有意义。"),
          type = "warning", duration = 15)
      }
      shiny::showNotification(i18n(sprintf("QC done: kept %s of %s cells.",
                                           format(res$after, big.mark = ","),
                                           format(res$before, big.mark = ",")),
                                   sprintf("质量控制完成：在 %s 个细胞中保留了 %s 个。",
                                           format(res$before, big.mark = ","),
                                           format(res$after, big.mark = ","))),
                              type = "message")
    })

    output$summary <- shiny::renderUI({
      if (is.na(res$before)) return(shiny::div(class = "omicone-placeholder",
                                               i18n("Set thresholds and click Compute & filter.",
                                                    "设定阈值后点击“计算并过滤”。")))
      removed <- res$before - res$after
      pct <- if (res$before > 0) 100 * removed / res$before else 0
      shiny::tagList(
        stat_tile(i18n("Before", "过滤前"), format(res$before, big.mark = ",")),
        stat_tile(i18n("Kept", "已保留"), format(res$after, big.mark = ",")),
        stat_tile(i18n("Removed", "已去除"),
                  sprintf("%s (%.1f%%)", format(removed, big.mark = ","), pct)),
        stat_tile(i18n("Thresholds", "阈值"),
                  if (identical(res$method, "manual")) i18n("manual", "手动")
                  else if (is.null(res$batch)) i18n("MAD, all cells", "MAD，全体细胞")
                  else i18n(sprintf("MAD, per %s (%d)", res$batch, nrow(res$thr)),
                            sprintf("MAD，按 %s（%d 个）", res$batch, nrow(res$thr))))
      )
    })

    output$insight <- shiny::renderUI({
      thr <- res$thr
      if (is.null(thr)) return(NULL)
      pct <- 100 * res$after / max(1, res$before)
      zero_mt <- sum(grepl("percent.mt", thr$mad_zero, fixed = TRUE))
      cav_en <- character(0)
      cav_zh <- character(0)
      if (identical(res$id_type, "ensembl")) {
        cav_en <- c(cav_en, "Feature names are Ensembl IDs, so the mitochondrial filter had no genes to measure.")
        cav_zh <- c(cav_zh, "基因名为 Ensembl ID，线粒体过滤没有可度量的基因。")
      } else if (zero_mt > 0) {
        cav_en <- c(cav_en, sprintf("The mitochondrial MAD was 0 in %d of %d threshold groups, so no adaptive mito cut was applied there (set a hard cap if needed).",
                                    zero_mt, nrow(thr)))
        cav_zh <- c(cav_zh, sprintf("在 %d / %d 个阈值组中线粒体比例的绝对中位差为 0，因此这些组未做自适应线粒体截断（必要时请设置硬上限）。",
                                    zero_mt, nrow(thr)))
      }
      if (!is.null(res$batch) && any(thr$n < 100)) {
        cav_en <- c(cav_en, sprintf("%d sample(s) have fewer than 100 cells; their MADs are imprecise.",
                                    sum(thr$n < 100)))
        cav_zh <- c(cav_zh, sprintf("%d 个样本少于 100 个细胞，其绝对中位差估计不够精确。",
                                    sum(thr$n < 100)))
      }
      insight_bar(
        paste(sprintf("Kept <b>%s</b> of %s cells (%.1f%%).", format(res$after, big.mark = ","),
                      format(res$before, big.mark = ","), pct),
              paste(cav_en, collapse = " ")),
        paste(sprintf("在 %s 个细胞中保留 <b>%s</b> 个（%.1f%%）。", format(res$before, big.mark = ","),
                      format(res$after, big.mark = ","), pct),
              paste(cav_zh, collapse = "")))
    })

    render_step_plot(output, input, "preview", function() {
      md <- res$md
      shiny::req(md, res$thr)
      qc_plot(md, res$keep, res$thr, res$batch)
    }, name = "qc")
  })
}
