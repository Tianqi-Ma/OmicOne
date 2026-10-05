#' WES module 5: Tumour mutational burden
#'
#' Non-synonymous mutations per megabase of captured sequence, per sample. A
#' biomarker in its own right and the covariate you most often need to adjust
#' for. The capture size has no default: it is the denominator of every value,
#' so the user states it (or uploads the kit's BED).
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_wes_tmb
NULL

#' @rdname mod_wes_tmb
#' @keywords internal
mod_wes_tmb_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = list(en = "Tumour mutational burden", zh = "肿瘤突变负荷 TMB"),
    what = list(
      en = "The number of non-synonymous somatic mutations per megabase of
            sequence your capture kit actually covered, calculated for each
            sample.",
      zh = "每个样本中，每兆碱基（Mb）实际捕获区域内的非同义体细胞突变数量。"),
    why  = list(
      en = "Raw mutation counts are not comparable between cohorts sequenced with
            different kits — dividing by the captured territory is the first step
            towards comparability. Callers, filters and which classes count as
            non-synonymous still differ between pipelines, so TMB values from
            different studies are rarely interchangeable.",
      zh = "不同试剂盒测序的队列之间，原始突变数不可比——除以实际捕获区域大小是走向可比的第一步。但不同流程的变异检出、过滤以及哪些类别计为非同义仍有差异，因此不同研究的 TMB 数值很少能直接互换。"),
    how  = list(
      en = "Enter your kit's <b>capture size</b> in Mb — it is required, because
            it scales every value. Exome kits differ: the TCGA MC3 reference
            used 35.8 Mb (Agilent SureSelect), other kits range from ~35 to over
            60 Mb. Better still, upload the kit's <b>target BED</b>: variants are
            then restricted to the targets and the territory is the BED's size
            after merging overlapping intervals.",
      zh = "请填写所用试剂盒的<b>捕获区域大小</b>（Mb）——这是必填项，因为它按比例影响每一个数值。不同外显子试剂盒差别很大：TCGA MC3 参考采用 35.8 Mb（Agilent SureSelect），其他试剂盒约 35 至 60 Mb 以上。更好的做法是上传试剂盒的<b>目标区域 BED</b>：变异先限定在目标区域内，捕获大小取合并重叠区间后的 BED 总长。"),
    read = list(
      en = "Each point is one sample's non-synonymous mutations divided by the
            capture size, sorted from lowest to highest; the dashed line is the
            cohort median. On the log axis, samples with no mutation are drawn
            as open triangles at the floor. The dotted 10 mut/Mb line is the
            FoundationOne CDx cut-off from KEYNOTE-158 — a panel assay — and is
            not calibrated for WES non-synonymous TMB. The <b>vs TCGA</b> tab
            puts your cohort beside the 33 TCGA MC3 cohorts.",
      zh = "每个点是一个样本的非同义突变数除以捕获区域大小，按从低到高排序；虚线为队列中位数。在对数轴上，没有突变的样本以空心三角画在底部。点线 10 mut/Mb 是 KEYNOTE-158 中 FoundationOne CDx（一种 panel 检测）的阈值，并未针对 WES 非同义 TMB 校准。<b>对比 TCGA</b> 页签把本队列与 33 个 TCGA MC3 队列并排展示。"),
    example = list(
      en = "With the TCGA LAML demo and 35.8 Mb, the median is about 0.25 mut/Mb
               and no sample reaches 10 mut/Mb — leukaemias sit at the bottom of
               the TCGA range.",
      zh = "使用 TCGA LAML 演示数据和 35.8 Mb 时，中位数约为 0.25 mut/Mb，没有样本达到 10 mut/Mb——白血病位于 TCGA 范围的最低端。")
  )
  controls <- shiny::tagList(
    label_with_help("Capture size (Mb) — required",
                    "The target territory of your capture kit, e.g. 35.8 = TCGA MC3 / Agilent SureSelect. Ignored when a BED is uploaded.",
                    label_zh = "捕获区域大小（Mb）——必填",
                    tip_zh = "捕获试剂盒的目标区域大小，例如 35.8 = TCGA MC3 / Agilent SureSelect。上传 BED 时以 BED 为准。"),
    shiny::numericInput(ns("capture"), NULL, value = NA, min = 0.1, max = 3200, step = 0.1),
    shiny::div(class = "omicone-muted",
               i18n("e.g. 35.8 = TCGA/SureSelect", "例如 35.8 = TCGA/SureSelect")),
    label_with_help("Target BED (optional)",
                    "The kit's target intervals (chrom, start, end). Overlaps are merged; chromosome names are matched to the MAF's style.",
                    label_zh = "目标区域 BED（可选）",
                    tip_zh = "试剂盒的目标区间（染色体、起点、终点）。重叠区间会被合并；染色体名会按 MAF 的写法匹配。"),
    shiny::fileInput(ns("bed"), NULL, accept = c(".bed", ".txt", ".tsv", ".gz")),
    shiny::checkboxInput(ns("log"), i18n("Log scale", "对数坐标"), value = TRUE),
    run_button(ns("run"), "Compute TMB", "计算 TMB")
  )
  step_container(id = id, 
    title     = list(en = "Tumour mutational burden", zh = "肿瘤突变负荷 TMB"),
    subtitle  = list(en = "Mutations per captured megabase, benchmarked against 33 TCGA cohorts.",
                     zh = "每捕获兆碱基的突变数，并与 33 个 TCGA 队列对照。"),
    explainer = explainer,
    controls  = controls,
    summary   = shiny::uiOutput(ns("summary")),
    preview   = shiny::tagList(
      shiny::uiOutput(ns("insight")),
      bslib::navset_card_tab(
        bslib::nav_panel(i18n("Distribution", "分布图"),
                         preview_plot_ui(ns("plot"), download = TRUE,
                                         guide = list(en = "The per-sample TMB distribution will be drawn here.",
                                                      zh = "运行后，这里将绘制每样本 TMB 分布图。"),
                                         caption = list(en = "Each point: one sample's non-synonymous mutations per captured Mb, sorted; dashed = median; dotted = 10 mut/Mb panel cut-off (not WES-calibrated); triangles = 0 mutations.",
                                                        zh = "每个点：一个样本每捕获 Mb 的非同义突变数，已排序；虚线＝中位数；点线＝10 mut/Mb panel 阈值（未针对 WES 校准）；三角＝0 突变。"))),
        bslib::nav_panel(i18n("vs TCGA", "对比 TCGA"),
                         preview_plot_ui(ns("tcga"), download = TRUE,
                                         guide = list(en = "Your cohort beside the 33 TCGA MC3 cohorts will be drawn here.",
                                                      zh = "运行后，这里将把本队列与 33 个 TCGA MC3 队列并排绘制。"),
                                         caption = list(en = "One column per cohort, samples sorted; red bar = cohort median; top numbers = samples drawn. Zero-mutation samples are left out of every column (tcgaCompare rm_zero = TRUE).",
                                                        zh = "每列一个队列，样本已排序；红线＝队列中位数；顶部数字＝绘制的样本数。所有列都不含 0 突变样本（tcgaCompare 的 rm_zero = TRUE）。"))),
        bslib::nav_panel(i18n("Per sample", "各样本"), shiny::uiOutput(ns("tbl_slot")))
      ))
  )
}

#' @rdname mod_wes_tmb
#' @keywords internal
mod_wes_tmb_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    ns  <- session$ns
    res <- step_results(rv, "wes", out = NULL, log = TRUE, bed_name = NULL,
                        bed_regions = NA_integer_)

    shiny::observeEvent(input$run, {
      shiny::req(rv$maf)
      if (!require_pkgs("maftools", "TMB")) return(NULL)
      log_scale <- isTRUE(input$log)
      cap <- suppressWarnings(as.numeric(input$capture %||% NA_real_))
      bed <- NULL
      bed_name <- NULL
      if (!is.null(input$bed)) {
        bed <- tryCatch(wes_read_bed(input$bed$datapath), error = function(e) {
          wes_notify(paste("BED:", conditionMessage(e)), paste("BED：", conditionMessage(e)))
          NULL
        })
        if (is.null(bed)) return(NULL)
        bed_name <- input$bed$name
      } else if (length(cap) != 1 || !is.finite(cap) || cap <= 0) {
        wes_notify("Enter the capture size in Mb (e.g. 35.8 = TCGA/SureSelect), or upload the kit's BED.",
                   "请填写捕获区域大小（Mb，例如 35.8 = TCGA/SureSelect），或上传试剂盒的 BED。",
                   duration = 12)
        return(NULL)
      }
      out <- with_progress_notify(
        wes_tmb(rv$maf, capture_size = cap, bed = bed, log_scale = log_scale, samples = rv$wes_sequenced),
        message = "Computing TMB...")
      if (is.null(out)) return(NULL)
      res$out <- out
      res$log <- log_scale
      res$bed_name <- bed_name
      res$bed_regions <- if (!is.null(bed)) nrow(bed) else NA_integer_
      mark_done(rv, "wes_tmb")
      log_step(log_rv, "WES TMB",
               params = list(captureSize = round(out$capture, 4), bed = bed_name %||% "(none)",
                             logScale = log_scale),
               code = wes_tmb_code(capture = out$capture, bed_name = bed_name,
                                   maf_prefix = wes_chr_has_prefix(rv$maf@data$Chromosome),
                                   log_scale = log_scale,
                                   cohort = wes_tcga_label(rv$maf_source)))
    })

    stats <- shiny::reactive({
      shiny::req(res$out)
      wes_tmb_stats(res$out$df)
    })

    output$summary <- shiny::renderUI({
      if (is.null(rv$maf)) return(wes_no_maf())
      if (is.null(res$out)) {
        return(wes_prompt("Enter the capture size (or upload a BED) and click <b>Compute TMB</b>.",
                          "填写捕获区域大小（或上传 BED）后点击<b>计算 TMB</b>。"))
      }
      s <- stats()
      shiny::tagList(
        stat_tile(i18n("Samples", "样本数"), wes_fmt(s$n)),
        stat_tile(i18n("Median TMB (mut/Mb)", "中位 TMB（mut/Mb）"), wes_fmt(s$median, 2)),
        stat_tile(i18n("≥ 10 mut/Mb", "≥ 10 mut/Mb"), sprintf("%d / %d", s$n_high, s$n)),
        stat_tile(i18n("Capture (Mb)", "捕获（Mb）"),
                  paste0(format(signif(res$out$capture, 4)),
                         if (!is.null(res$bed_name)) " (BED)" else ""))
      )
    })

    output$insight <- shiny::renderUI({
      if (is.null(res$out)) return(NULL)
      s <- stats()
      cap <- format(signif(res$out$capture, 4))
      zero_en <- if (s$n_zero) sprintf(" %d sample(s) have no non-synonymous mutation%s; they count as 0 here but are left out of the vs TCGA tab.",
                                       s$n_zero, if (!is.null(res$bed_name)) " inside the targets" else "")
                 else ""
      zero_zh <- if (s$n_zero) sprintf("%d 个样本%s没有非同义突变；此处按 0 计入，但在对比 TCGA 页签中不纳入。",
                                       s$n_zero, if (!is.null(res$bed_name)) "在目标区域内" else "")
                 else ""
      insight_bar(
        sprintf("Median <b>%.2f</b> mut/Mb (range %.2f–%.2f; non-synonymous mutations / %s Mb); <b>%d of %d</b> samples are ≥ 10 mut/Mb. That cut-off comes from FoundationOne CDx in KEYNOTE-158 and is not calibrated for WES non-synonymous TMB.%s",
                s$median, s$min, s$max, cap, s$n_high, s$n, zero_en),
        sprintf("中位数 <b>%.2f</b> mut/Mb（范围 %.2f–%.2f；非同义突变数 / %s Mb）；<b>%d / %d</b> 个样本 ≥ 10 mut/Mb。该阈值来自 KEYNOTE-158 中的 FoundationOne CDx，并未针对 WES 非同义 TMB 校准。%s",
                s$median, s$min, s$max, cap, s$n_high, s$n, zero_zh))
    })

    tmb_gg <- function() {
      shiny::req(res$out)
      wes_tmb_plot(res$out$df, res$out$capture, log_scale = res$log)
    }
    render_step_plot(output, input, "plot", tmb_gg, name = "wes_tmb", width = 10, height = 6)

    draw_tcga <- function() {
      shiny::req(rv$maf, res$out)
      wes_tcga_gg(res$out$maf, res$out$capture, cohort = wes_tcga_label(rv$maf_source), log_scale = res$log)
    }
    render_step_plot(output, input, "tcga", draw_tcga, name = "wes_tmb_vs_tcga", width = 13, height = 7)

    view <- shiny::reactive({
      shiny::req(res$out)
      df <- res$out$df
      df$total_perMB <- round(df$total_perMB, 4)
      names(df)[names(df) == "total"] <- "nonsyn_mutations"
      names(df)[names(df) == "total_perMB"] <- "TMB_per_Mb"
      df
    })
    tb <- wes_table(ns, "tbl", view, "wes_tmb_per_sample",
                    has_maf = function() !is.null(rv$maf),
                    ready = function() !is.null(res$out))
    output$tbl_slot <- tb$slot
    output$tbl <- tb$table
    output$tbl_dl <- tb$download
  })
}
