#' WES module 5: Tumour mutational burden
#'
#' Mutations per megabase of captured sequence, per sample. A biomarker in its
#' own right (high TMB predicts immunotherapy response in several tumour types)
#' and the covariate you most often need to adjust for.
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
      en = "The number of somatic mutations per megabase of sequence your capture
            kit actually covered, calculated for each sample.",
      zh = "每个样本中，每兆碱基（Mb）捕获区域内的体细胞突变数量。"),
    why  = list(
      en = "Raw mutation counts are not comparable between cohorts sequenced with
            different panels — dividing by the captured size is what makes them
            comparable. High TMB is an approved biomarker for checkpoint
            inhibitors in several tumour types.",
      zh = "不同 panel 测序的队列之间，原始突变数不可比——除以捕获区域大小才能比较。在多种肿瘤中，高 TMB 是免疫检查点抑制剂已获批的生物标志物。"),
    how  = list(
      en = "Set <b>capture size</b> to your kit's actual target size in Mb — the
            default 50 Mb is the usual whole-exome figure. Getting this wrong
            scales every value, so check your kit's documentation. The
            <b>vs TCGA</b> tab puts your cohort's median next to all 33 TCGA
            cohorts, so you can tell a genuinely high burden from an ordinary one.",
      zh = "把<b>捕获区域大小</b>设为你所用试剂盒的实际目标区域（Mb）——默认 50 Mb 是全外显子常见值。设错会让所有数值等比例偏移，请查阅试剂盒文档。<b>vs TCGA</b> 页签把本队列的中位 TMB 与全部 33 个 TCGA 队列并列展示，可以据此判断突变负荷是真的偏高还是普通水平。"),
    read = list(
      en = "Each point is a sample's non-synonymous mutations divided by the
            capture size. The y-axis is log-scaled: small visual gaps are big
            fold differences. As a rule of thumb, ~10 mut/Mb is the cutoff
            clinics use for likely immunotherapy response — most leukaemias sit
            far below it. The <b>vs TCGA</b> tab puts your cohort's median
            beside all 33 TCGA cohorts, so you can see whether it runs hot or
            cold for its tissue type.",
      zh = "每个点是一个样本的非同义突变数除以捕获区域大小。纵轴为对数刻度：图上一小段距离代表很大的倍数差异。经验上，约 10 mut/Mb 是临床用于预测免疫治疗响应的阈值——多数白血病远低于此。<b>对比 TCGA</b> 页签把本队列的中位数与全部 33 个 TCGA 队列并排，一眼看出它在同类肿瘤中偏高还是偏低。"),
    example = list(
      en = "Agilent SureSelect V6 covers ~60 Mb; IDT xGen Exome ~39 Mb; a
               targeted 500-gene panel might be ~1.5 Mb.",
      zh = "Agilent SureSelect V6 约覆盖 60 Mb；IDT xGen Exome 约 39 Mb；500 基因的靶向 panel 可能只有约 1.5 Mb。")
  )
  controls <- shiny::tagList(
    label_with_help("Capture size (Mb)",
                    "The target territory of your capture kit. 50 Mb is the conventional whole-exome default.",
                    label_zh = "捕获区域大小（Mb）",
                    tip_zh = "捕获试剂盒的目标区域大小。全外显子通常按 50 Mb 计算。"),
    shiny::numericInput(ns("capture"), NULL, value = 50, min = 0.1, max = 3200, step = 1),
    shiny::checkboxInput(ns("log"), i18n("Log scale", "对数坐标"), value = TRUE),
    run_button(ns("run"), "Compute TMB", "计算 TMB")
  )
  step_container(
    title     = list(en = "Tumour mutational burden", zh = "肿瘤突变负荷 TMB"),
    subtitle  = list(en = "Mutations per megabase, benchmarked against 33 TCGA cohorts.",
                     zh = "每兆碱基的突变数，并与 33 个 TCGA 队列对照。"),
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
                                       caption = list(en = "Each point: one sample's non-synonymous mutations per captured megabase (log axis).",
                                                      zh = "每个点：一个样本每捕获兆碱基的非同义突变数（对数轴）。"))),
      bslib::nav_panel(i18n("vs TCGA", "对比 TCGA"),  preview_plot_ui(ns("tcga"), download = TRUE)),
      bslib::nav_panel(i18n("Per sample", "各样本"),   shiny::uiOutput(ns("tbl_slot")))
    ))
  )
}

#' @rdname mod_wes_tmb
#' @keywords internal
mod_wes_tmb_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    ns  <- session$ns
    res <- shiny::reactiveValues(df = NULL, capture = NA_real_, log = TRUE)

    shiny::observeEvent(input$run, {
      shiny::req(rv$maf)
      if (!require_pkgs("maftools", "TMB")) return(NULL)
      cap <- suppressWarnings(as.numeric(input$capture %||% NA_real_))
      if (length(cap) != 1 || is.na(cap) || cap <= 0) {
        shiny::showNotification("Capture size must be a positive number of megabases.",
                                type = "error", duration = 10)
        return(NULL)
      }
      # maftools::tmb() draws as a side effect; capture the numbers here and
      # redraw in the plot output so the tab does not depend on run order.
      df <- with_progress_notify({
        grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
        wes_tmb(rv$maf, capture_size = cap, log_scale = isTRUE(input$log))
      }, message = "Computing TMB...")
      if (is.null(df)) return(NULL)
      res$df <- df; res$capture <- cap; res$log <- isTRUE(input$log)
      mark_done(rv, "wes_tmb")
      log_step(log_rv, "WES TMB",
               params = list(captureSize = cap, logScale = input$log),
               code = sprintf('tmb <- maftools::tmb(maf, captureSize = %s, logScale = %s)',
                              cap, input$log))
    })

    tmb_col <- function(df) {
      hit <- grep("per_MB|perMB", colnames(df), ignore.case = TRUE, value = TRUE)
      if (length(hit)) hit[1] else NULL
    }

    output$summary <- shiny::renderUI({
      if (is.null(rv$maf)) return(wes_no_maf())
      df <- res$df
      if (is.null(df)) {
        return(wes_prompt("Set the capture size and click <b>Compute TMB</b>.",
                          "设置捕获区域大小后点击<b>计算 TMB</b>。"))
      }
      cl <- tmb_col(df)
      v  <- if (!is.null(cl)) df[[cl]] else NA_real_
      bslib::layout_columns(
        col_widths = c(3, 3, 3, 3),
        stat_tile(i18n("Samples", "样本数"), format(nrow(df), big.mark = ",")),
        stat_tile(i18n("Median TMB", "中位 TMB"),
                  if (all(is.na(v))) "-" else sprintf("%.2f", stats::median(v, na.rm = TRUE))),
        stat_tile(i18n("Max TMB", "最高 TMB"),
                  if (all(is.na(v))) "-" else sprintf("%.2f", max(v, na.rm = TRUE))),
        stat_tile(i18n("Capture (Mb)", "捕获 (Mb)"), res$capture)
      )
    })

    output$insight <- shiny::renderUI({
      df <- res$df
      if (is.null(df)) return(NULL)
      cl <- tmb_col(df)
      v <- if (!is.null(cl)) df[[cl]] else NA_real_
      if (all(is.na(v))) return(NULL)
      med <- stats::median(v, na.rm = TRUE)
      cap <- if (is.na(res$capture)) "?" else format(res$capture)
      verdict_en <- if (isTRUE(med >= 10)) "above the ~10 mut/Mb line clinics use for likely immunotherapy response"
                    else "below the ~10 mut/Mb line clinics use for likely immunotherapy response"
      verdict_zh <- if (isTRUE(med >= 10)) "高于临床用于预测免疫治疗响应的 ~10 mut/Mb 参考线"
                    else "低于临床用于预测免疫治疗响应的 ~10 mut/Mb 参考线"
      insight_bar(
        sprintf("Median TMB <b>%.2f</b> mut/Mb (range %.2f–%.2f, capture %s Mb) — %s. The <b>vs TCGA</b> tab shows where this sits among the 33 TCGA cohorts.",
                med, min(v, na.rm = TRUE), max(v, na.rm = TRUE), cap, verdict_en),
        sprintf("中位 TMB <b>%.2f</b> mut/Mb（范围 %.2f–%.2f，捕获 %s Mb）——%s。<b>对比 TCGA</b> 页签可看它处于 33 个 TCGA 队列中的什么位置。",
                med, min(v, na.rm = TRUE), max(v, na.rm = TRUE), cap, verdict_zh))
    })

    draw_tmb <- with_text_boost(function() {
      shiny::req(rv$maf, res$df)
      maftools::tmb(maf = rv$maf, captureSize = res$capture, logScale = res$log)
    })
    output$plot <- render_base_plot(draw_tmb)
    register_figure_download(output, input, "plot", draw_tmb, "wes_tmb",
                             width = 10, height = 7)

    draw_tcga <- with_text_boost(function() {
      shiny::req(rv$maf, res$df)
      lab <- rv$maf_source %||% "This cohort"
      lab <- sub("\\.(maf|maf\\.gz|txt|tsv|csv)$", "", basename(lab),
                 ignore.case = TRUE)
      lab <- trimws(gsub("[()]", "", sub("(?i)\\bdemo\\b\\s*:?", "", lab,
                                         perl = TRUE)))
      if (!nzchar(lab)) lab <- "This cohort"
      if (nchar(lab) > 18) lab <- paste0(substr(lab, 1, 17), "~")
      maftools::tcgaCompare(maf = rv$maf, cohortName = lab,
                            capture_size = res$capture, logscale = res$log,
                            cohortFontSize = 1.0, axisFontSize = 1.15)
    })
    output$tcga <- render_base_plot(draw_tcga)
    register_figure_download(output, input, "tcga", draw_tcga,
                             "wes_tmb_vs_tcga", width = 12, height = 8)

    output$tbl_slot <- shiny::renderUI({
      if (is.null(res$df)) return(wes_no_maf())
      if (has_pkg("DT")) DT::dataTableOutput(ns("tbl"))
      else shiny::verbatimTextOutput(ns("tbl_txt"))
    })
    view <- shiny::reactive({
      df <- res$df; shiny::req(df)
      cl <- tmb_col(df)
      if (!is.null(cl)) df[[cl]] <- round(df[[cl]], 3)
      df[order(-df[[cl %||% colnames(df)[ncol(df)]]]), , drop = FALSE]
    })
    if (has_pkg("DT")) {
      output$tbl <- DT::renderDataTable(
        DT::datatable(view(), rownames = FALSE, filter = "top",
                      options = list(pageLength = 15, scrollX = TRUE)))
    } else {
      output$tbl_txt <- shiny::renderPrint(utils::head(view(), 20))
    }
  })
}
