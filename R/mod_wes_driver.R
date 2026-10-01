#' WES module 7: Drivers and gene interactions
#'
#' Two questions about which genes matter: which are mutated in a positionally
#' clustered way (the oncodrive signal of a driver), and which pairs of genes
#' tend to be mutated together — or never together — in the same patient.
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_wes_driver
NULL

#' @rdname mod_wes_driver
#' @keywords internal
mod_wes_driver_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = list(en = "Drivers & interactions", zh = "驱动基因与互作"),
    what = list(
      en = "<b>Oncodrive</b> finds genes whose mutations pile up at a few
            positions rather than spreading out. <b>Interactions</b> tests every
            pair of frequently mutated genes for co-occurrence or mutual
            exclusivity.",
      zh = "<b>Oncodrive</b> 找出突变集中在少数位点、而非均匀散布的基因。<b>互作分析</b>对每一对高频突变基因检验共现或互斥。"),
    why  = list(
      en = "A gene can be mutated often just because it is long. Positional
            clustering is evidence of selection, not size. Mutual exclusivity
            suggests two genes hit the same pathway, so one is enough.",
      zh = "一个基因突变频繁，可能只是因为它长。位点聚集是受选择的证据，而不是长度的结果。互斥则提示两个基因作用于同一通路，突变其一即可。"),
    how  = list(
      en = "Raise <b>minimum mutations</b> on a large cohort to cut noise. The
            interaction plot marks significance with dots: the darker the tile,
            the stronger the co-occurrence (green) or exclusivity (brown).",
      zh = "在大队列中提高<b>最小突变数</b>可以降噪。互作图用点标记显著性：格子颜色越深，共现（绿）或互斥（棕）越强。"),
    example = list(
      en = "In AML, <code>NPM1</code> and <code>FLT3</code> co-occur, while
               <code>TP53</code> is mutually exclusive with most of the
               <code>NPM1</code>-driven group.",
      zh = "在 AML 中，<code>NPM1</code> 与 <code>FLT3</code> 共现，而 <code>TP53</code> 与大多数由 <code>NPM1</code> 驱动的分组互斥。")
  )
  controls <- shiny::tagList(
    label_with_help("Minimum mutations per gene",
                    "Genes with fewer mutations than this are not tested by oncodrive.",
                    label_zh = "每基因最少突变数",
                    tip_zh = "突变数少于此值的基因不参与 oncodrive 检验。"),
    shiny::numericInput(ns("min_mut"), NULL, value = 5, min = 2, max = 50, step = 1),
    label_with_help("FDR cutoff", "Drivers at or below this FDR are highlighted.",
                    label_zh = "FDR 阈值", tip_zh = "FDR 不高于此值的驱动基因会被高亮。"),
    shiny::numericInput(ns("fdr"), NULL, value = 0.1, min = 0.001, max = 0.5, step = 0.01),
    shiny::uiOutput(ns("aa_ui")),
    label_with_help("Top genes for interactions",
                    "How many of the most mutated genes to test pairwise.",
                    label_zh = "互作检验的基因数",
                    tip_zh = "取多少个高频突变基因做两两检验。"),
    shiny::numericInput(ns("top"), NULL, value = 25, min = 5, max = 60, step = 1),
    run_button(ns("run"), "Find drivers", "检测驱动基因")
  )
  step_container(
    title     = list(en = "Drivers & interactions", zh = "驱动基因与互作"),
    explainer = explainer,
    controls  = controls,
    summary   = shiny::uiOutput(ns("summary")),
    preview   = bslib::navset_card_tab(
      bslib::nav_panel(i18n("Oncodrive", "Oncodrive"),   preview_plot_ui(ns("drv"), download = TRUE)),
      bslib::nav_panel(i18n("Interactions", "基因互作"), preview_plot_ui(ns("int"), download = TRUE)),
      bslib::nav_panel(i18n("Driver table", "驱动基因表"), shiny::uiOutput(ns("tbl_slot")))
    )
  )
}

#' @rdname mod_wes_driver
#' @keywords internal
mod_wes_driver_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    ns  <- session$ns
    res <- shiny::reactiveValues(drv = NULL, fdr = 0.1, top = 25)

    output$aa_ui <- shiny::renderUI({
      shiny::req(rv$maf)
      f <- wes_fields(rv$maf)
      shiny::tagList(
        label_with_help("Protein change column",
                        "Oncodrive needs amino-acid positions. Auto-detected.",
                        label_zh = "蛋白改变列",
                        tip_zh = "Oncodrive 需要氨基酸位置信息，会自动识别。"),
        shiny::selectInput(ns("aa_col"), NULL,
                           choices = c(stats::setNames("", "(auto)"), f),
                           selected = wes_guess_aa_col(rv$maf) %||% "")
      )
    })

    shiny::observeEvent(input$run, {
      shiny::req(rv$maf)
      if (!require_pkgs("maftools", "Driver detection")) return(NULL)
      aa <- if (nzchar(input$aa_col %||% "")) input$aa_col else NULL
      d <- with_progress_notify(
        wes_oncodrive(rv$maf, aa_col = aa, min_mut = max(2, as.integer(input$min_mut))),
        message = "Scoring positional clustering...")
      res$drv <- d
      res$fdr <- input$fdr
      res$top <- max(5, as.integer(input$top))
      mark_done(rv, "wes_driver")
      log_step(log_rv, "WES drivers",
               params = list(minMut = input$min_mut, fdr = input$fdr,
                             top_interactions = input$top),
               code = c(sprintf('drv <- maftools::oncodrive(maf, minMut = %s, pvalMethod = "zscore")',
                                input$min_mut),
                        sprintf('maftools::plotOncodrive(drv, fdrCutOff = %s)', input$fdr),
                        sprintf('maftools::somaticInteractions(maf, top = %s)', input$top)))
    })

    output$summary <- shiny::renderUI({
      if (is.null(rv$maf)) return(wes_no_maf())
      d <- res$drv
      if (is.null(d)) {
        return(wes_prompt("Click <b>Find drivers</b> to score positional clustering.",
                          "点击<b>检测驱动基因</b>以评估位点聚集程度。"))
      }
      df <- as.data.frame(d)
      sig <- if ("fdr" %in% colnames(df)) sum(df$fdr <= res$fdr, na.rm = TRUE) else NA_integer_
      topg <- if (nrow(df)) as.character(df$Hugo_Symbol[which.min(df$fdr)]) else "-"
      bslib::layout_columns(
        col_widths = c(4, 4, 4),
        stat_tile(i18n("Genes tested", "受检基因数"), format(nrow(df), big.mark = ",")),
        stat_tile(i18n("Below FDR", "低于 FDR"), if (is.na(sig)) "-" else sig),
        stat_tile(i18n("Strongest", "最显著"), topg)
      )
    })

    draw_drv <- function() {
      shiny::req(res$drv)
      n_sig <- tryCatch(sum(as.data.frame(res$drv)$fdr <= res$fdr, na.rm = TRUE),
                        error = function(e) 10)
      maftools::plotOncodrive(res = res$drv, fdrCutOff = res$fdr, useFraction = TRUE,
                              labelSize = adaptive_cex(max(1, n_sig), base = 0.85,
                                                       n_ref = 12, lo = 0.6, hi = 0.95))
    }
    output$drv <- render_base_plot(draw_drv)
    register_figure_download(output, input, "drv", draw_drv, "wes_oncodrive",
                             width = 10, height = 6)

    draw_int <- with_text_boost(function() {
      shiny::req(rv$maf, res$drv)   # gate on the run having happened
      maftools::somaticInteractions(maf = rv$maf, top = res$top,
                                    pvalue = c(0.05, 0.1),
                                    fontSize = adaptive_cex(res$top, n_ref = 25,
                                                            lo = 0.8, hi = 1.2),
                                    countsFontSize = 0.9)
    })
    output$int <- render_base_plot(draw_int)
    register_figure_download(output, input, "int", draw_int, "wes_interactions",
                             width = function() max(8, min(14, 2 + 0.3 * res$top)),
                             height = 8)

    output$tbl_slot <- shiny::renderUI({
      if (is.null(res$drv)) return(wes_no_maf())
      if (has_pkg("DT")) DT::dataTableOutput(ns("tbl"))
      else shiny::verbatimTextOutput(ns("tbl_txt"))
    })
    view <- shiny::reactive({
      d <- res$drv; shiny::req(d)
      df <- as.data.frame(d)
      num <- vapply(df, is.numeric, logical(1))
      df[num] <- lapply(df[num], function(x) signif(x, 3))
      if ("fdr" %in% colnames(df)) df <- df[order(df$fdr), , drop = FALSE]
      df
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
