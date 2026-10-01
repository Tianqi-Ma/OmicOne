#' WES module 2: Cohort summary
#'
#' maftools' dashboard view of the whole cohort: variant classifications, variant
#' types, SNV classes, per-sample burden, and the top mutated genes. The first
#' look at a MAF, before any hypothesis.
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_wes_summary
NULL

#' @rdname mod_wes_summary
#' @keywords internal
mod_wes_summary_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = list(en = "Cohort summary", zh = "队列概览"),
    what = list(
      en = "A dashboard of the whole cohort: what kinds of mutations there are,
            how many each sample carries, and which genes are hit most often.",
      zh = "整个队列的仪表盘：都有哪些类型的突变、每个样本携带多少、哪些基因被击中得最频繁。"),
    why  = list(
      en = "This is the sanity check before any analysis. A cohort where every
            sample has thousands of variants, or where one sample has ten times
            the rest, usually means a calling or filtering problem — not biology.",
      zh = "这是任何分析之前的合理性检查。如果每个样本都有上千个变异，或某个样本比其他样本多十倍，通常意味着变异检出或过滤有问题，而不是生物学差异。"),
    how  = list(
      en = "<b>Remove outliers</b> keeps one hypermutated sample from flattening
            the boxplot. Turn the <b>dashboard</b> off for a plain stacked
            barplot of variant classifications only.",
      zh = "<b>剔除离群值</b>可避免某个超突变样本把箱线图压平。关闭<b>仪表盘</b>则只显示变异分类的堆叠柱状图。"),
    example = list(
      en = "In TCGA LAML most samples carry only a handful of variants and
               <code>FLT3</code>, <code>DNMT3A</code>, <code>NPM1</code> top the
               gene list — that is what a clean cohort looks like.",
      zh = "在 TCGA LAML 中，多数样本只携带少数变异，基因列表由 <code>FLT3</code>、<code>DNMT3A</code>、<code>NPM1</code> 领衔——干净的队列就该是这个样子。")
  )
  controls <- shiny::tagList(
    label_with_help("Top genes", "How many genes to show in the summary panel.",
                    label_zh = "显示基因数", tip_zh = "概览面板中显示多少个基因。"),
    shiny::numericInput(ns("top"), NULL, value = 10, min = 3, max = 30, step = 1),
    label_with_help("Per-sample statistic", "The line drawn over the burden boxplot.",
                    label_zh = "每样本统计量", tip_zh = "叠加在突变负荷箱线图上的统计线。"),
    shiny::selectInput(ns("stat"), NULL,
                       c("Median" = "median", "Mean" = "mean", "None" = "none"),
                       selected = "median"),
    shiny::checkboxInput(ns("rm_outlier"),
                         i18n("Remove outliers from the boxplot", "从箱线图中剔除离群值"),
                         value = TRUE),
    shiny::checkboxInput(ns("dashboard"),
                         i18n("Full dashboard", "完整仪表盘"), value = TRUE),
    run_button(ns("run"), "Draw summary", "绘制概览")
  )
  step_container(
    title     = list(en = "Cohort summary", zh = "队列概览"),
    explainer = explainer,
    controls  = controls,
    summary   = shiny::uiOutput(ns("summary")),
    preview   = bslib::navset_card_tab(
      bslib::nav_panel(i18n("Dashboard", "仪表盘"), preview_plot_ui(ns("plot"), download = TRUE)),
      bslib::nav_panel(i18n("Gene frequencies", "基因频率"), shiny::uiOutput(ns("tbl_slot")))
    )
  )
}

#' @rdname mod_wes_summary
#' @keywords internal
mod_wes_summary_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    ns   <- session$ns
    opts <- shiny::reactiveVal(NULL)

    shiny::observeEvent(input$run, {
      shiny::req(rv$maf)
      if (!require_pkgs("maftools", "Cohort summary")) return(NULL)
      opts(list(top = max(3, as.integer(input$top)),
                stat = input$stat, rm_outlier = isTRUE(input$rm_outlier),
                dashboard = isTRUE(input$dashboard)))
      mark_done(rv, "wes_summary")
      log_step(log_rv, "WES cohort summary",
               params = list(top = input$top, addStat = input$stat,
                             rmOutlier = input$rm_outlier),
               code = sprintf('maftools::plotmafSummary(maf, top = %s, addStat = "%s", rmOutlier = %s, dashboard = %s)',
                              input$top, input$stat, input$rm_outlier, input$dashboard))
    })

    output$summary <- shiny::renderUI({
      if (is.null(rv$maf)) return(wes_no_maf())
      if (is.null(opts())) {
        return(wes_prompt("Set the options and click <b>Draw summary</b>.",
                          "设置选项后点击<b>绘制概览</b>。"))
      }
      ov <- wes_overview(rv$maf)
      fmt <- function(x) if (is.null(x) || is.na(x)) "-" else format(round(x), big.mark = ",")
      bslib::layout_columns(
        col_widths = c(4, 4, 4),
        stat_tile(i18n("Variants", "变异总数"), fmt(ov$variants)),
        stat_tile(i18n("Median / sample", "中位数/样本"), fmt(ov$median_per_sample)),
        stat_tile(i18n("Most mutated", "最高频基因"),
                  if (is.na(ov$top_gene)) "-"
                  else sprintf("%s (%.0f%%)", ov$top_gene, ov$top_pct))
      )
    })

    draw_summary <- function() {
      o <- opts(); shiny::req(rv$maf, o)
      n_samples <- tryCatch(nrow(maftools::getSampleSummary(rv$maf)),
                            error = function(e) NA_integer_)
      maftools::plotmafSummary(maf = rv$maf, rmOutlier = o$rm_outlier,
                               addStat = if (identical(o$stat, "none")) NULL else o$stat,
                               dashboard = o$dashboard, titvRaw = FALSE,
                               top = o$top,
                               fs = adaptive_cex(n_samples, base = 1.2, n_ref = 200,
                                                 lo = 1.0, hi = 1.25),
                               textSize = 1.1, titleSize = c(1.25, 1.05))
    }
    output$plot <- render_base_plot(draw_summary)
    register_figure_download(output, input, "plot", draw_summary,
                             "wes_cohort_summary", width = 12, height = 9)

    output$tbl_slot <- shiny::renderUI({
      if (is.null(rv$maf)) return(wes_no_maf())
      if (has_pkg("DT")) DT::dataTableOutput(ns("tbl"))
      else shiny::verbatimTextOutput(ns("tbl_txt"))
    })
    gene_df <- shiny::reactive({ shiny::req(rv$maf); wes_gene_table(rv$maf, n = 200) })
    if (has_pkg("DT")) {
      output$tbl <- DT::renderDataTable(
        DT::datatable(gene_df(), rownames = FALSE, filter = "top",
                      options = list(pageLength = 15, scrollX = TRUE)))
    } else {
      output$tbl_txt <- shiny::renderPrint(utils::head(gene_df(), 20))
    }
  })
}
