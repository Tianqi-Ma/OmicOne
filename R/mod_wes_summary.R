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
      en = "This is the sanity check before any analysis. Thousands of variants in
            every sample usually mean unfiltered germline variants or artefacts.
            A few samples far above the rest are often real biology —
            POLE-mutant or mismatch-repair-deficient tumours are hypermutated —
            so look at their signatures before calling them artefacts.",
      zh = "这是任何分析之前的合理性检查。如果每个样本都有上千个变异，通常意味着残留了胚系变异或假阳性。少数样本远高于其余样本则往往是真实生物学——POLE 突变或错配修复缺陷的肿瘤本就是超突变——在判定其为假阳性之前，请先看它们的突变特征。"),
    how  = list(
      en = "<b>Remove outliers</b> keeps one hypermutated sample from flattening
            the boxplot (it only changes the drawing, not the data). Turn the
            <b>dashboard</b> off for a plain stacked barplot of variant
            classifications only.",
      zh = "<b>剔除离群值</b>可避免某个超突变样本把箱线图压平（只影响绘图，不改动数据）。关闭<b>仪表盘</b>则只显示变异分类的堆叠柱状图。"),
    read = list(
      en = "Six panels. Top row: which consequences (<i>Missense</i>,
            <i>Nonsense</i>…), variant types (SNP/INS/DEL) and base changes
            dominate — C>T is the common ageing (clock-like) background; a C>A
            excess has several sources (tobacco, 8-oxoG oxidative damage,
            including during library preparation). Bottom row: variants per
            sample (the dashed line is your chosen statistic), the burden as a
            boxplot, and the top genes with the share of samples hit.",
      zh = "共六个面板。上排：主要的突变后果（<i>错义</i>、<i>无义</i>……）、变异类型（SNP/INS/DEL）与碱基替换——C>T 是常见的衰老（时钟样）背景；C>A 偏多则有多种来源（烟草、8-oxoG 氧化损伤，包括建库过程中产生的）。下排：每样本变异数（虚线为你选择的统计量）、负荷箱线图、以及高频基因及其样本占比。"),
    example = list(
      en = "In TCGA LAML the median sample carries 9 non-synonymous variants and
               <code>FLT3</code> (27%), <code>DNMT3A</code> (25%) and
               <code>NPM1</code> (17%) top the gene list — the profile of a
               low-burden leukaemia.",
      zh = "在 TCGA LAML 中，样本的非同义变异中位数为 9 个，基因列表由 <code>FLT3</code>（27%）、<code>DNMT3A</code>（25%）、<code>NPM1</code>（17%）领衔——这是低突变负荷白血病的典型样貌。")
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
    subtitle  = list(en = "The cohort at a glance: mutation types, per-sample burden, top genes.",
                     zh = "队列全景一瞥：突变类型、每样本负荷、高频基因。"),
    explainer = explainer,
    controls  = controls,
    summary   = shiny::uiOutput(ns("summary")),
    preview   = shiny::tagList(
      shiny::uiOutput(ns("insight")),
      bslib::navset_card_tab(
        bslib::nav_panel(i18n("Dashboard", "仪表盘"),
                         preview_plot_ui(ns("plot"), download = TRUE,
                                         guide = list(en = "A six-panel dashboard of the cohort's mutation landscape will be drawn here.",
                                                      zh = "运行后，这里将绘制队列突变全景的六联仪表盘。"),
                                         caption = list(en = "Consequences, variant types, base changes, non-synonymous variants per sample, and the most mutated genes (share of samples).",
                                                        zh = "突变后果、变异类型、碱基替换、每样本非同义变异数与高频基因（样本占比）。"))),
        bslib::nav_panel(i18n("Gene frequencies", "基因频率"), shiny::uiOutput(ns("tbl_slot")))
      ))
  )
}

#' @rdname mod_wes_summary
#' @keywords internal
mod_wes_summary_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    ns   <- session$ns
    opts <- step_result(rv, "wes")

    draw_with <- function(maf, o) {
      n_samples <- tryCatch(nrow(maftools::getSampleSummary(maf)),
                            error = function(e) NA_integer_)
      maftools::plotmafSummary(maf = maf, rmOutlier = o$rm_outlier,
                               addStat = if (identical(o$stat, "none")) NULL else o$stat,
                               dashboard = o$dashboard, titvRaw = FALSE,
                               top = o$top,
                               fs = adaptive_cex(n_samples, base = 1.2, n_ref = 200,
                                                 lo = 1.0, hi = 1.25),
                               textSize = 1.1, titleSize = c(1.25, 1.05))
    }

    shiny::observeEvent(input$run, {
      shiny::req(rv$maf)
      if (!require_pkgs("maftools", "Cohort summary")) return(NULL)
      o <- list(top = wes_int(input$top, 10, lo = 3, hi = 30),
                stat = input$stat %||% "median",
                rm_outlier = isTRUE(input$rm_outlier),
                dashboard = isTRUE(input$dashboard))
      maf <- rv$maf
      ok <- with_progress_notify(wes_dry_run(function() draw_with(maf, o)),
                                 message = "Drawing the cohort summary...")
      if (is.null(ok)) return(NULL)
      opts(o)
      mark_done(rv, "wes_summary")
      log_step(log_rv, "WES cohort summary",
               params = list(top = o$top, addStat = o$stat, rmOutlier = o$rm_outlier,
                             dashboard = o$dashboard),
               code = wes_code("maftools::plotmafSummary",
                               list(maf = quote(maf), rmOutlier = o$rm_outlier,
                                    addStat = if (identical(o$stat, "none")) NULL else o$stat,
                                    dashboard = o$dashboard, titvRaw = FALSE, top = o$top)))
    })

    stats <- shiny::reactive({
      shiny::req(rv$maf)
      wes_overview(rv$maf)
    })

    output$summary <- shiny::renderUI({
      if (is.null(rv$maf)) return(wes_no_maf())
      if (is.null(opts())) {
        return(wes_prompt("Set the options and click <b>Draw summary</b>.",
                          "设置选项后点击<b>绘制概览</b>。"))
      }
      ov <- stats()
      shiny::tagList(
        stat_tile(i18n("Non-synonymous variants", "非同义变异数"), wes_fmt(ov$variants)),
        stat_tile(i18n("Median / sample", "中位数/样本"), wes_fmt(ov$median_per_sample)),
        stat_tile(i18n("Most mutated", "最高频基因"),
                  if (is.na(ov$top_gene)) "-"
                  else sprintf("%s (%.0f%%)", ov$top_gene, ov$top_pct))
      )
    })

    output$insight <- shiny::renderUI({
      if (is.null(rv$maf) || is.null(opts())) return(NULL)
      ov <- stats()
      top <- if (is.na(ov$top_gene)) "-" else sprintf("%s (%.0f%%)", ov$top_gene, ov$top_pct)
      insight_bar(
        sprintf("The cohort carries <b>%s</b> non-synonymous variants — median <b>%s</b> per sample, led by <b>%s</b>. Samples far above the median can be real hypermutators (POLE, mismatch-repair deficiency) or filtering problems; check their signatures before excluding them.",
                wes_fmt(ov$variants), wes_fmt(ov$median_per_sample), top),
        sprintf("该队列共 <b>%s</b> 个非同义变异，每样本中位 <b>%s</b> 个，最高频基因为 <b>%s</b>。远高于中位数的样本可能是真实的超突变（POLE、错配修复缺陷），也可能是过滤问题；排除前请先查看其突变特征。",
                wes_fmt(ov$variants), wes_fmt(ov$median_per_sample), top))
    })

    draw_summary <- function() {
      o <- opts()
      shiny::req(rv$maf, o)
      draw_with(rv$maf, o)
    }
    output$plot <- render_base_plot(draw_summary)
    register_figure_download(output, input, "plot", draw_summary,
                             "wes_cohort_summary", width = 12, height = 9)

    gene_df <- shiny::reactive({
      shiny::req(rv$maf)
      wes_gene_table(rv$maf, n = Inf)
    })
    tb <- wes_table(ns, "tbl", gene_df, "wes_gene_frequencies",
                    has_maf = function() !is.null(rv$maf),
                    ready = function() !is.null(opts()))
    output$tbl_slot <- tb$slot
    output$tbl <- tb$table
    output$tbl_dl <- tb$download
  })
}
