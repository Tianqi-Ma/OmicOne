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
            exclusivity (Fisher's exact test, BH-adjusted).",
      zh = "<b>Oncodrive</b> 找出突变集中在少数位点、而非均匀散布的基因。<b>互作分析</b>对每一对高频突变基因检验共现或互斥（Fisher 精确检验，BH 校正）。"),
    why  = list(
      en = "A gene can be mutated often just because it is long. Positional
            clustering is evidence of selection, not size — but only of the
            oncogene kind: tumour suppressors hit by truncating mutations
            scattered along the gene do not cluster and are missed, so use a
            frequency-based method (dNdScv, MutSigCV) for those. Mutual
            exclusivity suggests two genes hit the same pathway.",
      zh = "一个基因突变频繁，可能只是因为它长。位点聚集是受选择的证据，而不是长度的结果——但只针对癌基因类：截短突变散布全长的抑癌基因不会聚集，会被漏掉，这类基因请用基于频率的方法（dNdScv、MutSigCV）。互斥则提示两个基因作用于同一通路。"),
    how  = list(
      en = "Raise <b>minimum mutations</b> on a large cohort to cut noise.
            Interactions are confounded by burden: hypermutated samples carry
            many genes at once and inflate co-occurrence, so check the pairs
            against TMB before reading them as biology.",
      zh = "在大队列中提高<b>最小突变数</b>可以降噪。互作结果会受突变负荷混杂：超突变样本同时携带许多基因的突变，会虚增共现，解读为生物学之前请先对照 TMB 检查这些基因对。"),
    read = list(
      en = "<b>Oncodrive</b>: one bubble per tested gene — x is the fraction of its
            mutations inside clusters, y is −log10(FDR), bubble size the number
            of clusters; genes below your FDR cutoff are red and labelled.
            <b>Interactions</b>: every tile is a gene pair coloured by
            −log10(FDR) — green co-occurs, brown is exclusive; * marks FDR < 0.05
            and · FDR < 0.1. The <b>Interaction table</b> lists every pair with
            its counts.",
      zh = "<b>Oncodrive</b>：每个气泡是一个受检基因——横轴为其突变落在簇内的比例，纵轴为 −log10(FDR)，气泡大小为簇数；低于 FDR 阈值的基因标红并标注。<b>互作</b>：每个格子是一对基因，按 −log10(FDR) 着色——绿色为共现，棕色为互斥；* 表示 FDR < 0.05，· 表示 FDR < 0.1。<b>互作结果表</b>列出每一对基因及其计数。"),
    example = list(
      en = "In TCGA LAML, <code>IDH1</code>, <code>IDH2</code> and
               <code>NPM1</code> are the top positional-clustering hits; among
               the interactions, <code>FLT3</code>–<code>NPM1</code> and
               <code>DNMT3A</code>–<code>NPM1</code> co-occur.",
      zh = "在 TCGA LAML 中，<code>IDH1</code>、<code>IDH2</code>、<code>NPM1</code> 是位点聚集最显著的基因；互作结果中，<code>FLT3</code>–<code>NPM1</code> 与 <code>DNMT3A</code>–<code>NPM1</code> 共现。")
  )
  controls <- shiny::tagList(
    label_with_help("Minimum mutations per gene",
                    "Genes with fewer mutations than this are not tested by oncodrive.",
                    label_zh = "每基因最少突变数",
                    tip_zh = "突变数少于此值的基因不参与 oncodrive 检验。"),
    shiny::numericInput(ns("min_mut"), NULL, value = 5, min = 2, max = 50, step = 1),
    label_with_help("FDR cutoff", "Drivers with FDR below this are highlighted.",
                    label_zh = "FDR 阈值", tip_zh = "FDR 低于此值的驱动基因会被高亮。"),
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
    subtitle  = list(en = "Positional-clustering driver calls and pairwise gene interactions.",
                     zh = "基于位点聚集识别驱动基因，并检验基因两两互作。"),
    explainer = explainer,
    controls  = controls,
    summary   = shiny::uiOutput(ns("summary")),
    preview   = shiny::tagList(
      shiny::uiOutput(ns("insight")),
      bslib::navset_card_tab(
        bslib::nav_panel(i18n("Oncodrive", "Oncodrive"),
                         preview_plot_ui(ns("drv"), download = TRUE,
                                         guide = list(en = "Every tested gene's clustering score and FDR will be drawn here.",
                                                      zh = "运行后，这里将绘制每个受检基因的聚集程度与 FDR。"),
                                         caption = list(en = "Bubble = tested gene; x = fraction of its mutations in clusters; y = −log10(FDR); size = number of clusters. All tested genes are drawn; red = FDR below the cutoff.",
                                                        zh = "气泡＝受检基因；横轴＝突变落在簇内的比例；纵轴＝−log10(FDR)；大小＝簇数。所有受检基因均绘出；红色＝FDR 低于阈值。"))),
        bslib::nav_panel(i18n("Interactions", "基因互作"),
                         preview_plot_ui(ns("int"), download = TRUE,
                                         guide = list(en = "Pairwise co-occurrence / exclusivity of the top genes will be drawn here.",
                                                      zh = "运行后，这里将绘制高频基因两两之间的共现/互斥。"),
                                         caption = list(en = "Tile = gene pair, colour = −log10(BH FDR), green co-occurring / brown exclusive; * FDR < 0.05, · FDR < 0.1. Burden confounds co-occurrence.",
                                                        zh = "格子＝基因对，颜色＝−log10(BH FDR)，绿色共现/棕色互斥；* FDR < 0.05，· FDR < 0.1。突变负荷会混杂共现。"))),
        bslib::nav_panel(i18n("Driver table", "驱动基因表"), shiny::uiOutput(ns("tbl_slot"))),
        bslib::nav_panel(i18n("Interaction table", "互作结果表"), shiny::uiOutput(ns("int_slot")))
      ))
  )
}

#' @rdname mod_wes_driver
#' @keywords internal
mod_wes_driver_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    ns  <- session$ns
    res <- step_results(rv, "wes", drv = NULL, bg_note = character(0), int = NULL,
                        fdr = 0.1, top = 25)

    output$aa_ui <- shiny::renderUI({
      shiny::req(rv$maf)
      wes_col_select(ns, "aa_col",
                     label = list(en = "Protein change column", zh = "蛋白改变列"),
                     tip = list(en = "Oncodrive needs amino-acid positions. Auto-detected.",
                                zh = "Oncodrive 需要氨基酸位置信息，会自动识别。"),
                     choices = wes_fields(rv$maf), selected = wes_guess_aa_col(rv$maf) %||% "",
                     none = "(auto)")
    })

    shiny::observeEvent(input$run, {
      shiny::req(rv$maf)
      if (!require_pkgs("maftools", "Driver detection")) return(NULL)
      aa <- if (nzchar(input$aa_col %||% "")) input$aa_col else NULL
      min_mut <- wes_int(input$min_mut, 5, lo = 2, hi = 50)
      fdr <- wes_prob(input$fdr, 0.1)
      top <- wes_int(input$top, 25, lo = 5, hi = 60)
      maf <- rv$maf
      d <- with_progress_notify(wes_oncodrive(maf, aa_col = aa, min_mut = min_mut),
                                message = "Scoring positional clustering...")
      if (is.null(d)) return(NULL)
      it <- with_progress_notify(wes_interactions(maf, top = top),
                                 message = "Testing gene pairs...")
      if (is.null(it)) return(NULL)
      res$drv <- d$res
      res$bg_note <- d$bg_note
      res$int <- it
      res$fdr <- fdr
      res$top <- top
      mark_done(rv, "wes_driver")
      log_step(log_rv, "WES drivers",
               params = list(minMut = min_mut, fdr = fdr, AACol = aa %||% "(auto)",
                             top_interactions = top),
               code = c(wes_code("maftools::oncodrive",
                                 list(maf = quote(maf), AACol = aa, minMut = min_mut,
                                      pvalMethod = "zscore"), assign = "drv"),
                        wes_code("maftools::plotOncodrive",
                                 list(res = quote(drv), fdrCutOff = fdr, useFraction = TRUE)),
                        wes_code("maftools::somaticInteractions",
                                 list(maf = quote(maf), top = top, pvalue = c(0.05, 0.1),
                                      plotPadj = if (wes_has_arg("somaticInteractions",
                                                                 "plotPadj")) TRUE),
                                 assign = "pairs")))
    })

    stats <- shiny::reactive({
      shiny::req(res$drv)
      df <- res$drv
      q <- if ("fdr" %in% names(df)) df$fdr else rep(NA_real_, nrow(df))
      it <- res$int
      sig_pairs <- if (!is.null(it) && "pAdj" %in% names(it)) sum(it$pAdj < 0.05, na.rm = TRUE)
                   else NA_integer_
      list(tested = nrow(df), sig = sum(q < res$fdr, na.rm = TRUE),
           top = if (nrow(df) && any(!is.na(q))) as.character(df$Hugo_Symbol[which.min(q)])
                 else "-",
           pairs = if (!is.null(it)) nrow(it) else NA_integer_, sig_pairs = sig_pairs)
    })

    output$summary <- shiny::renderUI({
      if (is.null(rv$maf)) return(wes_no_maf())
      if (is.null(res$drv)) {
        return(wes_prompt("Click <b>Find drivers</b> to score positional clustering.",
                          "点击<b>检测驱动基因</b>以评估位点聚集程度。"))
      }
      s <- stats()
      shiny::tagList(
        stat_tile(i18n("Genes tested", "受检基因数"), wes_fmt(s$tested)),
        stat_tile(sprintf("FDR < %g", res$fdr), s$sig),
        stat_tile(i18n("Strongest", "最显著"), s$top),
        stat_tile(i18n("Pairs at FDR < 0.05", "FDR < 0.05 的基因对"),
                  sprintf("%s / %s", wes_fmt(s$sig_pairs), wes_fmt(s$pairs)))
      )
    })

    output$insight <- shiny::renderUI({
      if (is.null(res$drv)) return(NULL)
      s <- stats()
      bg_en <- ""
      bg_zh <- ""
      if (length(res$bg_note)) {
        bg_en <- " Too few synonymous variants to estimate a background: maftools used its preset one (mean 0.279, SD 0.13), so the p-values rest on that assumption."
        bg_zh <- "同义变异太少，无法估计背景分布：maftools 使用了预设背景（均值 0.279，标准差 0.13），因此 p 值依赖这一假设。"
      }
      if (s$sig > 0) {
        insight_bar(
          sprintf("<b>%d</b> of %s tested genes pass FDR < %g — strongest: <b>%s</b>. Positional clustering finds oncogene-type hotspots; scattered truncating tumour suppressors need a frequency-based test (dNdScv).%s",
                  s$sig, wes_fmt(s$tested), res$fdr, s$top, bg_en),
          sprintf("%s 个受检基因中有 <b>%d</b> 个通过 FDR < %g——最显著为 <b>%s</b>。位点聚集法找的是癌基因型热点；截短突变散布的抑癌基因需要基于频率的检验（dNdScv）。%s",
                  wes_fmt(s$tested), s$sig, res$fdr, s$top, bg_zh))
      } else {
        insight_bar(
          sprintf("No gene passes FDR < %g among the %s tested. That is common on small cohorts; consider a frequency-based method such as dNdScv rather than loosening the cutoff.%s",
                  res$fdr, wes_fmt(s$tested), bg_en),
          sprintf("受检的 %s 个基因无一通过 FDR < %g。小队列中这很常见；可考虑 dNdScv 等基于频率的方法，而不是放宽阈值。%s",
                  wes_fmt(s$tested), res$fdr, bg_zh))
      }
    })

    draw_drv <- function() {
      shiny::req(res$drv)
      s <- stats()
      maftools::plotOncodrive(res = data.table::as.data.table(res$drv), fdrCutOff = res$fdr,
                              useFraction = TRUE,
                              labelSize = adaptive_cex(max(1, s$sig), base = 0.85,
                                                       n_ref = 12, lo = 0.6, hi = 0.95))
    }
    output$drv <- render_base_plot(draw_drv)
    register_figure_download(output, input, "drv", draw_drv, "wes_oncodrive",
                             width = 10, height = 6)

    draw_int <- with_text_boost(function() {
      shiny::req(rv$maf, res$int)
      args <- list(maf = rv$maf, top = res$top, pvalue = c(0.05, 0.1),
                   fontSize = adaptive_cex(res$top, n_ref = 25, lo = 0.8, hi = 1.2),
                   countsFontSize = 0.9)
      if (wes_has_arg("somaticInteractions", "plotPadj")) args$plotPadj <- TRUE
      invisible(do.call(maftools::somaticInteractions, args))
    })
    output$int <- render_base_plot(draw_int)
    register_figure_download(output, input, "int", draw_int, "wes_interactions",
                             width = function() max(8, min(14, 2 + 0.3 * res$top)),
                             height = 8)

    has_maf <- function() !is.null(rv$maf)
    drv_view <- shiny::reactive({
      shiny::req(res$drv)
      df <- res$drv
      num <- vapply(df, is.numeric, logical(1))
      df[num] <- lapply(df[num], function(x) signif(x, 3))
      if ("fdr" %in% colnames(df)) df <- df[order(df$fdr), , drop = FALSE]
      df
    })
    tb <- wes_table(ns, "tbl", drv_view, "wes_oncodrive", has_maf,
                    function() !is.null(res$drv))
    output$tbl_slot <- tb$slot
    output$tbl <- tb$table
    output$tbl_dl <- tb$download

    int_view <- shiny::reactive({
      shiny::req(res$int)
      df <- res$int
      keep <- intersect(c("gene1", "gene2", "Event", "oddsRatio", "pValue", "pAdj", "11",
                          "10", "01", "00", "event_ratio"), names(df))
      df <- df[, keep, drop = FALSE]
      num <- vapply(df, is.numeric, logical(1))
      df[num] <- lapply(df[num], function(x) signif(x, 3))
      if ("pAdj" %in% names(df)) df <- df[order(df$pAdj), , drop = FALSE]
      df
    })
    ti <- wes_table(ns, "int_tbl", int_view, "wes_interactions", has_maf,
                    function() !is.null(res$int))
    output$int_slot <- ti$slot
    output$int_tbl <- ti$table
    output$int_tbl_dl <- ti$download
  })
}
