#' WES module 10: Cohort comparison
#'
#' Split the cohort by a clinical column and test, gene by gene, whether the two
#' groups are mutated at different rates — the mutation-level equivalent of a
#' differential expression test. Significance is the BH q-value over all genes
#' tested (mafCompare's `adjPval`).
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_wes_compare
NULL

#' @rdname mod_wes_compare
#' @keywords internal
mod_wes_compare_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = list(en = "Cohort comparison", zh = "队列比较"),
    what = list(
      en = "Pick a clinical column and two of its levels. Every gene mutated in
            enough samples of either group is tested with Fisher's exact test
            for a difference in the share of mutated samples, and the p-values
            are BH-adjusted across all genes tested.",
      zh = "选择一个临床列及其中两个水平。对在任一组中突变样本数足够多的每个基因，用 Fisher 精确检验比较两组的突变样本比例，并对所有受检基因的 p 值做 BH 校正。"),
    why  = list(
      en = "This is how you turn a cohort into a comparison: responders versus
            non-responders, primary versus metastatic, treated versus naive.
            The forest plot shows the odds ratio and its confidence interval per
            gene, so you can see effect size and not just a q-value.",
      zh = "这就是把一个队列变成一次比较的方法：应答 vs 不应答、原发 vs 转移、治疗过 vs 初治。森林图给出每个基因的比值比及其置信区间，让你看到效应量而不只是 q 值。"),
    how  = list(
      en = "<b>Minimum mutated samples</b> is applied within each group: a gene is
            tested when at least that many samples of either group carry it.
            Raising it avoids testing genes hit in one or two patients, which
            add tests (and so lower power after BH) without any chance of
            reaching significance.",
      zh = "<b>最少突变样本数</b>按组内计算：任一组中至少有这么多样本携带该基因突变时，才对它做检验。提高该值可避免检验只在一两个患者中突变的基因——它们几乎不可能显著，却会增加检验次数（BH 校正后降低检验效能）。"),
    read = list(
      en = "The forest plot gives one row per gene passing your FDR cutoff, with
            the raw counts in the label: the dot is the odds ratio (log axis),
            the whiskers its 95% CI. Right of 1 means mutated more often in the
            first group. An open dot means a zero cell (0 mutated, or all
            mutated, in a group): its odds ratio is drawn with a +1
            pseudo-count, while the q-value comes from the uncorrected test.
            The <b>Frequencies</b> tab shows the same genes as side-by-side bars.",
      zh = "森林图每行一个通过 FDR 阈值的基因，标签中给出原始计数：圆点是比值比（对数轴），横须是 95% 置信区间。1 以右表示在第一组中突变更多。空心圆点表示存在零格（某组中无人突变或全部突变）：其比值比按 +1 伪计数绘制，q 值则来自未校正的检验。<b>频率对比</b>页签以并排柱条展示同一批基因。"),
    example = list(
      en = "In TCGA LAML, M2 (44 samples) versus M3 (21) tests 7 genes:
               <code>NPM1</code> 8 vs 0 mutated gives raw p = 0.046 but q = 0.20 —
               not significant once the multiple tests are accounted for.",
      zh = "在 TCGA LAML 中比较 M2（44 个样本）与 M3（21 个），共检验 7 个基因：<code>NPM1</code> 突变 8 vs 0，原始 p = 0.046，但 q = 0.20——考虑多重检验后并不显著。")
  )
  controls <- shiny::tagList(
    shiny::uiOutput(ns("feat_ui")),
    shiny::uiOutput(ns("level_ui")),
    label_with_help("Minimum mutated samples (in either group)",
                    "A gene is tested if at least this many samples of either group carry it.",
                    label_zh = "最少突变样本数（任一组内）",
                    tip_zh = "任一组中至少有这么多样本携带该基因突变时，才对其做检验。"),
    shiny::numericInput(ns("min_mut"), NULL, value = 5, min = 1, max = 50, step = 1),
    label_with_help("FDR cutoff", "Genes with a BH q-value below this are drawn and counted.",
                    label_zh = "FDR 阈值", tip_zh = "BH q 值低于此值的基因会被绘制并计数。"),
    shiny::numericInput(ns("fdr"), NULL, value = 0.1, min = 0.001, max = 0.25, step = 0.01),
    run_button(ns("run"), "Compare cohorts", "比较队列")
  )
  step_container(id = id, 
    title     = list(en = "Cohort comparison", zh = "队列比较"),
    subtitle  = list(en = "Fisher tests (BH-adjusted) for genes that differ between two groups.",
                     zh = "用 Fisher 检验（BH 校正）找出两组间突变比例不同的基因。"),
    explainer = explainer,
    controls  = controls,
    summary   = shiny::uiOutput(ns("summary")),
    preview   = shiny::tagList(
      shiny::uiOutput(ns("insight")),
      bslib::navset_card_tab(
        bslib::nav_panel(i18n("Forest plot", "森林图"),
                         preview_plot_ui(ns("forest"), download = TRUE,
                                         guide = list(en = "Odds ratios of the genes passing the FDR cutoff will be drawn here.",
                                                      zh = "运行后，这里将绘制通过 FDR 阈值的基因的比值比。"),
                                         caption = list(en = "Dot = odds ratio (log axis), whiskers = 95% CI; right of 1 = mutated more often in group 1; open dot = zero cell, drawn with a +1 pseudo-count. Rows: genes at FDR below the cutoff.",
                                                        zh = "圆点＝比值比（对数轴），横须＝95% CI；1 以右＝在第 1 组中突变更多；空心点＝存在零格，按 +1 伪计数绘制。行：FDR 低于阈值的基因。"))),
        bslib::nav_panel(i18n("Frequencies", "频率对比"),
                         preview_plot_ui(ns("cobar"), download = TRUE,
                                         guide = list(en = "Side-by-side mutation frequencies of the significant genes will be drawn here.",
                                                      zh = "运行后，这里将并排绘制显著基因在两组中的突变频率。"),
                                         caption = list(en = "Share of samples mutated in each group, for the genes passing the FDR cutoff only.",
                                                        zh = "各组中携带突变的样本比例，仅含通过 FDR 阈值的基因。"))),
        bslib::nav_panel(i18n("Results", "结果表"), shiny::uiOutput(ns("tbl_slot")))
      ))
  )
}

#' @rdname mod_wes_compare
#' @keywords internal
mod_wes_compare_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    ns  <- session$ns
    res <- step_results(rv, "wes", cmp = NULL, fdr = 0.1, l1 = NULL, l2 = NULL,
                        feature = NULL)

    clin_cd <- shiny::reactive({
      shiny::req(rv$maf)
      tryCatch(as.data.frame(maftools::getClinicalData(rv$maf)), error = function(e) NULL)
    })

    output$feat_ui <- shiny::renderUI({
      shiny::req(rv$maf)
      cd <- clin_cd()
      cols <- if (is.null(cd)) character(0)
              else setdiff(categorical_cols(cd, min_levels = 2, max_levels = 50, allow_na = TRUE),
                           "Tumor_Sample_Barcode")
      if (!length(cols)) {
        return(shiny::div(class = "omicone-status-empty",
                          i18n("This MAF has no categorical clinical column. Load a clinical table on the Import step.",
                               "该 MAF 没有分类型临床列。请在导入步骤加载临床表。")))
      }
      wes_col_select(ns, "feature",
                     label = list(en = "Split by", zh = "分组依据"),
                     tip = list(en = "The clinical column defining the two groups.",
                                zh = "用于定义两个分组的临床列。"),
                     choices = cols, selected = keep_selected(shiny::isolate(input$feature), cols))
    })

    feature_levels <- shiny::reactive({
      cd <- clin_cd()
      shiny::req(cd, input$feature, input$feature %in% names(cd))
      wes_feature_levels(cd, input$feature, min_n = 2)
    })

    output$level_ui <- shiny::renderUI({
      lv <- tryCatch(feature_levels(), error = function(e) character(0))
      if (length(lv) < 2) {
        return(shiny::div(class = "omicone-status-empty",
                          i18n("That column needs at least two levels with 2+ samples each.",
                               "该列需要至少两个水平，且每个水平至少有 2 个样本。")))
      }
      shiny::tagList(
        shiny::selectInput(ns("l1"), i18n("Group 1", "组 1"), lv, selected = lv[1]),
        shiny::selectInput(ns("l2"), i18n("Group 2", "组 2"), lv, selected = lv[2])
      )
    })

    shiny::observeEvent(input$run, {
      shiny::req(rv$maf, input$feature, input$l1, input$l2)
      if (!require_pkgs("maftools", "Cohort comparison")) return(NULL)
      feature <- input$feature
      l1 <- input$l1
      l2 <- input$l2
      if (identical(l1, l2)) {
        wes_notify("Pick two different groups.", "请选择两个不同的分组。", duration = 8)
        return(NULL)
      }
      min_mut <- wes_int(input$min_mut, 5, lo = 1, hi = 50)
      fdr <- wes_prob(input$fdr, 0.1)
      out <- with_progress_notify(
        wes_compare_cohorts(rv$maf, feature, l1, l2, min_mut = min_mut),
        message = "Comparing cohorts...")
      if (is.null(out)) return(NULL)
      res$cmp  <- out
      res$fdr  <- fdr
      res$l1   <- l1
      res$l2   <- l2
      res$feature <- feature
      mark_done(rv, "wes_compare")
      q <- function(x) deparse(x)
      log_step(log_rv, "WES cohort comparison",
               params = list(feature = feature, group1 = l1, group2 = l2, minMut = min_mut,
                             fdr = fdr, n1 = out$n1, n2 = out$n2),
               code = c("cd <- as.data.frame(maftools::getClinicalData(maf))",
                        sprintf("s1 <- as.character(cd$Tumor_Sample_Barcode[which(as.character(cd[[%s]]) == %s)])",
                                q(feature), q(l1)),
                        sprintf("s2 <- as.character(cd$Tumor_Sample_Barcode[which(as.character(cd[[%s]]) == %s)])",
                                q(feature), q(l2)),
                        "m1 <- maftools::subsetMaf(maf = maf, tsb = s1)",
                        "m2 <- maftools::subsetMaf(maf = maf, tsb = s2)",
                        wes_code("maftools::mafCompare",
                                 list(m1 = quote(m1), m2 = quote(m2), m1Name = l1, m2Name = l2,
                                      minMut = min_mut), assign = "cmp"),
                        "cmp$SampleSummary   # samples each test used",
                        sprintf("sig <- cmp$results[adjPval < %s]   # BH q-values", fdr),
                        wes_code("maftools::forestPlot",
                                 list(mafCompareRes = quote(cmp), fdr = fdr)),
                        wes_code("maftools::coBarplot",
                                 list(m1 = quote(m1), m2 = quote(m2), m1Name = l1, m2Name = l2,
                                      genes = quote(sig$Hugo_Symbol)))))
    })

    stats <- shiny::reactive({
      c0 <- res$cmp
      shiny::req(c0)
      d <- as.data.frame(c0$res$results)
      fd <- wes_forest_data(c0$res, res$fdr, c0$n1, c0$n2)
      list(tested = nrow(d), sig = nrow(fd), forest = fd, n1 = c0$n1, n2 = c0$n2,
           genes = as.character(fd$Hugo_Symbol),
           min_p = if (nrow(d)) min(d$pval, na.rm = TRUE) else NA_real_)
    })

    output$summary <- shiny::renderUI({
      if (is.null(rv$maf)) return(wes_no_maf())
      if (is.null(res$cmp)) {
        return(wes_prompt("Choose two groups and click <b>Compare cohorts</b>.",
                          "选择两个分组后点击<b>比较队列</b>。"))
      }
      s <- stats()
      shiny::tagList(
        stat_tile(sprintf("%s (n)", res$l1), wes_fmt(s$n1)),
        stat_tile(sprintf("%s (n)", res$l2), wes_fmt(s$n2)),
        stat_tile(i18n("Genes tested", "受检基因数"), wes_fmt(s$tested)),
        stat_tile(sprintf("FDR < %g", res$fdr), s$sig)
      )
    })

    output$insight <- shiny::renderUI({
      if (is.null(res$cmp)) return(NULL)
      s <- stats()
      small_en <- if (min(s$n1, s$n2) < 10) " One group has fewer than 10 samples, so only very large differences can reach significance." else ""
      small_zh <- if (min(s$n1, s$n2) < 10) "有一组少于 10 个样本，只有非常大的差异才可能显著。" else ""
      genes <- if (s$sig) paste0(": ", paste(utils::head(s$genes, 5), collapse = ", ")) else ""
      min_p <- if (is.na(s$min_p)) "-" else format(signif(s$min_p, 3))
      insight_bar(
        sprintf("<b>%d</b> of %d tested genes differ between <b>%s</b> (n = %s) and <b>%s</b> (n = %s) at FDR < %g (Fisher, BH)%s. Smallest raw p = %s.%s",
                s$sig, s$tested, res$l1, wes_fmt(s$n1), res$l2, wes_fmt(s$n2), res$fdr,
                genes, min_p, small_en),
        sprintf("%d 个受检基因中有 <b>%d</b> 个在 <b>%s</b>（n = %s）与 <b>%s</b>（n = %s）之间突变比例不同（FDR < %g，Fisher 检验，BH 校正）%s。最小原始 p = %s。%s",
                s$tested, s$sig, res$l1, wes_fmt(s$n1), res$l2, wes_fmt(s$n2), res$fdr,
                sub("^: ", "：", genes), min_p, small_zh))
    })

    forest_gg <- function() {
      shiny::req(res$cmp)
      s <- stats()
      wes_forest_plot(s$forest, res$l1, res$l2, s$n1, s$n2, res$fdr)
    }
    render_step_plot(output, input, "forest", forest_gg, name = "wes_forest",
                     width = 9, height = function() max(4, 2 + 0.35 * stats()$sig))

    draw_cobar <- function() {
      c0 <- res$cmp
      shiny::req(c0)
      s <- stats()
      if (!s$sig) {
        stop(sprintf("No gene passes FDR < %g, so there is nothing to compare side by side. The Results tab lists every tested gene.",
                     res$fdr))
      }
      wes_cobar_gg(c0$m1, c0$m2, res$l1, res$l2, s$genes)
    }
    render_step_plot(output, input, "cobar", draw_cobar, name = "wes_cobarplot", width = 10,
                     height = function() max(4.5, 2.2 + 0.4 * length(stats()$genes)))

    view <- shiny::reactive({
      shiny::req(res$cmp)
      df <- as.data.frame(res$cmp$res$results)
      num <- vapply(df, is.numeric, logical(1))
      df[num] <- lapply(df[num], function(x) signif(x, 3))
      if ("adjPval" %in% colnames(df)) df <- df[order(df$adjPval, df$pval), , drop = FALSE]
      df
    })
    tb <- wes_table(ns, "tbl", view, "wes_cohort_comparison",
                    has_maf = function() !is.null(rv$maf),
                    ready = function() !is.null(res$cmp))
    output$tbl_slot <- tb$slot
    output$tbl <- tb$table
    output$tbl_dl <- tb$download
  })
}
