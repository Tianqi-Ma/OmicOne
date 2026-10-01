#' WES module 10: Cohort comparison
#'
#' Split the cohort by a clinical column and test, gene by gene, whether the two
#' groups are mutated at different rates — the mutation-level equivalent of a
#' differential expression test.
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
      en = "Pick a clinical column and two of its levels; every gene mutated often
            enough is then tested with Fisher's exact test for a difference in
            mutation frequency between the two groups.",
      zh = "选择一个临床列及其中两个水平；随后对每个突变频率足够高的基因，用 Fisher 精确检验比较两组之间的突变频率差异。"),
    why  = list(
      en = "This is how you turn a cohort into a comparison: responders versus
            non-responders, primary versus metastatic, treated versus naive.
            The forest plot shows the odds ratio and its confidence interval per
            gene, so you can see effect size and not just a p-value.",
      zh = "这就是把一个队列变成一次比较的方法：应答 vs 不应答、原发 vs 转移、治疗过 vs 初治。森林图给出每个基因的比值比及其置信区间，让你看到效应量而不只是 p 值。"),
    how  = list(
      en = "Raise <b>minimum mutations</b> to avoid testing genes hit in one or
            two patients — with a few hundred genes tested, those produce
            spurious hits. The co-barplot next to the forest plot shows the raw
            frequencies the test was run on.",
      zh = "提高<b>最小突变数</b>，避免检验只在一两个患者中突变的基因——在检验几百个基因的情况下，这类基因很容易出现假阳性。森林图旁的并排柱状图展示了检验所依据的原始频率。"),
    example = list(
      en = "Comparing <code>FAB_classification</code> M0 against M3 in TCGA LAML
               recovers the known M3 pattern.",
      zh = "在 TCGA LAML 中比较 <code>FAB_classification</code> 的 M0 与 M3，可以还原已知的 M3 突变模式。")
  )
  controls <- shiny::tagList(
    shiny::uiOutput(ns("feat_ui")),
    shiny::uiOutput(ns("level_ui")),
    label_with_help("Minimum mutations per gene",
                    "A gene must be mutated at least this many times overall to be tested.",
                    label_zh = "每基因最少突变数",
                    tip_zh = "一个基因在整体上至少要有这么多突变才会被检验。"),
    shiny::numericInput(ns("min_mut"), NULL, value = 5, min = 2, max = 50, step = 1),
    label_with_help("p-value cutoff", "Genes below this are drawn on the forest plot.",
                    label_zh = "p 值阈值", tip_zh = "低于此值的基因会画在森林图上。"),
    shiny::numericInput(ns("pval"), NULL, value = 0.05, min = 0.001, max = 0.5, step = 0.01),
    run_button(ns("run"), "Compare cohorts", "比较队列")
  )
  step_container(
    title     = list(en = "Cohort comparison", zh = "队列比较"),
    explainer = explainer,
    controls  = controls,
    summary   = shiny::uiOutput(ns("summary")),
    preview   = bslib::navset_card_tab(
      bslib::nav_panel(i18n("Forest plot", "森林图"), preview_plot_ui(ns("forest"), download = TRUE)),
      bslib::nav_panel(i18n("Frequencies", "频率对比"), preview_plot_ui(ns("cobar"), download = TRUE)),
      bslib::nav_panel(i18n("Results", "结果表"), shiny::uiOutput(ns("tbl_slot")))
    )
  )
}

#' @rdname mod_wes_compare
#' @keywords internal
mod_wes_compare_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    ns  <- session$ns
    res <- shiny::reactiveValues(cmp = NULL, pval = 0.05, l1 = NULL, l2 = NULL)

    output$feat_ui <- shiny::renderUI({
      shiny::req(rv$maf)
      cols <- wes_clinical_cols(rv$maf)
      if (!length(cols)) {
        return(shiny::div(class = "omicone-status-empty",
                          i18n("This MAF has no clinical columns. Load a clinical table on the Import step.",
                               "该 MAF 无临床列。请在导入步骤加载临床表。")))
      }
      shiny::tagList(
        label_with_help("Split by", "The clinical column defining the two groups.",
                        label_zh = "分组依据", tip_zh = "用于定义两个分组的临床列。"),
        shiny::selectInput(ns("feature"), NULL, cols)
      )
    })

    feature_levels <- shiny::reactive({
      shiny::req(rv$maf, input$feature)
      cd <- as.data.frame(maftools::getClinicalData(rv$maf))
      v <- as.character(cd[[input$feature]])
      tab <- sort(table(v[!is.na(v) & nzchar(v)]), decreasing = TRUE)
      names(tab)[tab >= 2]
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
      if (identical(input$l1, input$l2)) {
        shiny::showNotification("Pick two different groups.", type = "error", duration = 8)
        return(NULL)
      }
      out <- with_progress_notify(
        wes_compare_cohorts(rv$maf, input$feature, input$l1, input$l2,
                            min_mut = max(2, as.integer(input$min_mut))),
        message = "Comparing cohorts...")
      if (is.null(out)) return(NULL)
      res$cmp  <- out
      res$pval <- input$pval
      res$l1   <- input$l1
      res$l2   <- input$l2
      mark_done(rv, "wes_compare")
      log_step(log_rv, "WES cohort comparison",
               params = list(feature = input$feature, group1 = input$l1,
                             group2 = input$l2, minMut = input$min_mut),
               code = sprintf('maftools::mafCompare(m1, m2, m1Name = "%s", m2Name = "%s", minMut = %s)',
                              input$l1, input$l2, input$min_mut))
    })

    output$summary <- shiny::renderUI({
      if (is.null(rv$maf)) return(wes_no_maf())
      c0 <- res$cmp
      if (is.null(c0)) {
        return(wes_prompt("Choose two groups and click <b>Compare cohorts</b>.",
                          "选择两个分组后点击<b>比较队列</b>。"))
      }
      n_sig <- tryCatch(sum(as.data.frame(c0$res$results)$pval <= res$pval, na.rm = TRUE),
                        error = function(e) NA_integer_)
      bslib::layout_columns(
        col_widths = c(3, 3, 3, 3),
        stat_tile(res$l1, format(c0$n1, big.mark = ",")),
        stat_tile(res$l2, format(c0$n2, big.mark = ",")),
        stat_tile(i18n("Genes tested", "受检基因数"),
                  tryCatch(nrow(as.data.frame(c0$res$results)), error = function(e) "-")),
        stat_tile(i18n("Significant", "显著基因数"), if (is.na(n_sig)) "-" else n_sig)
      )
    })

    draw_forest <- with_text_boost(function() {
      shiny::req(res$cmp)
      maftools::forestPlot(mafCompareRes = res$cmp$res, pVal = res$pval)
    })
    output$forest <- render_base_plot(draw_forest)
    register_figure_download(output, input, "forest", draw_forest,
                             "wes_forest", width = 9, height = 8)

    draw_cobar <- with_text_boost(function() {
      c0 <- res$cmp; shiny::req(c0)
      maftools::coBarplot(m1 = c0$m1, m2 = c0$m2, m1Name = res$l1, m2Name = res$l2)
    })
    output$cobar <- render_base_plot(draw_cobar)
    register_figure_download(output, input, "cobar", draw_cobar,
                             "wes_cobarplot", width = 10, height = 8)

    output$tbl_slot <- shiny::renderUI({
      if (is.null(res$cmp)) return(wes_no_maf())
      if (has_pkg("DT")) DT::dataTableOutput(ns("tbl"))
      else shiny::verbatimTextOutput(ns("tbl_txt"))
    })
    view <- shiny::reactive({
      shiny::req(res$cmp)
      df <- as.data.frame(res$cmp$res$results)
      num <- vapply(df, is.numeric, logical(1))
      df[num] <- lapply(df[num], function(x) signif(x, 3))
      if ("pval" %in% colnames(df)) df <- df[order(df$pval), , drop = FALSE]
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
