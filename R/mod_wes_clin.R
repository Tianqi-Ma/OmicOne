#' WES module 9: Clinical enrichment, pathways and drugs
#'
#' Three ways of asking "so what?": which genes are mutated more often in one
#' clinical group, which known oncogenic pathways the cohort's mutations fall
#' into, and which mutated genes are druggable.
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_wes_clin
NULL

#' @rdname mod_wes_clin
#' @keywords internal
mod_wes_clin_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = list(en = "Clinical / pathway / drug", zh = "临床 / 通路 / 药物"),
    what = list(
      en = "<b>Enrichment</b> tests every gene against the levels of one clinical
            variable. <b>Pathways</b> collapses genes into the ten canonical
            oncogenic pathways. <b>Drugs</b> looks the mutated genes up in a
            drug-gene interaction database.",
      zh = "<b>富集</b>针对某个临床变量的各个水平检验每个基因。<b>通路</b>把基因归入十条经典致癌通路。<b>药物</b>在药物-基因相互作用数据库中检索这些突变基因。"),
    why  = list(
      en = "A gene list is not a finding. Tying mutations to a clinical grouping,
            to a pathway, or to a drug is what turns the cohort into something you
            can write about or act on.",
      zh = "一份基因列表本身不是结论。把突变与临床分组、与通路、或与药物联系起来，才能让这个队列变成可以写进文章、或可据以决策的东西。"),
    how  = list(
      en = "Enrichment needs a <b>clinical column</b> attached at import — if the
            dropdown is empty, load a clinical table on the Import step. Pathways
            and drugs need no clinical data at all.",
      zh = "富集分析需要在导入时附带<b>临床列</b>——如果下拉框是空的，请在导入步骤加载临床表。通路和药物分析不需要任何临床数据。"),
    read = list(
      en = "<b>Enrichment</b>: one bar per gene whose mutation frequency differs
            across the clinical groups — the longer the bar, the smaller the
            p-value; bars past the line survive your cutoff. <b>Pathways</b>:
            the share of samples hitting each of the ten canonical oncogenic
            pathways, counting any member gene. <b>Drugs</b>: mutated genes with
            known drug interactions, ranked by how many compounds target them —
            the top of that list is where repurposing hypotheses start.",
      zh = "<b>富集</b>：每根条代表一个突变频率在临床分组间有差异的基因——条越长 p 值越小；超过竖线的基因通过了你的阈值。<b>通路</b>：命中十条经典致癌通路各自的样本占比（任一成员基因命中即计入）。<b>药物</b>：具有已知药物相互作用的突变基因，按对应化合物数量排序——排在前列的就是老药新用假设的起点。"),
    example = list(
      en = "In TCGA LAML, enrichment on <code>FAB_classification</code> recovers
               the M3-specific mutation pattern without being told which samples
               are M3.",
      zh = "在 TCGA LAML 中，按 <code>FAB_classification</code> 做富集，无需事先告知哪些样本是 M3，就能还原出 M3 特有的突变模式。")
  )
  controls <- shiny::tagList(
    shiny::uiOutput(ns("feat_ui")),
    label_with_help("Significance cutoff", "Enrichment results at or below this p-value are drawn.",
                    label_zh = "显著性阈值", tip_zh = "p 值不高于此值的富集结果会被绘制。"),
    shiny::numericInput(ns("pval"), NULL, value = 0.05, min = 0.001, max = 0.5, step = 0.01),
    run_button(ns("run"), "Run analyses", "运行分析")
  )
  step_container(
    title     = list(en = "Clinical / pathway / drug", zh = "临床 / 通路 / 药物"),
    subtitle  = list(en = "Tie mutations to clinical groups, pathways, and drugs.",
                     zh = "把突变与临床分组、通路和药物联系起来。"),
    explainer = explainer,
    controls  = controls,
    summary   = shiny::uiOutput(ns("summary")),
    preview   = shiny::tagList(
      shiny::uiOutput(ns("insight")),
      bslib::navset_card_tab(
      bslib::nav_panel(i18n("Clinical enrichment", "临床富集"),
                       preview_plot_ui(ns("enr"), download = TRUE,
                                       guide = list(en = "Genes enriched in one clinical group will be drawn here.",
                                                    zh = "运行后，这里将绘制在某个临床分组中富集的基因。"),
                                       caption = list(en = "Bar length ∝ −log10(p): how strongly each gene's mutation rate differs across clinical groups.",
                                                      zh = "条长 ∝ −log10(p)：每个基因的突变频率在临床分组间的差异强度。"))),
      bslib::nav_panel(i18n("Pathways", "通路"),   preview_plot_ui(ns("path"), download = TRUE)),
      bslib::nav_panel(i18n("Drugs", "药物"),      preview_plot_ui(ns("drug"), download = TRUE)),
      bslib::nav_panel(i18n("Enrichment table", "富集结果表"), shiny::uiOutput(ns("tbl_slot")))
    ))
  )
}

#' @rdname mod_wes_clin
#' @keywords internal
mod_wes_clin_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    ns  <- session$ns
    res <- shiny::reactiveValues(enr = NULL, feature = NULL, pval = 0.05, ran = FALSE)

    output$feat_ui <- shiny::renderUI({
      shiny::req(rv$maf)
      cols <- wes_clinical_cols(rv$maf)
      if (!length(cols)) {
        return(shiny::div(class = "omicone-status-empty",
                          i18n("No clinical columns in this MAF — pathways and drugs still work.",
                               "该 MAF 无临床列——通路与药物分析仍可运行。")))
      }
      shiny::tagList(
        label_with_help("Clinical feature", "The grouping to test each gene against.",
                        label_zh = "临床变量", tip_zh = "用于检验每个基因的分组变量。"),
        shiny::selectInput(ns("feature"), NULL, cols)
      )
    })

    shiny::observeEvent(input$run, {
      shiny::req(rv$maf)
      if (!require_pkgs("maftools", "Clinical / pathway / drug")) return(NULL)
      res$pval <- input$pval
      res$ran  <- TRUE
      res$feature <- input$feature
      if (!is.null(input$feature) && nzchar(input$feature)) {
        res$enr <- with_progress_notify(
          wes_clinical_enrichment(rv$maf, input$feature),
          message = "Testing clinical enrichment...")
      } else {
        res$enr <- NULL
      }
      mark_done(rv, "wes_clin")
      log_step(log_rv, "WES clinical / pathway / drug",
               params = list(feature = input$feature %||% "(none)", pvalue = input$pval),
               code = c(if (!is.null(input$feature) && nzchar(input$feature))
                          sprintf('enr <- maftools::clinicalEnrichment(maf, clinicalFeature = "%s")',
                                  input$feature),
                        'maftools::pathways(maf)',
                        'maftools::drugInteractions(maf)'))
    })

    output$summary <- shiny::renderUI({
      if (is.null(rv$maf)) return(wes_no_maf())
      if (!isTRUE(res$ran)) {
        return(wes_prompt("Click <b>Run analyses</b>.", "点击<b>运行分析</b>。"))
      }
      n_sig <- tryCatch({
        d <- as.data.frame(res$enr$groupwise_comparision)
        sum(d$p_value <= res$pval, na.rm = TRUE)
      }, error = function(e) NA_integer_)
      bslib::layout_columns(
        col_widths = c(4, 4, 4),
        stat_tile(i18n("Feature", "临床变量"), res$feature %||% "-"),
        stat_tile(i18n("Enriched genes", "富集基因数"),
                  if (is.na(n_sig)) "-" else n_sig),
        stat_tile(i18n("p cutoff", "p 阈值"), res$pval)
      )
    })

    output$insight <- shiny::renderUI({
      if (!isTRUE(res$ran)) return(NULL)
      if (is.null(res$enr)) {
        return(insight_bar(
          "Pathway and drug analyses are done; no clinical feature was attached, so enrichment testing was skipped.",
          "通路与药物分析已完成；因未附带临床变量，富集检验被跳过。"))
      }
      n_sig <- tryCatch({
        d <- as.data.frame(res$enr$groupwise_comparision)
        sum(d$p_value <= res$pval, na.rm = TRUE)
      }, error = function(e) NA_integer_)
      if (is.na(n_sig)) return(NULL)
      insight_bar(
        sprintf("<b>%d</b> genes differ across <b>%s</b> at p ≤ %g. The <b>Pathways</b> tab shows which canonical pathways the cohort hits; <b>Drugs</b> ranks repurposing candidates.",
                n_sig, res$feature, res$pval),
        sprintf("<b>%d</b> 个基因在 <b>%s</b> 各组间的突变频率不同（p ≤ %g）。<b>通路</b>页签展示队列命中的经典通路；<b>药物</b>页签列出老药新用候选。",
                n_sig, res$feature, res$pval))
    })

    draw_enr <- with_text_boost(function() {
      shiny::req(res$ran)
      if (is.null(res$enr)) {
        stop("No clinical feature selected. Attach a clinical table on the Import step to use this tab.")
      }
      maftools::plotEnrichmentResults(enrich_res = res$enr, pVal = res$pval)
    })
    output$enr <- render_base_plot(draw_enr)
    register_figure_download(output, input, "enr", draw_enr, "wes_enrichment",
                             width = 10, height = 7)

    draw_path <- with_text_boost(function() {
      shiny::req(rv$maf, res$ran)
      wes_pathways(rv$maf)
    })
    output$path <- render_base_plot(draw_path)
    register_figure_download(output, input, "path", draw_path, "wes_pathways",
                             width = 10, height = 7)

    draw_drug <- with_text_boost(function() {
      shiny::req(rv$maf, res$ran)
      maftools::drugInteractions(maf = rv$maf, fontSize = 0.95)
    })
    output$drug <- render_base_plot(draw_drug)
    register_figure_download(output, input, "drug", draw_drug, "wes_drugs",
                             width = 10, height = 8)

    output$tbl_slot <- shiny::renderUI({
      if (is.null(res$enr)) return(wes_no_maf())
      if (has_pkg("DT")) DT::dataTableOutput(ns("tbl"))
      else shiny::verbatimTextOutput(ns("tbl_txt"))
    })
    view <- shiny::reactive({
      shiny::req(res$enr)
      df <- as.data.frame(res$enr$groupwise_comparision)
      num <- vapply(df, is.numeric, logical(1))
      df[num] <- lapply(df[num], function(x) signif(x, 3))
      if ("p_value" %in% colnames(df)) df <- df[order(df$p_value), , drop = FALSE]
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
