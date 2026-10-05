#' WES module 9: Clinical enrichment, pathways and drugs
#'
#' Three ways of asking "so what?": which genes are mutated more often in one
#' clinical group, which known oncogenic pathways the cohort's mutations fall
#' into, and which mutated genes fall into druggable categories.
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
      en = "<b>Enrichment</b> tests every frequently mutated gene in each level of
            one clinical variable against all other samples (Fisher's exact
            test, BH-adjusted across genes and levels). <b>Pathways</b> collapses
            genes into the ten canonical oncogenic signalling pathways
            (Sanchez-Vega 2018). <b>Drugs</b> looks the 20 most mutated genes up
            in DGIdb's druggable-gene categories.",
      zh = "<b>富集</b>针对某个临床变量的每个水平，把各高频突变基因与其余全部样本比较（Fisher 精确检验，跨基因与水平做 BH 校正）。<b>通路</b>把基因归入十条经典致癌信号通路（Sanchez-Vega 2018）。<b>药物</b>在 DGIdb 的可成药基因类别中检索突变最多的 20 个基因。"),
    why  = list(
      en = "A gene list is not a finding. Tying mutations to a clinical grouping,
            to a pathway, or to a drug category is what turns the cohort into
            something you can write about — with the multiple-testing burden
            stated.",
      zh = "一份基因列表本身不是结论。把突变与临床分组、通路或可成药类别联系起来，才能让这个队列变成可以写进文章的东西——同时要交代多重检验的负担。"),
    how  = list(
      en = "Enrichment needs a <b>clinical column</b> attached at import; only
            columns with 2–10 levels and at least 3 samples per level are
            offered. Samples with a missing value are left out. Genes must be
            mutated in more than <b>minimum mutated samples</b> to be tested.
            Pathways and drugs need no clinical data.",
      zh = "富集分析需要在导入时附带<b>临床列</b>；只列出有 2–10 个水平、且每个水平至少 3 个样本的列。取值缺失的样本不纳入。基因的突变样本数需大于<b>最少突变样本数</b>才会被检验。通路与药物分析不需要临床数据。"),
    read = list(
      en = "<b>Enrichment</b>: one bar pair per gene enriched in a group at your
            FDR (OR > 1 only): the upper bar is the share of that group's
            samples carrying a mutation, the lower bar the share in all other
            samples; whiskers are 95% binomial CIs and the labels give the raw
            counts. <b>Pathways</b>: per pathway, the fraction of its genes that
            are mutated and the fraction of samples with any mutated member.
            <b>Drugs</b>: DGIdb druggable categories among the top 20 genes —
            bar length = number of those genes in the category. This is gene
            druggability, not a list of approved drugs.",
      zh = "<b>富集</b>：每个在某组中达到 FDR 阈值且 OR > 1 的基因画一对柱：上方柱为该组样本中携带突变的比例，下方柱为其余样本中的比例；误差线为 95% 二项置信区间，标签给出原始计数。<b>通路</b>：每条通路中被突变的基因比例，以及有任一成员突变的样本比例。<b>药物</b>：前 20 个基因所属的 DGIdb 可成药类别——柱长＝属于该类别的基因数。这表示基因的可成药性，不是已获批药物的清单。"),
    example = list(
      en = "In TCGA LAML, enrichment on <code>FAB_classification</code> tests 8
               levels in the 192 samples with a FAB value; <code>IDH1</code> in M1
               (11 of 44 vs 7 of 148) and <code>TP53</code> in M7 (3 of 3) are the
               only hits at FDR < 0.05.",
      zh = "在 TCGA LAML 中，按 <code>FAB_classification</code> 做富集，在有 FAB 取值的 192 个样本中共检验 8 个水平；只有 M1 中的 <code>IDH1</code>（44 例中 11 例 vs 148 例中 7 例）与 M7 中的 <code>TP53</code>（3 例中 3 例）在 FDR < 0.05 时显著。")
  )
  controls <- shiny::tagList(
    shiny::uiOutput(ns("feat_ui")),
    label_with_help("FDR cutoff", "Enrichment results with a BH q-value below this are counted and drawn.",
                    label_zh = "FDR 阈值", tip_zh = "BH q 值低于此值的富集结果会被计数并绘制。"),
    shiny::numericInput(ns("fdr"), NULL, value = 0.05, min = 0.001, max = 0.25, step = 0.01),
    label_with_help("Minimum mutated samples", "A gene is tested only if more than this many samples carry it.",
                    label_zh = "最少突变样本数", tip_zh = "只有突变样本数大于此值的基因才会被检验。"),
    shiny::numericInput(ns("min_mut"), NULL, value = 5, min = 1, max = 50, step = 1),
    run_button(ns("run"), "Run analyses", "运行分析")
  )
  step_container(
    title     = list(en = "Clinical / pathway / drug", zh = "临床 / 通路 / 药物"),
    subtitle  = list(en = "Tie mutations to clinical groups, pathways, and druggable categories.",
                     zh = "把突变与临床分组、通路和可成药类别联系起来。"),
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
                                         caption = list(en = "Bars = share of samples with a mutation in the group (top) versus all other samples (bottom); whiskers = 95% binomial CI; genes at your FDR with OR > 1.",
                                                        zh = "柱＝该组中携带突变的样本比例（上）对其余样本（下）；误差线＝95% 二项 CI；仅含达到 FDR 阈值且 OR > 1 的基因。"))),
        bslib::nav_panel(i18n("Pathways", "通路"),
                         preview_plot_ui(ns("path"), download = TRUE,
                                         guide = list(en = "The ten oncogenic signalling pathways hit in the cohort will be drawn here.",
                                                      zh = "运行后，这里将绘制该队列命中的十条致癌信号通路。"),
                                         caption = list(en = "Per pathway: fraction of its genes mutated (n/N genes) and fraction of samples with any mutated member.",
                                                        zh = "每条通路：被突变基因的比例（n/N 个基因），以及有任一成员突变的样本比例。"))),
        bslib::nav_panel(i18n("Drugs", "药物"),
                         preview_plot_ui(ns("drug"), download = TRUE,
                                         guide = list(en = "DGIdb druggable categories of the top mutated genes will be drawn here.",
                                                      zh = "运行后，这里将绘制高频突变基因所属的 DGIdb 可成药类别。"),
                                         caption = list(en = "DGIdb druggable categories among the 20 most mutated genes; bar = number of genes in the category, label = up to five of them.",
                                                        zh = "突变最多的 20 个基因所属的 DGIdb 可成药类别；柱＝该类别中的基因数，标签＝其中最多五个基因。"))),
        bslib::nav_panel(i18n("Enrichment table", "富集结果表"), shiny::uiOutput(ns("tbl_slot")))
      ))
  )
}

#' @rdname mod_wes_clin
#' @keywords internal
mod_wes_clin_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    ns  <- session$ns
    res <- step_results(rv, "wes", enr = NULL, feature = NULL, fdr = 0.05, ran = FALSE,
                        path_err = NULL, drug_err = NULL)

    features <- shiny::reactive({
      shiny::req(rv$maf)
      cd <- tryCatch(as.data.frame(maftools::getClinicalData(rv$maf)),
                     error = function(e) NULL)
      if (is.null(cd)) character(0) else wes_enrichment_features(cd)
    })

    output$feat_ui <- shiny::renderUI({
      shiny::req(rv$maf)
      cols <- features()
      if (!length(cols)) {
        return(shiny::div(class = "omicone-status-empty",
                          i18n("No clinical column with 2–10 levels of at least 3 samples each — pathways and drugs still work.",
                               "没有水平数为 2–10 且每个水平至少 3 个样本的临床列——通路与药物分析仍可运行。")))
      }
      wes_col_select(ns, "feature",
                     label = list(en = "Clinical feature", zh = "临床变量"),
                     tip = list(en = "Columns with 2–10 levels, each carried by at least 3 samples.",
                                zh = "水平数为 2–10、且每个水平至少有 3 个样本的列。"),
                     choices = cols, selected = keep_selected(shiny::isolate(input$feature), cols))
    })

    shiny::observeEvent(input$run, {
      shiny::req(rv$maf)
      if (!require_pkgs("maftools", "Clinical / pathway / drug")) return(NULL)
      maf <- rv$maf
      feature <- input$feature
      if (!is.null(feature) && !feature %in% features()) feature <- NULL
      fdr <- wes_prob(input$fdr, 0.05)
      min_mut <- wes_int(input$min_mut, 5, lo = 1, hi = 50)
      enr <- NULL
      if (!is.null(feature) && nzchar(feature)) {
        enr <- with_progress_notify(wes_clinical_enrichment(maf, feature, min_mut = min_mut),
                                    message = "Testing clinical enrichment...")
        if (is.null(enr)) return(NULL)
      }
      path_err <- tryCatch({
        wes_dry_run(function() wes_pathways(maf))
        NULL
      }, error = function(e) conditionMessage(e))
      drug_err <- tryCatch({
        wes_dry_run(function() maftools::drugInteractions(maf = maf, fontSize = 0.95))
        NULL
      }, error = function(e) conditionMessage(e))
      res$enr <- enr
      res$feature <- if (!is.null(enr)) feature else NULL
      res$fdr <- fdr
      res$path_err <- path_err
      res$drug_err <- drug_err
      res$ran <- TRUE
      mark_done(rv, "wes_clin")
      path_code <- if (exists("pathways", where = asNamespace("maftools"), inherits = FALSE)) {
        wes_code("maftools::pathways", list(maf = quote(maf), plotType = "bar"))
      } else {
        wes_code("maftools:::OncogenicPathways", list(maf = quote(maf)))
      }
      enr_code <- if (!is.null(enr)) {
        c("cd <- as.data.frame(maftools::getClinicalData(maf))",
          sprintf("anno <- cd[!is.na(cd[[%1$s]]) & nzchar(trimws(as.character(cd[[%1$s]]))), c(\"Tumor_Sample_Barcode\", %1$s)]",
                  deparse(feature)),
          wes_code("maftools::clinicalEnrichment",
                   list(maf = quote(maf), clinicalFeature = feature, annotationDat = quote(anno),
                        minMut = min_mut), assign = "enr"),
          sprintf("enr$groupwise_comparision[fdr < %s]   # BH across genes x levels", fdr),
          "enr_sig <- enr",
          sprintf("enr_sig$groupwise_comparision <- enr$groupwise_comparision[fdr < %s]", fdr),
          wes_code("maftools::plotEnrichmentResults",
                   list(enrich_res = quote(enr_sig), pVal = 0.05)))
      }
      log_step(log_rv, "WES clinical / pathway / drug",
               params = list(feature = feature %||% "(none)", fdr = fdr, minMut = min_mut),
               code = c(enr_code, path_code,
                        wes_code("maftools::drugInteractions", list(maf = quote(maf)))))
    })

    stats <- shiny::reactive({
      shiny::req(res$ran)
      if (is.null(res$enr)) return(NULL)
      sig <- wes_enrichment_sig(res$enr, res$fdr)
      sig$n_used <- attr(res$enr, "n_used")
      sig$n_missing <- attr(res$enr, "n_missing")
      sig$levels <- tryCatch(nrow(res$enr$cf_sizes), error = function(e) NA_integer_)
      sig
    })

    output$summary <- shiny::renderUI({
      if (is.null(rv$maf)) return(wes_no_maf())
      if (!isTRUE(res$ran)) {
        return(wes_prompt("Click <b>Run analyses</b>.", "点击<b>运行分析</b>。"))
      }
      s <- stats()
      shiny::tagList(
        stat_tile(i18n("Feature", "临床变量"), res$feature %||% "-"),
        stat_tile(i18n("Genes tested", "受检基因数"), if (is.null(s)) "-" else s$tested_genes),
        stat_tile(sprintf("FDR < %g", res$fdr), if (is.null(s)) "-" else s$sig_genes),
        stat_tile(i18n("Samples used", "纳入样本数"), if (is.null(s)) "-" else wes_fmt(s$n_used))
      )
    })

    output$insight <- shiny::renderUI({
      if (!isTRUE(res$ran)) return(NULL)
      if (is.null(res$enr)) {
        return(insight_bar(
          "Pathway and drug look-ups are done; no clinical feature was selected, so enrichment testing was skipped.",
          "通路与药物检索已完成；未选择临床变量，因此跳过了富集检验。"))
      }
      s <- stats()
      miss_en <- if (isTRUE(s$n_missing > 0)) sprintf(" %d sample(s) with a missing %s were left out.",
                                                      s$n_missing, res$feature) else ""
      miss_zh <- if (isTRUE(s$n_missing > 0)) sprintf("%d 个 %s 缺失的样本未纳入。",
                                                      s$n_missing, res$feature) else ""
      genes_txt <- if (length(s$genes)) {
        paste0(": ", paste(utils::head(s$genes, 5), collapse = ", "))
      } else ""
      n_tests <- nrow(wes_enrichment_table(res$enr))
      insight_bar(
        sprintf("<b>%d</b> of %d tested genes are enriched in at least one level of <b>%s</b> at FDR < %g (BH over %d gene × level tests)%s. Each level is compared with all other samples; small levels give wide CIs.%s",
                s$sig_genes, s$tested_genes, res$feature, res$fdr, n_tests, genes_txt, miss_en),
        sprintf("%d 个受检基因中有 <b>%d</b> 个在 <b>%s</b> 的至少一个水平中富集（FDR < %g，对 %d 个基因 × 水平检验做 BH 校正）%s。每个水平均与其余全部样本比较；样本少的水平置信区间很宽。%s",
                s$tested_genes, s$sig_genes, res$feature, res$fdr, n_tests,
                sub("^: ", "：", genes_txt), miss_zh))
    })

    draw_enr <- with_text_boost(function() {
      shiny::req(res$ran)
      if (is.null(res$enr)) {
        stop("No clinical feature selected. Attach a clinical table on the Import step to use this tab.")
      }
      s <- stats()
      if (!nrow(s$drawn)) {
        stop(sprintf("No gene is enriched (OR > 1) in any level at FDR < %g. The Enrichment table lists every test.",
                     res$fdr))
      }
      sub <- res$enr
      sub$groupwise_comparision <- data.table::as.data.table(s$sig)
      # pVal fixed at 0.05: maftools also uses it as the CI alpha, so the
      # whiskers stay 95% whatever FDR cutoff the user chose
      maftools::plotEnrichmentResults(enrich_res = sub, pVal = 0.05)
    })
    output$enr <- render_base_plot(draw_enr)
    register_figure_download(output, input, "enr", draw_enr, "wes_enrichment",
                             width = 10, height = 7)

    draw_path <- with_text_boost(function() {
      shiny::req(rv$maf, res$ran)
      if (!is.null(res$path_err)) stop(res$path_err)
      wes_pathways(rv$maf)
    })
    output$path <- render_base_plot(draw_path)
    register_figure_download(output, input, "path", draw_path, "wes_pathways",
                             width = 10, height = 7)

    draw_drug <- with_text_boost(function() {
      shiny::req(rv$maf, res$ran)
      if (!is.null(res$drug_err)) stop(res$drug_err)
      invisible(maftools::drugInteractions(maf = rv$maf, fontSize = 0.95))
    })
    output$drug <- render_base_plot(draw_drug)
    register_figure_download(output, input, "drug", draw_drug, "wes_drugs",
                             width = 10, height = 8)

    view <- shiny::reactive({
      shiny::req(res$enr)
      df <- wes_enrichment_table(res$enr)
      num <- vapply(df, is.numeric, logical(1))
      df[num] <- lapply(df[num], function(x) signif(x, 3))
      df
    })
    tb <- wes_table(ns, "tbl", view, "wes_clinical_enrichment",
                    has_maf = function() !is.null(rv$maf),
                    ready = function() !is.null(res$enr),
                    not_ready = list(en = "Select a clinical feature and run this step to show the enrichment table.",
                                     zh = "选择临床变量并运行本步骤后显示富集结果表。"))
    output$tbl_slot <- tb$slot
    output$tbl <- tb$table
    output$tbl_dl <- tb$download
  })
}
