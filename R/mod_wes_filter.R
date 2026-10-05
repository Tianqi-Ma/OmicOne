#' Module: WES variant filters
#'
#' Remove calls a careful analysis would not keep (FILTER != PASS, low depth,
#' few supporting reads, low allele fraction, common germline polymorphisms),
#' count what each filter removed, and flag hypermutated samples. Every later
#' WES step reads the filtered MAF. The step always starts from the imported
#' MAF, so re-running it with new thresholds never filters twice.
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_wes_filter
NULL

#' @rdname mod_wes_filter
#' @keywords internal
mod_wes_filter_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = list(en = "Variant filters", zh = "变异过滤"),
    what = list(
      en = "Drop low-confidence calls from the MAF and see how many each filter removes; flag samples
            with an extreme number of mutations.",
      zh = "从 MAF 中去掉置信度低的突变，查看每个过滤条件去掉了多少；并标记突变数异常多的样本。"),
    why = list(
      en = "Frequencies, oncoplots, TMB and signatures are only as good as the calls behind them. A
            MAF exported straight from a caller can still contain calls the caller itself flagged,
            calls supported by two reads, and germline polymorphisms that a missing or thin normal
            let through. Hypermutated samples (MMR or POLE deficiency, some treatments, artefacts)
            can dominate every cohort-level result.",
      zh = "突变频率、oncoplot、TMB 和突变特征的可靠性取决于背后的突变检测结果。直接从检测工具导出的 MAF 里，可能仍有工具自己标记为不通过的突变、只有两条 reads 支持的突变，以及因为没有正常对照或对照太浅而漏过来的胚系多态性。超突变样本（MMR 或 POLE 缺陷、某些治疗、技术假象）会主导所有队列层面的结果。"),
    how = list(
      en = "Filters run in this order; leave a box empty to switch a filter off. Typical values:
            depth &ge; 10–20, alt reads &ge; 3–5, VAF &ge; 0.05 (lower for subclonal or low-purity
            work), population AF &le; 0.001. A filter whose column is missing from the MAF is skipped
            and reported. A variant with an empty value is kept (the filter cannot judge it).",
      zh = "过滤按这里的顺序进行；某项留空即关闭该过滤。常用取值：深度 &ge; 10–20，支持 reads &ge; 3–5，VAF &ge; 0.05（研究亚克隆或低纯度样本时可更低），人群频率 &le; 0.001。MAF 中缺少对应列的过滤会被跳过并注明。某个值为空的突变会被保留（无法判断）。"),
    read = list(
      en = "The funnel shows non-synonymous variants left after each filter, with the number removed.
            The VAF tab shows which part of the distribution was cut. The hypermutator tab sorts
            samples by mutation count; points above the dashed fence are flagged. Flagged samples
            stay in: rerun key analyses without them as a sensitivity check.",
      zh = "漏斗图显示每一步过滤后剩下的非同义突变数，以及被去掉的数量。VAF 页显示分布中被切掉的部分。超突变页按突变数给样本排序，虚线以上的点被标记。被标记的样本仍保留在数据中：可去掉它们重做关键分析，作为敏感性检验。"),
    example = list(
      en = "The TCGA-LAML demo carries only a VAF column: VAF &ge; 0.05 removes 26 of its 1,732
            non-synonymous calls, the other filters are skipped, and no sample is hypermutated.",
      zh = "TCGA-LAML 演示数据只有 VAF 一列：VAF &ge; 0.05 去掉 1,732 个非同义突变中的 26 个，其余过滤被跳过，没有样本被标记为超突变。")
  )
  controls <- shiny::tagList(
    shiny::uiOutput(ns("cols_ui")),
    shiny::checkboxInput(ns("pass"), i18n("Keep only FILTER = PASS", "仅保留 FILTER = PASS"), value = TRUE),
    label_with_help("Min tumour depth", "Reads covering the position in the tumour (t_depth, or t_ref_count + t_alt_count).",
                    "最小肿瘤深度", "肿瘤中覆盖该位置的 reads 数（t_depth，或 t_ref_count + t_alt_count）。"),
    shiny::numericInput(ns("min_depth"), NULL, value = 10, min = 0, step = 1),
    label_with_help("Min alt reads", "Reads supporting the variant (t_alt_count).",
                    "最少支持 reads", "支持该突变的 reads 数（t_alt_count）。"),
    shiny::numericInput(ns("min_alt"), NULL, value = 3, min = 0, step = 1),
    label_with_help("Min VAF", "Variant allele fraction in the tumour, as a fraction (0.05 = 5%). A column in percent is detected and rescaled.",
                    "最小 VAF", "肿瘤中的变异等位基因频率，用小数表示（0.05 = 5%）。百分数形式的列会被自动识别并换算。"),
    shiny::numericInput(ns("min_vaf"), NULL, value = 0.05, min = 0, max = 1, step = 0.01),
    label_with_help("Max population AF", "Frequency in gnomAD / ExAC / 1000 Genomes. Somatic calls above 0.001 are most likely germline. A variant absent from the database counts as 0.",
                    "最大人群频率", "在 gnomAD / ExAC / 千人基因组中的频率。体细胞突变若高于 0.001，多半是胚系变异。数据库中没有的变异按 0 计。"),
    shiny::numericInput(ns("max_pop"), NULL, value = 0.001, min = 0, max = 1, step = 0.001),
    label_with_help("Max alt reads in the normal", "Reads supporting the variant in the matched normal (n_alt_count); a somatic call should have almost none. Empty = off.",
                    "正常样本中最多支持 reads", "配对正常样本中支持该突变的 reads 数（n_alt_count）；体细胞突变在正常样本中应几乎没有。留空 = 关闭。"),
    shiny::numericInput(ns("max_nalt"), NULL, value = NA, min = 0, step = 1),
    run_button(ns("run"), "Apply filters", "应用过滤")
  )
  step_container(
    title     = list(en = "Variant filters", zh = "变异过滤"),
    subtitle  = list(en = "Keep confident somatic calls; see what each filter removes.",
                     zh = "保留可信的体细胞突变；查看每个过滤条件去掉了什么。"),
    explainer = explainer,
    controls  = controls,
    summary   = shiny::uiOutput(ns("summary")),
    preview   = shiny::tagList(
      shiny::uiOutput(ns("insight")),
      bslib::navset_card_tab(
        bslib::nav_panel(i18n("Funnel", "漏斗"),
          preview_plot_ui(ns("funnel"), download = TRUE, scene = "wes_filter",
            guide = list(en = "The number of variants left after each filter will be drawn here.",
                         zh = "运行后，这里将绘制每一步过滤后剩下的突变数。"),
            caption = list(en = "Non-synonymous variants left after each filter, in the order applied; (-n) = removed by that filter.",
                           zh = "按过滤顺序，每一步之后剩下的非同义突变数；（-n）= 该步去掉的数量。"))),
        bslib::nav_panel(i18n("VAF", "VAF"),
          preview_plot_ui(ns("vaf"), download = TRUE,
            caption = list(en = "Tumour VAF of the non-synonymous variants: kept vs removed by any filter; dashed = the VAF threshold.",
                           zh = "非同义突变的肿瘤 VAF：保留与被任一过滤去掉的；虚线 = VAF 阈值。"))),
        bslib::nav_panel(i18n("Hypermutators", "超突变"),
          preview_plot_ui(ns("hyper"), download = TRUE,
            caption = list(en = "Non-synonymous variants per sample after filtering (log scale); above the dashed fence = flagged.",
                           zh = "过滤后每个样本的非同义突变数（对数坐标）；虚线以上 = 被标记。"))),
        bslib::nav_panel(i18n("Table", "表格"), shiny::uiOutput(ns("tbl_slot")))
      )
    )
  )
}

#' @rdname mod_wes_filter
#' @keywords internal
mod_wes_filter_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    res <- step_results(rv, "wes", out = NULL, hyper = NULL)

    # the imported MAF: this step never filters an already filtered MAF
    source_maf <- shiny::reactive({
      shiny::req(rv$maf)
      step_input(rv, "wes_filter", rv$maf)
    })

    output$cols_ui <- shiny::renderUI({
      if (is.null(rv$maf)) return(wes_no_maf())
      cols <- wes_filter_cols(wes_fields(source_maf()))
      row <- function(en, zh, col) {
        shiny::tags$li(i18n(en, zh), ": ",
                       if (is.null(col)) shiny::tags$span(class = "omicone-muted", i18n("not in this MAF", "MAF 中没有"))
                       else shiny::tags$code(col))
      }
      shiny::div(class = "omicone-sheet",
        shiny::tags$b(i18n("Columns found", "识别到的列")),
        shiny::tags$ul(class = "omicone-help-text",
          row("FILTER", "FILTER", cols$filter),
          row("Depth", "深度", cols$depth %||% if (!is.null(cols$ref) && !is.null(cols$alt)) paste(cols$ref, "+", cols$alt)),
          row("Alt reads", "支持 reads", cols$alt),
          row("VAF", "VAF", cols$vaf %||% if (!is.null(cols$alt) && (!is.null(cols$depth) || !is.null(cols$ref))) "alt / depth"),
          row("Population AF", "人群频率", cols$pop),
          row("Normal alt reads", "正常样本支持 reads", cols$normal_alt)))
    })

    shiny::observeEvent(input$run, {
      shiny::req(rv$maf)
      if (!require_pkgs("maftools", "Variant filters")) return(NULL)
      th <- list(pass = isTRUE(input$pass),
                 min_depth = num_input(input$min_depth, 0),
                 min_alt = num_input(input$min_alt, 0),
                 min_vaf = num_input(input$min_vaf, 0, 1),
                 max_pop = num_input(input$max_pop, 0, 1),
                 max_normal_alt = num_input(input$max_nalt, 0))
      maf_in <- source_maf()
      out <- with_progress_notify(
        do.call(wes_filter_apply, c(list(maf = maf_in), th)),
        message = "Filtering variants...")
      if (is.null(out)) return(NULL)
      hm <- tryCatch(wes_hypermutators(out$maf), error = function(e) NULL)
      res$out <- out
      res$hyper <- hm
      rv$maf <- out$maf
      mark_done(rv, "wes_filter")
      fun <- out$funnel
      log_step(log_rv, "WES variant filters",
               params = c(Filter(function(v) !(is.numeric(v) && is.na(v)), out$thresholds),
                          list(skipped = if (length(out$unavailable)) paste(out$unavailable, collapse = ", ") else "none",
                               nonsyn_before = fun$nonsyn_left[1], nonsyn_after = fun$nonsyn_left[nrow(fun)])),
               code = wes_filter_code(out))
      wes_notify(sprintf("Kept %s of %s non-synonymous variants.", wes_fmt(fun$nonsyn_left[nrow(fun)]),
                         wes_fmt(fun$nonsyn_left[1])),
                 sprintf("保留了 %s / %s 个非同义突变。", wes_fmt(fun$nonsyn_left[nrow(fun)]),
                         wes_fmt(fun$nonsyn_left[1])), type = "message", duration = 5)
    })

    output$summary <- shiny::renderUI({
      if (is.null(rv$maf)) return(wes_no_maf())
      out <- res$out
      if (is.null(out)) {
        return(wes_prompt("Set the thresholds and click <b>Apply filters</b> (optional step: skip it to use the MAF as imported).",
                          "设置阈值后点击<b>应用过滤</b>（可选步骤：跳过则使用导入时的 MAF）。"))
      }
      f <- out$funnel
      hm <- res$hyper
      shiny::tagList(
        stat_tile(i18n("Non-synonymous kept", "保留的非同义突变"),
                  sprintf("%s / %s", wes_fmt(f$nonsyn_left[nrow(f)]), wes_fmt(f$nonsyn_left[1]))),
        stat_tile(i18n("Filters applied", "应用的过滤"), nrow(f) - 1),
        stat_tile(i18n("Samples with no variant left", "无突变剩余的样本"), length(out$lost)),
        stat_tile(i18n("Hypermutated (flag)", "超突变（标记）"),
                  if (is.null(hm)) "-" else sum(hm$hypermutated))
      )
    })

    output$insight <- shiny::renderUI({
      out <- res$out
      if (is.null(out)) return(NULL)
      f <- out$funnel
      hm <- res$hyper
      big <- f[-1, , drop = FALSE]
      big <- if (nrow(big)) big[which.max(big$removed_nonsyn), ] else NULL
      skipped <- out$unavailable
      lab <- c(pass = "FILTER", min_depth = "depth", min_alt = "alt reads", min_vaf = "VAF",
               max_pop = "population AF", max_normal_alt = "normal alt reads")
      nh <- if (is.null(hm)) 0 else sum(hm$hypermutated)
      insight_bar(
        sprintf("%s of %s non-synonymous variants kept (%.1f%%)%s.%s%s%s",
                wes_fmt(f$nonsyn_left[nrow(f)]), wes_fmt(f$nonsyn_left[1]),
                100 * f$nonsyn_left[nrow(f)] / f$nonsyn_left[1],
                if (!is.null(big) && big$removed_nonsyn > 0) sprintf("; the %s filter removed the most (%s)", big$filter, wes_fmt(big$removed_nonsyn)) else "",
                if (length(skipped)) sprintf(" Skipped, column not in the MAF: %s.", paste(lab[skipped], collapse = ", ")) else "",
                if (length(out$lost)) sprintf(" %d sample(s) have no variant left: still sequenced, so wild-type for every gene and TMB 0.", length(out$lost)) else "",
                if (nh) sprintf(" %d hypermutated sample(s) flagged (kept): repeat key analyses without them.", nh) else " No sample is hypermutated."),
        sprintf("保留了 %s / %s 个非同义突变（%.1f%%）%s。%s%s%s",
                wes_fmt(f$nonsyn_left[nrow(f)]), wes_fmt(f$nonsyn_left[1]),
                100 * f$nonsyn_left[nrow(f)] / f$nonsyn_left[1],
                if (!is.null(big) && big$removed_nonsyn > 0) sprintf("；去掉最多的是 %s（%s 个）", big$filter, wes_fmt(big$removed_nonsyn)) else "",
                if (length(skipped)) sprintf("因 MAF 中缺少对应列而跳过：%s。", paste(lab[skipped], collapse = "、")) else "",
                if (length(out$lost)) sprintf("%d 个样本已没有任何突变：它们仍是已测序样本，所有基因都算野生型、TMB 为 0。", length(out$lost)) else "",
                if (nh) sprintf("标记了 %d 个超突变样本（仍保留）：建议去掉它们重做关键分析。", nh) else "没有超突变样本。"))
    })

    render_step_plot(output, input, "funnel", function() {
      shiny::req(res$out)
      wes_filter_funnel_plot(res$out$funnel)
    }, name = "wes_filter_funnel", width = 9, height = 5)

    render_step_plot(output, input, "vaf", function() {
      shiny::req(res$out)
      wes_filter_vaf_plot(res$out$values, res$out$thresholds$min_vaf)
    }, name = "wes_filter_vaf", width = 9, height = 5.5)

    render_step_plot(output, input, "hyper", function() {
      shiny::req(res$hyper)
      wes_hyper_plot(res$hyper)
    }, name = "wes_hypermutators", width = 9, height = 5.5)

    view <- shiny::reactive({
      out <- res$out
      shiny::req(out)
      f <- out$funnel
      f[, c("step", "filter", "column", "removed", "removed_nonsyn", "not_evaluable", "nonsyn_left",
            "all_left", "samples_left")]
    })
    tb <- wes_table(ns, "tbl", view, "wes_filter_funnel",
                    has_maf = function() !is.null(rv$maf),
                    ready = function() !is.null(res$out))
    output$tbl_slot <- tb$slot
    output$tbl <- tb$table
    output$tbl_dl <- tb$download
  })
}
