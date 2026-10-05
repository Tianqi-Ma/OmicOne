#' Module: Pseudobulk differential expression
#'
#' Compare conditions (responders vs non-responders, treated vs untreated)
#' within each cell type, with the sample as the unit of replication: counts of
#' one cell type are summed per sample, then tested across samples with edgeR
#' (quasi-likelihood), limma-voom or DESeq2. This is the between-condition test;
#' the Markers step is cell-level marker discovery and must not be used for it.
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_pseudobulk
NULL

#' @rdname mod_pseudobulk
#' @keywords internal
mod_pseudobulk_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = list(en = "Pseudobulk differential expression", zh = "Pseudobulk 差异表达"),
    what = list(
      en = "For each cell type, add up the counts of its cells within every sample, then compare the
            samples of two conditions with a bulk RNA-seq model.",
      zh = "对每种细胞类型，先把每个样本内该类型细胞的计数加起来，再用 bulk RNA-seq 模型比较两种条件下的样本。"),
    why = list(
      en = "Cells from one patient are not independent replicates. Testing thousands of cells as if
            they were gives tiny, misleading p-values, and the result mostly reflects which patient has
            the most cells. Summing per sample makes the patient the replicate, which is what a claim
            about responders or treatment needs (Squair et al. 2021, <i>Nat Commun</i>).",
      zh = "来自同一患者的细胞不是独立重复。把成千上万个细胞当作重复来检验，会得到极小且具有误导性的 p 值，结果主要反映哪个患者细胞多。按样本求和后，患者才是重复单位，这才能支持关于疗效或处理效应的结论（Squair 等 2021，<i>Nat Commun</i>）。"),
    how = list(
      en = "Pick the <b>sample</b> column, the <b>condition</b> and the two groups to compare. <b>edgeR</b>
            (quasi-likelihood) is the recommended default; <b>limma-voom</b> is similar and fast;
            <b>DESeq2</b> is offered for familiarity. A sample contributes to a cell type only if it has at
            least <b>min cells</b> of it. For paired designs (same patient before / after), add the patient
            as the <b>covariate</b>. You need at least 2 samples per group; 3 or more is far better.",
      zh = "选择<b>样本</b>列、<b>条件</b>列以及要比较的两组。推荐默认使用 <b>edgeR</b>（quasi-likelihood）；<b>limma-voom</b> 结果相近且更快；也提供 <b>DESeq2</b>。一个样本只有在某细胞类型达到<b>最少细胞数</b>时才参与该类型的检验。配对设计（同一患者治疗前后）请把患者设为<b>协变量</b>。每组至少需要 2 个样本，3 个以上好得多。"),
    read = list(
      en = "The overview shows how many genes change in each cell type (FDR within each cell type). The
            volcano shows one cell type: x = log2 fold change (compared group vs reference), y = -log10 p;
            coloured points pass the FDR threshold. A cell type with few samples has little power: no
            significant gene is not evidence of no change.",
      zh = "概览图显示每种细胞类型中有多少基因发生变化（FDR 在每种细胞类型内计算）。火山图显示单一细胞类型：x = log2 倍数变化（比较组相对参照组），y = -log10 p；有颜色的点通过 FDR 阈值。样本少的细胞类型检验力低：没有显著基因不代表没有变化。"),
    example = list(
      en = "Tumour-infiltrating CD8 T cells from 6 responders and 6 non-responders: summed per patient,
            edgeR may find interferon-response genes up in responders, with 12 replicates rather than
            20,000 cells.",
      zh = "6 名应答者和 6 名非应答者的肿瘤浸润 CD8 T 细胞：按患者求和后，edgeR 可能发现应答者中干扰素应答基因上调，重复数是 12 个患者，而不是 20,000 个细胞。")
  )
  tbl_out <- if (has_pkg("DT")) DT::dataTableOutput(ns("tbl")) else shiny::verbatimTextOutput(ns("tbl"))
  smp_out <- if (has_pkg("DT")) DT::dataTableOutput(ns("smp")) else shiny::verbatimTextOutput(ns("smp"))
  controls <- shiny::tagList(
    sample_design_controls(ns, levels = TRUE, covariate = TRUE),
    label_with_help("Method",
                    "edgeR quasi-likelihood F test (recommended), limma-voom (moderated t), or DESeq2 (Wald test; needs the DESeq2 package).",
                    "方法",
                    "edgeR quasi-likelihood F 检验（推荐）、limma-voom（moderated t）或 DESeq2（Wald 检验；需要 DESeq2 包）。"),
    shiny::selectInput(ns("method"), NULL,
                       choices = c("edgeR (quasi-likelihood)" = "edgeR", "limma-voom" = "limma", "DESeq2" = "DESeq2"),
                       selected = "edgeR"),
    label_with_help("Min cells per sample",
                    "A sample enters a cell type's test only if it has at least this many cells of that type; sums of a handful of cells are too noisy.",
                    "每样本最少细胞数",
                    "某样本只有在该细胞类型的细胞数不少于此值时才参与检验；几个细胞的加和噪声太大。"),
    shiny::numericInput(ns("min_cells"), NULL, value = 10, min = 1, step = 1),
    label_with_help("FDR threshold",
                    "Genes with a BH-adjusted p-value below this are called changed (display and counts only; the table has every gene).",
                    "FDR 阈值",
                    "BH 校正 p 值低于此值的基因视为有变化（仅影响显示和计数；表格包含所有基因）。"),
    shiny::numericInput(ns("fdr"), NULL, value = 0.05, min = 0.001, max = 0.5, step = 0.01),
    run_button(ns("run"), "Run pseudobulk DE", "运行 pseudobulk 差异分析")
  )
  step_container(id = id, 
    title     = list(en = "Pseudobulk differential expression", zh = "Pseudobulk 差异表达"),
    subtitle  = list(en = "Compare conditions within each cell type, with samples as replicates.",
                     zh = "在每种细胞类型内比较条件，以样本为重复单位。"),
    explainer = explainer,
    controls  = controls,
    summary   = shiny::uiOutput(ns("summary")),
    preview   = shiny::tagList(
      shiny::uiOutput(ns("insight")),
      bslib::navset_card_tab(
        bslib::nav_panel(i18n("Overview", "概览"),
          preview_plot_ui(ns("preview"), download = TRUE,
            guide = list(en = "The number of changed genes per cell type will be drawn here.",
                         zh = "运行后，这里将绘制每种细胞类型中变化基因的数量。"),
            caption = list(en = "Bars = genes passing the FDR threshold, higher (right) or lower (left) in the compared group; FDR within each cell type.",
                           zh = "条形 = 通过 FDR 阈值的基因数，右侧为比较组中更高，左侧为更低；FDR 在每种细胞类型内计算。"))),
        bslib::nav_panel(i18n("Volcano", "火山图"),
          shiny::uiOutput(ns("volcano_pick")),
          preview_plot_ui(ns("volcano"), download = TRUE,
            caption = list(en = "One point per gene: log2 fold change (compared vs reference) against -log10 p; colour = passes the FDR threshold.",
                           zh = "每个点是一个基因：log2 倍数变化（比较组相对参照组）对 -log10 p；颜色 = 通过 FDR 阈值。"))),
        bslib::nav_panel(i18n("Samples", "样本"),
          shiny::div(class = "omicone-table", smp_out)),
        bslib::nav_panel(i18n("Table", "表格"),
          shiny::div(class = "omicone-table",
                     shiny::downloadButton(ns("csv"), i18n("Download full table (CSV)", "下载完整表格（CSV）"),
                                           class = "btn-sm"),
                     tbl_out))
      )
    )
  )
}

#' @rdname mod_pseudobulk
#' @keywords internal
mod_pseudobulk_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    res <- step_results(rv, "sc", out = NULL, pb_meta = NULL, design = NULL)
    design <- sample_design_server(input, output, session, rv, levels = TRUE, covariate = TRUE)

    shiny::observeEvent(input$run, {
      shiny::req(rv$obj)
      d <- design()
      method <- input$method
      min_cells <- int_input(input$min_cells, 1, 1e6)
      fdr <- num_input(input$fdr, 1e-6, 1)
      if (is.na(min_cells) || is.na(fdr)) {
        shiny::showNotification(i18n("Enter a whole number of cells and an FDR between 0 and 1.",
                                     "请输入整数的细胞数，以及 0 到 1 之间的 FDR。"), type = "error")
        return(NULL)
      }
      md <- obj_meta(rv$obj)
      if (!all(c(d$sample, d$condition, d$group) %in% names(md)) || is.null(d$ref) || is.null(d$alt)) {
        shiny::showNotification(i18n("Choose the sample, condition, groups and cell-type columns first.",
                                     "请先选择样本、条件、比较组和细胞类型列。"), type = "error")
        return(NULL)
      }
      if (identical(d$ref, d$alt)) {
        shiny::showNotification(i18n("The reference and the compared group must differ.",
                                     "参照组和比较组不能相同。"), type = "error")
        return(NULL)
      }
      if (!require_pkgs(c("edgeR", if (method == "limma") "limma", if (method == "DESeq2") "DESeq2"),
                        "Pseudobulk DE")) return(NULL)
      samples <- tryCatch(sample_table(md, d$sample, d$condition, d$covariate),
                          error = function(e) conditionMessage(e))
      if (is.character(samples)) {
        shiny::showNotification(samples, type = "error", duration = 10)
        return(NULL)
      }
      out <- with_progress_notify({
        pb <- pseudobulk_counts(obj_counts(rv$obj), md, d$sample, d$group, min_cells = min_cells)
        r <- pseudobulk_de_all(pb, samples, d$ref, d$alt, method = method,
                               covariate = !is.null(d$covariate), fdr = fdr)
        r$pb_meta <- pb$meta
        r$pb_dropped <- pb$dropped
        r$samples <- samples
        r
      }, message = "Summing counts per sample and testing each cell type...")
      if (is.null(out)) return(NULL)
      if (is.null(out$table)) {
        shiny::showNotification(
          i18n(sprintf("No cell type could be tested: %s", out$skipped$reason[1] %||% "too few samples"),
               sprintf("没有可检验的细胞类型：%s", out$skipped$reason[1] %||% "样本太少")),
          type = "error", duration = 10)
        return(NULL)
      }
      res$out <- out
      res$design <- c(d, list(method = method, min_cells = min_cells, fdr = fdr))
      mark_done(rv, "pseudobulk")
      log_step(log_rv, "Pseudobulk DE",
               params = list(sample = d$sample, condition = d$condition, compare = paste(d$alt, "vs", d$ref),
                             cell_type = d$group, covariate = d$covariate %||% "none", method = method,
                             min_cells = min_cells),
               code = pseudobulk_log_code(d$sample, d$group, d$condition, d$ref, d$alt, method,
                                          min_cells, covariate = d$covariate))
      n_sig <- sum(out$overview$up + out$overview$down)
      shiny::showNotification(
        i18n(sprintf("Tested %d cell type(s); %s gene(s) pass FDR < %g.", nrow(out$overview),
                     format(n_sig, big.mark = ","), fdr),
             sprintf("检验了 %d 种细胞类型；%s 个基因 FDR < %g。", nrow(out$overview),
                     format(n_sig, big.mark = ","), fdr)),
        type = "message")
    })

    output$summary <- shiny::renderUI({
      out <- res$out
      if (is.null(out)) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("Choose the sample and condition columns, then click <b>Run pseudobulk DE</b>.",
                               "选择样本列和条件列后，点击<b>运行 pseudobulk 差异分析</b>。")))
      }
      ov <- out$overview
      shiny::tagList(
        stat_tile(i18n("Comparison", "比较"), sprintf("%s vs %s", out$alt, out$ref)),
        stat_tile(i18n("Samples (ref / compared)", "样本（参照 / 比较）"),
                  sprintf("%d / %d", sum(out$samples$condition == out$ref), sum(out$samples$condition == out$alt))),
        stat_tile(i18n("Cell types tested", "检验的细胞类型"),
                  sprintf("%d%s", nrow(ov), if (nrow(out$skipped)) sprintf(" (+%d skipped)", nrow(out$skipped)) else "")),
        stat_tile(i18n(sprintf("Genes FDR < %g", out$fdr), sprintf("FDR < %g 的基因", out$fdr)),
                  format(sum(ov$up + ov$down), big.mark = ",")),
        stat_tile(i18n("Method", "方法"), out$method)
      )
    })

    output$insight <- shiny::renderUI({
      out <- res$out
      if (is.null(out)) return(NULL)
      ov <- out$overview
      top <- ov[order(-(ov$up + ov$down)), , drop = FALSE][1, ]
      small <- ov$cell_type[pmin(ov$samples_ref, ov$samples_alt) < 3]
      insight_bar(
        sprintf("%s: %s gene(s) change in %s (%d up, %d down, FDR < %g)%s%s. Samples are the replicates; FDR is controlled within each cell type.",
                paste(out$alt, "vs", out$ref), format(top$up + top$down, big.mark = ","), top$cell_type,
                top$up, top$down, out$fdr,
                if (length(small)) sprintf("; fewer than 3 samples per group in %s, so treat those as exploratory",
                                           paste(utils::head(small, 4), collapse = ", ")) else "",
                if (nrow(out$skipped)) sprintf("; %d cell type(s) skipped (see Samples)", nrow(out$skipped)) else ""),
        sprintf("%s：%s 中有 %s 个基因变化（上调 %d、下调 %d，FDR < %g）%s%s。以样本为重复单位；FDR 在每种细胞类型内控制。",
                paste(out$alt, "vs", out$ref), top$cell_type, format(top$up + top$down, big.mark = ","),
                top$up, top$down, out$fdr,
                if (length(small)) sprintf("；%s 每组不足 3 个样本，结果仅供探索", paste(utils::head(small, 4), collapse = "、")) else "",
                if (nrow(out$skipped)) sprintf("；%d 种细胞类型被跳过（见「样本」页）", nrow(out$skipped)) else ""))
    })

    render_step_plot(output, input, "preview", function() {
      out <- res$out
      shiny::req(out, nrow(out$overview) > 0)
      pseudobulk_overview_plot(out$overview, fdr = out$fdr, ref = out$ref, alt = out$alt)
    }, name = "pseudobulk_overview", width = 9, height = 6)

    output$volcano_pick <- shiny::renderUI({
      out <- res$out
      shiny::req(out)
      ov <- out$overview
      ch <- ov$cell_type[order(-(ov$up + ov$down))]
      shiny::div(class = "omicone-inline-pick",
                 shiny::selectInput(session$ns("volcano_ct"), i18n("Cell type", "细胞类型"), choices = ch,
                                    selected = keep_selected(shiny::isolate(input$volcano_ct), ch, ch[1]),
                                    width = "260px"))
    })

    render_step_plot(output, input, "volcano", function() {
      out <- res$out
      shiny::req(out, input$volcano_ct)
      df <- out$table[out$table$cell_type == input$volcano_ct, , drop = FALSE]
      shiny::req(nrow(df) > 0)
      pseudobulk_volcano(df, fdr = out$fdr,
                         title = sprintf("%s: %s vs %s", input$volcano_ct, out$alt, out$ref))
    }, name = "pseudobulk_volcano", width = 8, height = 6.5)

    output$smp <- render_tbl_wrap(function() {
      out <- res$out
      shiny::req(out)
      used <- out$pb_meta
      used$used <- rep(TRUE, nrow(used))
      dropped <- out$pb_dropped
      dropped$used <- rep(FALSE, nrow(dropped))
      m <- rbind(used, dropped)
      m$condition <- out$samples$condition[match(m$sample, out$samples$sample)]
      m <- m[order(m$group, m$condition, m$sample), c("group", "sample", "condition", "n_cells", "used"), drop = FALSE]
      names(m) <- c("cell_type", "sample", "condition", "n_cells", "enough_cells")
      skipped <- out$skipped
      if (nrow(skipped)) {
        m$note <- ifelse(m$cell_type %in% skipped$cell_type,
                         paste("not tested:", skipped$reason[match(m$cell_type, skipped$cell_type)]), "")
      }
      m
    })

    output$tbl <- render_tbl_wrap(function() {
      out <- res$out
      shiny::req(out)
      utils::head(out$table[order(out$table$fdr), , drop = FALSE], 5000)
    })

    output$csv <- shiny::downloadHandler(
      filename = function() sprintf("pseudobulk_%s.csv", format(Sys.time(), "%Y%m%d_%H%M%S")),
      content = function(file) {
        out <- res$out
        utils::write.csv(if (is.null(out)) data.frame() else out$table, file, row.names = FALSE)
      }
    )
  })
}
