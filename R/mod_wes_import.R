#' WES module 1: Import MAF
#'
#' Read a Mutation Annotation Format file (and, optionally, a per-sample
#' clinical table) into a maftools MAF object, which becomes `rv$maf` — the
#' working object for the whole WES pipeline, the way `rv$obj` is for
#' single-cell. The clinical table is also kept as read in `rv$wes_clin_raw`:
#' `read.maf()` silently drops clinical rows whose sample has no variant at
#' all, and the survival step needs those patients back as wild-type.
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_wes_import
NULL

#' @rdname mod_wes_import
#' @keywords internal
mod_wes_import_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = list(en = "Import MAF", zh = "导入 MAF"),
    what = list(
      en = "Load somatic variant calls in MAF (Mutation Annotation Format) — one
            row per mutation per sample — plus an optional clinical table.",
      zh = "载入 MAF（Mutation Annotation Format）格式的体细胞变异结果——每行是某个样本的一个突变——以及可选的临床信息表。"),
    why  = list(
      en = "Every later step reads this one object. Attaching clinical data now
            means the cohort comparison, enrichment and survival steps can use it
            without you loading anything again.",
      zh = "之后的每一步都读取这一个对象。现在就接上临床数据，后面的队列比较、富集和生存分析就无需再次加载。"),
    how  = list(
      en = "<b>Just exploring?</b> Choose <b>Demo data</b> — maftools ships a TCGA
            LAML cohort (193 samples, with a 200-patient clinical table) that loads
            instantly and offline. For your own data, a MAF from
            Mutect2/Strelka/VarScan via <code>vcf2maf</code> or the GDC works
            as-is; column names are matched without regard to case. Only the
            classes listed as non-synonymous are counted as mutations;
            <i>Silent</i> (and, unless you tick it, <i>Splice_Region</i>) variants
            are kept aside.",
      zh = "<b>只是想体验？</b>选择<b>演示数据</b>——maftools 自带一份 TCGA LAML 队列（193 个样本，附 200 名患者的临床表），可离线即时加载。用自己的数据时，Mutect2/Strelka/VarScan 经 <code>vcf2maf</code> 转换的 MAF、或 GDC 下载的 MAF 都可直接使用；列名不区分大小写。只有非同义类别才计为突变；<i>Silent</i>（以及未勾选时的 <i>Splice_Region</i>）变异另行保存，不计入。"),
    read = list(
      en = "The <b>Per-sample</b> tab lists one row per tumour: non-synonymous
            variants plus a breakdown by consequence. A sample far above the
            median can be genuinely hypermutated (POLE, mismatch-repair
            deficiency) or carry unfiltered artefacts — check its signature
            before excluding it. <b>Per-gene</b> flips it around: one row per
            gene, with how many samples and variants hit it. <b>Clinical</b>
            shows the attached table as read, with a column saying whether each
            sample has a MAF record.",
      zh = "<b>各样本</b>页签每行一个肿瘤样本：非同义变异数及按后果的分类统计。远高于中位数的样本可能是真实的超突变（POLE、错配修复缺陷），也可能含未过滤的假阳性——排除前请先看它的突变特征。<b>各基因</b>页签则反过来：每行一个基因，列出命中它的样本数与变异数。<b>临床数据</b>页签原样展示读入的临床表，并有一列标明每个样本是否有 MAF 记录。"),
    example = list(
      en = "The clinical table needs a sample-id column matching the MAF's
               <code>Tumor_Sample_Barcode</code> (pick it in the dropdown), plus
               whatever else you have (<code>FAB_classification</code>,
               <code>days_to_last_followup</code>,
               <code>Overall_Survival_Status</code>). In the demo, 7 of the 200
               clinical patients have no MAF record.",
      zh = "临床表需要一列与 MAF 的 <code>Tumor_Sample_Barcode</code> 对应的样本编号（在下拉框中选择），其余列随意（如 <code>FAB_classification</code>、<code>days_to_last_followup</code>、<code>Overall_Survival_Status</code>）。演示数据的 200 名临床患者中有 7 名没有 MAF 记录。")
  )

  controls <- shiny::tagList(
    label_with_help("Data source",
                    "The demo is maftools' bundled TCGA LAML cohort: instant, offline, no download.",
                    label_zh = "数据来源",
                    tip_zh = "演示数据是 maftools 自带的 TCGA LAML 队列：即时、离线、无需下载。"),
    shiny::radioButtons(ns("source"), NULL,
                        c("Demo data (TCGA LAML)" = "demo", "Upload files" = "upload"),
                        selected = "demo"),
    shiny::conditionalPanel(
      sprintf("input['%s'] == 'upload'", ns("source")),
      label_with_help("MAF file", "A .maf or .maf.gz file (or any tab-delimited MAF-shaped table).",
                      label_zh = "MAF 文件", tip_zh = ".maf 或 .maf.gz 文件（或任何 MAF 结构的制表符分隔表格）。"),
      shiny::fileInput(ns("maf_file"), NULL, accept = c(".maf", ".gz", ".txt", ".tsv")),
      label_with_help("Clinical table (optional)",
                      "CSV (comma) or TSV/TXT (tab), chosen from the file extension. Sample ids are read as text, so 001 stays 001.",
                      label_zh = "临床表（可选）",
                      tip_zh = "CSV（逗号）或 TSV/TXT（制表符），按扩展名自动判断。样本编号按文本读取，001 不会变成 1。"),
      shiny::fileInput(ns("clin_file"), NULL, accept = c(".csv", ".tsv", ".txt")),
      shiny::uiOutput(ns("clin_id_ui"))
    ),
    label_with_help("Reading options",
                    "Splice_Region: variants within a few bases of a splice site; maftools counts them as non-synonymous only when asked. TCGA barcodes: keep the first 12 characters (patient level). FLAGS: 20 long genes (TTN, MUC16 ...) mutated by chance in most exomes.",
                    label_zh = "读取选项",
                    tip_zh = "Splice_Region：距剪接位点几个碱基以内的变异；只有勾选时 maftools 才把它计为非同义。TCGA 条码：只保留前 12 个字符（患者层面）。FLAGS：20 个因基因长而在多数外显子组中偶然突变的基因（TTN、MUC16 等）。"),
    shiny::checkboxInput(ns("splice_region"),
                         i18n("Count Splice_Region as non-synonymous", "把 Splice_Region 计为非同义"),
                         value = FALSE),
    shiny::checkboxInput(ns("is_tcga"),
                         i18n("TCGA barcodes: keep 12 characters", "TCGA 条码：保留前 12 位"),
                         value = FALSE),
    shiny::checkboxInput(ns("rm_flags"),
                         i18n("Remove the 20 FLAGS genes", "移除 20 个 FLAGS 基因"),
                         value = FALSE),
    run_button(ns("run"), "Load MAF", "加载 MAF")
  )

  step_container(
    title     = list(en = "Import MAF", zh = "导入 MAF"),
    subtitle  = list(en = "Load the mutation file every later step reads from.",
                     zh = "载入突变数据——之后的每一步都从这里读取。"),
    explainer = explainer,
    controls  = controls,
    summary   = shiny::uiOutput(ns("summary")),
    preview   = shiny::tagList(
      shiny::uiOutput(ns("insight")),
      bslib::navset_card_tab(
        bslib::nav_panel(i18n("Per-sample", "各样本"), shiny::uiOutput(ns("samp_slot"))),
        bslib::nav_panel(i18n("Per-gene", "各基因"),   shiny::uiOutput(ns("gene_slot"))),
        bslib::nav_panel(i18n("Clinical", "临床数据"), shiny::uiOutput(ns("clin_slot")))
      ))
  )
}

#' @rdname mod_wes_import
#' @keywords internal
mod_wes_import_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # The uploaded clinical table, read once so its id column can be chosen
    clin_upload <- shiny::reactive({
      f <- input$clin_file
      shiny::req(f)
      tryCatch(wes_read_clinical(f$datapath, f$name),
               error = function(e) {
                 wes_notify(paste("Clinical table:", conditionMessage(e)),
                            paste("临床表：", conditionMessage(e)), type = "warning")
                 NULL
               })
    })

    output$clin_id_ui <- shiny::renderUI({
      df <- clin_upload()
      shiny::req(df)
      cols <- names(df)
      guess <- wes_find_col(cols, "Tumor_Sample_Barcode")
      if (is.null(guess)) {
        like <- grep("barcode|sample|^id$|_id$", cols, ignore.case = TRUE, value = TRUE)
        guess <- if (length(like)) like[1] else cols[1]
      }
      wes_col_select(ns, "clin_id",
                     label = list(en = "Sample-id column", zh = "样本编号列"),
                     tip = list(en = "The clinical column whose values match the MAF's Tumor_Sample_Barcode.",
                                zh = "临床表中与 MAF 的 Tumor_Sample_Barcode 取值对应的列。"),
                     choices = cols,
                     selected = keep_selected(shiny::isolate(input$clin_id), cols, guess))
    })

    shiny::observeEvent(input$run, {
      if (!require_pkgs("maftools", "WES import")) return(NULL)
      src <- input$source %||% "demo"
      splice <- isTRUE(input$splice_region)
      is_tcga <- isTRUE(input$is_tcga)
      rm_flags <- isTRUE(input$rm_flags)

      clin <- NULL
      id_from <- NULL
      clin_name <- NULL
      if (identical(src, "demo")) {
        p <- wes_demo_paths()
        if (!nzchar(p$maf) || !file.exists(p$maf)) {
          wes_notify("maftools' bundled example was not found in your installation.",
                     "未在已安装的 maftools 中找到自带示例数据。")
          return(NULL)
        }
        maf_path <- p$maf
        if (nzchar(p$clinical) && file.exists(p$clinical)) {
          clin <- wes_read_clinical(p$clinical, p$clinical)
          clin_name <- basename(p$clinical)
          id_from <- "Tumor_Sample_Barcode"
        }
        label <- "Demo: TCGA LAML"
        fname <- "tcga_laml.maf.gz"
      } else {
        if (is.null(input$maf_file)) {
          wes_notify("Choose a MAF file first.", "请先选择 MAF 文件。", type = "warning")
          return(NULL)
        }
        maf_path <- input$maf_file$datapath
        label <- input$maf_file$name
        fname <- input$maf_file$name
        if (!is.null(input$clin_file)) {
          id_from <- input$clin_id
          clin <- tryCatch(wes_read_clinical(input$clin_file$datapath, input$clin_file$name,
                                             id_col = id_from),
                           error = function(e) NULL)
          clin_name <- input$clin_file$name
          if (is.null(clin) || is.null(id_from) || !id_from %in% names(clin)) {
            wes_notify("Pick the sample-id column of the clinical table; loading the MAF without it.",
                       "请选择临床表的样本编号列；本次仅加载 MAF。", type = "warning", duration = 12)
            clin <- NULL
            clin_name <- NULL
          }
        }
      }
      if (!is.null(clin)) {
        text_cols <- attr(clin, "text_cols")
        if (!identical(id_from, "Tumor_Sample_Barcode")) {
          names(clin)[names(clin) == "Tumor_Sample_Barcode"] <- "Tumor_Sample_Barcode_orig"
          names(clin)[names(clin) == id_from] <- "Tumor_Sample_Barcode"
        }
        attr(clin, "text_cols") <- text_cols
      }

      hdr <- tryCatch(wes_maf_header(maf_path), error = function(e) list(skip = 0L,
                                                                         renamed = character(0)))
      # clinical ids that differ from the MAF only in case / spaces take the
      # MAF's spelling for read.maf(); rv$wes_clin_raw keeps them as read
      clin_for_maf <- clin
      aligned <- 0L
      if (!is.null(clin)) {
        tsb <- tryCatch(wes_maf_barcodes(maf_path, hdr, is_tcga), error = function(e) NULL)
        if (!is.null(tsb)) {
          al <- wes_align_ids(clin$Tumor_Sample_Barcode, tsb)
          clin_for_maf$Tumor_Sample_Barcode <- al$ids
          aligned <- al$n
        }
      }
      vc <- if (splice) wes_vc_nonsyn(TRUE) else NULL
      maf <- with_progress_notify(
        wes_read_maf(maf_path, clinical = clin_for_maf, vc_nonSyn = vc, is_tcga = is_tcga,
                     rm_flags = rm_flags),
        message = "Reading MAF...")
      if (is.null(maf)) return(NULL)

      # a new cohort: forget every WES result computed on the previous one
      start_epoch(rv, "wes", log_rv)
      rv$maf <- maf
      rv$maf_source <- label
      if (!is.null(clin)) {
        attr(clin, "file") <- clin_name
        attr(clin, "aligned") <- aligned
      }
      rv$wes_clin_raw <- clin
      mark_done(rv, "wes_import")
      read_args <- list(isTCGA = is_tcga, rmFlags = rm_flags, vc_nonSyn = vc)
      log_step(log_rv, "WES import",
               params = list(source = src, file = fname, clinical = clin_name %||% "(none)",
                             clinical_id = id_from %||% "(none)",
                             splice_region_nonsyn = splice, isTCGA = is_tcga,
                             rmFlags = rm_flags),
               code = wes_import_code(maf_file = fname, clin_file = clin_name,
                                      demo = identical(src, "demo"), id_from = id_from,
                                      text_cols = attr(clin, "text_cols"),
                                      renamed = hdr$renamed, skip = hdr$skip,
                                      sep = if (!is.null(clin_name)) wes_sep_for(clin_name) else "\t",
                                      read_args = read_args, aligned = aligned))
      ov <- wes_overview(maf)
      wes_notify(sprintf("MAF loaded: %s samples, %s mutated genes.",
                         wes_fmt(ov$samples), wes_fmt(ov$genes)),
                 sprintf("MAF 已加载：%s 个样本，%s 个突变基因。",
                         wes_fmt(ov$samples), wes_fmt(ov$genes)),
                 type = "message", duration = 5)
    })

    # every number shown in the pills and the insight, computed once
    stats <- shiny::reactive({
      maf <- rv$maf
      shiny::req(maf)
      ov <- wes_overview(maf)
      clin <- rv$wes_clin_raw
      cl <- NULL
      if (!is.null(clin) && "Tumor_Sample_Barcode" %in% names(clin)) {
        ids <- as.character(clin$Tumor_Sample_Barcode)
        maf_ids <- wes_samples(maf)
        cl <- list(rows = nrow(clin),
                   matched = sum(wes_norm_id(ids) %in% wes_norm_id(maf_ids)),
                   aligned = attr(clin, "aligned") %||% 0L,
                   no_maf = sum(!wes_norm_id(ids) %in% wes_norm_id(maf_ids)))
      }
      list(ov = ov, clin = cl)
    })

    output$summary <- shiny::renderUI({
      if (is.null(rv$maf)) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("No MAF yet. Tip: keep <b>Demo data</b> selected and click <b>Load MAF</b> to try it instantly.",
                               "尚无 MAF。提示：保持选中<b>演示数据</b>并点击<b>加载 MAF</b> 即可立即试用。")))
      }
      s <- stats()
      shiny::tagList(
        stat_tile(i18n("Samples", "样本数"), wes_fmt(s$ov$samples)),
        stat_tile(i18n("Mutated genes", "突变基因数"), wes_fmt(s$ov$genes)),
        stat_tile(i18n("Non-synonymous variants", "非同义变异数"), wes_fmt(s$ov$variants)),
        stat_tile(i18n("Median / sample", "中位数/样本"), wes_fmt(s$ov$median_per_sample)),
        if (!is.null(s$clin)) {
          stat_tile(i18n("Clinical rows in MAF", "临床行有 MAF 记录"),
                    sprintf("%s / %s", wes_fmt(s$clin$matched), wes_fmt(s$clin$rows)))
        }
      )
    })

    output$insight <- shiny::renderUI({
      if (is.null(rv$maf)) {
        return(explain_scene("wes_import",
                             "Upload a MAF (or load the demo cohort) on the left to begin.",
                             "在左侧上传 MAF 文件（或加载演示队列）开始。"))
      }
      s <- stats()
      ov <- s$ov
      top <- if (is.na(ov$top_gene)) "-" else sprintf("%s (%.0f%%)", ov$top_gene, ov$top_pct)
      clin_en <- ""
      clin_zh <- ""
      if (!is.null(s$clin)) {
        clin_en <- sprintf(" Clinical table: %s rows, %s match a MAF sample; %s have no MAF record (read.maf keeps only samples with variants — the survival step can count them as sequenced wild-type).",
                           wes_fmt(s$clin$rows), wes_fmt(s$clin$matched), wes_fmt(s$clin$no_maf))
        clin_zh <- sprintf("临床表：%s 行，其中 %s 行与 MAF 样本匹配；%s 行没有 MAF 记录（read.maf 只保留有变异的样本——生存分析步骤可把它们作为已测序的野生型计入）。",
                           wes_fmt(s$clin$rows), wes_fmt(s$clin$matched), wes_fmt(s$clin$no_maf))
        if (s$clin$aligned > 0) {
          clin_en <- paste0(clin_en, sprintf(" %s of these ids differ from the MAF only in case or spaces and were aligned to the MAF's spelling.",
                                             wes_fmt(s$clin$aligned)))
          clin_zh <- paste0(clin_zh, sprintf("其中 %s 个编号与 MAF 仅在大小写或空格上不同，已按 MAF 的写法对齐。",
                                             wes_fmt(s$clin$aligned)))
        }
      }
      insight_bar(
        sprintf("Loaded <b>%s</b> samples with <b>%s</b> non-synonymous variants across <b>%s</b> genes — median %s per sample (silent variants are kept aside, not counted). Top gene: <b>%s</b>.%s",
                wes_fmt(ov$samples), wes_fmt(ov$variants), wes_fmt(ov$genes),
                wes_fmt(ov$median_per_sample), top, clin_en),
        sprintf("已载入 <b>%s</b> 个样本、<b>%s</b> 个非同义变异、<b>%s</b> 个基因——每样本中位 %s 个（同义变异另存，不计入）。最高频基因：<b>%s</b>。%s",
                wes_fmt(ov$samples), wes_fmt(ov$variants), wes_fmt(ov$genes),
                wes_fmt(ov$median_per_sample), top, clin_zh))
    })

    has_maf <- function() !is.null(rv$maf)
    samp_df <- shiny::reactive({
      shiny::req(rv$maf)
      as.data.frame(maftools::getSampleSummary(rv$maf))
    })
    gene_df <- shiny::reactive({
      shiny::req(rv$maf)
      as.data.frame(maftools::getGeneSummary(rv$maf))
    })
    clin_df <- shiny::reactive({
      shiny::req(rv$maf, rv$wes_clin_raw)
      d <- rv$wes_clin_raw
      d$in_MAF <- wes_norm_id(d$Tumor_Sample_Barcode) %in% wes_norm_id(wes_samples(rv$maf))
      d
    })

    t_samp <- wes_table(ns, "samp_tbl", samp_df, "wes_sample_summary", has_maf,
                        function() TRUE)
    output$samp_slot <- t_samp$slot
    output$samp_tbl <- t_samp$table
    output$samp_tbl_dl <- t_samp$download

    t_gene <- wes_table(ns, "gene_tbl", gene_df, "wes_gene_summary", has_maf,
                        function() TRUE)
    output$gene_slot <- t_gene$slot
    output$gene_tbl <- t_gene$table
    output$gene_tbl_dl <- t_gene$download

    t_clin <- wes_table(ns, "clin_tbl", clin_df, "wes_clinical", has_maf,
                        function() !is.null(rv$wes_clin_raw),
                        not_ready = list(en = "No clinical table was attached at import.",
                                         zh = "导入时未附带临床表。"))
    output$clin_slot <- t_clin$slot
    output$clin_tbl <- t_clin$table
    output$clin_tbl_dl <- t_clin$download
  })
}
