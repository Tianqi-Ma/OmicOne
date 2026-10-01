#' WES module 1: Import MAF
#'
#' Read a Mutation Annotation Format file (and, optionally, a per-sample
#' clinical table) into a maftools MAF object, which becomes `rv$maf` — the
#' working object for the whole WES pipeline, the way `rv$obj` is for
#' single-cell.
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
            LAML cohort (193 samples) that loads instantly and offline. For your
            own data, a MAF from Mutect2/Strelka/VarScan via <code>vcf2maf</code>
            or the GDC works as-is.",
      zh = "<b>只是想体验？</b>选择<b>演示数据</b>——maftools 自带一份 TCGA LAML 队列（193 个样本），可离线即时加载。用自己的数据时，Mutect2/Strelka/VarScan 经 <code>vcf2maf</code> 转换的 MAF、或 GDC 下载的 MAF 都可直接使用。"),
    read = list(
      en = "The <b>Per-sample</b> tab lists one row per tumour: total variants
            plus a breakdown by consequence — a sample with ten times the median
            is usually a calling artefact. <b>Per-gene</b> flips it around: one
            row per gene, with how many samples and variants hit it.
            <b>Clinical</b> shows the attached table as loaded; check the row
            count there matches the number of samples.",
      zh = "<b>各样本</b>页签每行一个肿瘤样本：变异总数及按后果的分类统计——若某样本是中位数的十倍，多半是检出假阳性。<b>各基因</b>页签则反过来：每行一个基因，列出命中它的样本数与变异数。<b>临床数据</b>页签原样展示载入的临床表；请核对其行数是否等于样本数。"),
    example = list(
      en = "The clinical table needs a <code>Tumor_Sample_Barcode</code> column
               matching the MAF, plus whatever else you have
               (<code>FAB_classification</code>, <code>days_to_last_followup</code>,
               <code>Overall_Survival_Status</code>).",
      zh = "临床表需要一列 <code>Tumor_Sample_Barcode</code> 与 MAF 对应，其余列随意（如 <code>FAB_classification</code>、<code>days_to_last_followup</code>、<code>Overall_Survival_Status</code>）。")
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
                      "CSV/TSV with a Tumor_Sample_Barcode column matching the MAF.",
                      label_zh = "临床表（可选）",
                      tip_zh = "CSV/TSV，需含与 MAF 对应的 Tumor_Sample_Barcode 列。"),
      shiny::fileInput(ns("clin_file"), NULL, accept = c(".csv", ".tsv", ".txt")),
      shiny::selectInput(ns("clin_sep"), i18n("Clinical separator", "临床表分隔符"),
                         c("Tab" = "\t", "Comma" = ","), selected = "\t")
    ),
    run_button(ns("run"), "Load MAF", "加载 MAF")
  )

  step_container(
    title     = list(en = "Import MAF", zh = "导入 MAF"),
    subtitle  = list(en = "Load the mutation file every later step reads from.",
                     zh = "载入突变数据——之后的每一步都从这里读取。"),
    explainer = explainer,
    controls  = controls,
    summary   = shiny::uiOutput(ns("summary")),
    preview   = bslib::navset_card_tab(
      bslib::nav_panel(i18n("Per-sample", "各样本"),   shiny::uiOutput(ns("samp_slot"))),
      bslib::nav_panel(i18n("Per-gene", "各基因"),     shiny::uiOutput(ns("gene_slot"))),
      bslib::nav_panel(i18n("Clinical", "临床数据"),   shiny::uiOutput(ns("clin_slot")))
    )
  )
}

#' @rdname mod_wes_import
#' @keywords internal
mod_wes_import_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    shiny::observeEvent(input$run, {
      if (!require_pkgs("maftools", "WES import")) return(NULL)

      if (identical(input$source, "demo")) {
        p <- wes_demo_paths()
        if (!nzchar(p$maf) || !file.exists(p$maf)) {
          shiny::showNotification(
            "maftools' bundled example was not found in your installation.",
            type = "error", duration = 12)
          return(NULL)
        }
        maf_path <- p$maf
        clin <- if (nzchar(p$clinical) && file.exists(p$clinical)) {
          utils::read.delim(p$clinical, sep = "\t", stringsAsFactors = FALSE)
        } else NULL
        label <- "Demo: TCGA LAML"
        fname <- "tcga_laml.maf.gz"
      } else {
        shiny::req(input$maf_file)
        maf_path <- input$maf_file$datapath
        clin <- if (!is.null(input$clin_file)) {
          tryCatch(utils::read.delim(input$clin_file$datapath, sep = input$clin_sep,
                                     stringsAsFactors = FALSE, check.names = FALSE),
                   error = function(e) {
                     shiny::showNotification(paste("Clinical table:", conditionMessage(e)),
                                             type = "warning", duration = 10)
                     NULL
                   })
        } else NULL
        label <- input$maf_file$name
        fname <- input$maf_file$name
      }

      if (!is.null(clin) && !"Tumor_Sample_Barcode" %in% colnames(clin)) {
        shiny::showNotification(
          "The clinical table has no 'Tumor_Sample_Barcode' column; loading the MAF without it.",
          type = "warning", duration = 12)
        clin <- NULL
      }

      maf <- with_progress_notify(wes_read_maf(maf_path, clinical = clin),
                                  message = "Reading MAF...")
      if (is.null(maf)) return(NULL)

      rv$maf <- maf
      rv$maf_source <- label
      mark_done(rv, "wes_import")
      log_step(log_rv, "WES import",
               params = list(source = input$source, file = fname,
                             clinical = !is.null(clin)),
               code = sprintf('maf <- maftools::read.maf("%s"%s)', fname,
                              if (!is.null(clin)) ", clinicalData = clin" else ""))
      ov <- wes_overview(maf)
      shiny::showNotification(
        sprintf("MAF loaded: %s samples, %s mutated genes.",
                format(ov$samples, big.mark = ","), format(ov$genes, big.mark = ",")),
        type = "message")
    })

    output$summary <- shiny::renderUI({
      maf <- rv$maf
      if (is.null(maf)) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("No MAF yet. Tip: keep <b>Demo data</b> selected and click <b>Load MAF</b> to try it instantly.",
                               "尚无 MAF。提示：保持选中<b>演示数据</b>并点击<b>加载 MAF</b> 即可立即试用。")))
      }
      ov <- wes_overview(maf)
      fmt <- function(x) if (is.null(x) || is.na(x)) "-" else format(round(x), big.mark = ",")
      bslib::layout_columns(
        col_widths = c(3, 3, 3, 3),
        stat_tile(i18n("Samples", "样本数"), fmt(ov$samples)),
        stat_tile(i18n("Mutated genes", "突变基因数"), fmt(ov$genes)),
        stat_tile(i18n("Variants", "变异总数"), fmt(ov$variants)),
        stat_tile(i18n("Median / sample", "中位数/样本"), fmt(ov$median_per_sample))
      )
    })

    tbl_slot <- function(out_id) {
      shiny::renderUI({
        if (is.null(rv$maf)) return(wes_no_maf())
        if (has_pkg("DT")) DT::dataTableOutput(ns(out_id))
        else shiny::verbatimTextOutput(ns(paste0(out_id, "_txt")))
      })
    }
    output$samp_slot <- tbl_slot("samp_tbl")
    output$gene_slot <- tbl_slot("gene_tbl")
    output$clin_slot <- tbl_slot("clin_tbl")

    samp_df <- shiny::reactive({
      shiny::req(rv$maf); as.data.frame(maftools::getSampleSummary(rv$maf))
    })
    gene_df <- shiny::reactive({
      shiny::req(rv$maf); wes_gene_table(rv$maf, n = 200)
    })
    clin_df <- shiny::reactive({
      shiny::req(rv$maf)
      cd <- as.data.frame(maftools::getClinicalData(rv$maf))
      if (!ncol(cd)) data.frame(note = "No clinical data attached.") else cd
    })

    if (has_pkg("DT")) {
      output$samp_tbl <- DT::renderDataTable(
        DT::datatable(samp_df(), rownames = FALSE, filter = "top",
                      options = list(pageLength = 15, scrollX = TRUE)))
      output$gene_tbl <- DT::renderDataTable(
        DT::datatable(gene_df(), rownames = FALSE, filter = "top",
                      options = list(pageLength = 15, scrollX = TRUE)))
      output$clin_tbl <- DT::renderDataTable(
        DT::datatable(clin_df(), rownames = FALSE, filter = "top",
                      options = list(pageLength = 15, scrollX = TRUE)))
    } else {
      output$samp_tbl_txt <- shiny::renderPrint(utils::head(samp_df(), 20))
      output$gene_tbl_txt <- shiny::renderPrint(utils::head(gene_df(), 20))
      output$clin_tbl_txt <- shiny::renderPrint(utils::head(clin_df(), 20))
    }
  })
}
