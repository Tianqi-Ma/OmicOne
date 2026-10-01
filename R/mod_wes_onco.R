#' WES module 3: Oncoplot
#'
#' The waterfall / oncoprint view: genes as rows, samples as columns, coloured by
#' mutation type, sorted so the mutually exclusive structure of the cohort shows
#' up as a staircase. The single most-used figure in a mutation paper.
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_wes_onco
NULL

#' @rdname mod_wes_onco
#' @keywords internal
mod_wes_onco_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = list(en = "Oncoplot", zh = "Oncoplot"),
    what = list(
      en = "Genes as rows, samples as columns, each tile coloured by the kind of
            mutation. Empty tile = that sample has no mutation in that gene.",
      zh = "基因为行、样本为列，每个格子按突变类型着色。空白格表示该样本在该基因上没有突变。"),
    why  = list(
      en = "It shows at a glance which genes define the cohort, how often they are
            hit, and — from the staircase pattern — which genes tend not to be
            mutated in the same patient.",
      zh = "一眼就能看出哪些基因定义了这个队列、被击中的频率，以及从阶梯状排布中看出哪些基因倾向于不在同一个患者中共同突变。"),
    how  = list(
      en = "Start with the <b>top N</b> genes. Once you have a hypothesis, paste a
            specific <b>gene list</b> instead. Add <b>clinical annotations</b> to
            colour the bar above the plot by subtype, sex, treatment, and tick
            <b>sort by annotation</b> to group samples by it.",
      zh = "先用 <b>Top N</b> 基因。有了假设之后，改为粘贴具体的<b>基因列表</b>。加上<b>临床注释</b>可以按亚型、性别、治疗给顶部的注释条着色，勾选<b>按注释排序</b>可让样本按注释分组。"),
    example = list(
      en = "Top 20 genes in TCGA LAML, annotated by
               <code>FAB_classification</code>: the M3 samples separate out on
               <code>PML</code>/<code>RARA</code> fusion background.",
      zh = "TCGA LAML 的 Top 20 基因，按 <code>FAB_classification</code> 注释：M3 样本会在 <code>PML</code>/<code>RARA</code> 融合背景上被区分出来。")
  )
  controls <- shiny::tagList(
    label_with_help("Gene selection", "Top N by mutation frequency, or a list you type.",
                    label_zh = "基因选择", tip_zh = "按突变频率取 Top N，或自行输入基因列表。"),
    shiny::radioButtons(ns("gene_mode"), NULL,
                        c("Top N genes" = "top", "My gene list" = "list"),
                        selected = "top", inline = TRUE),
    shiny::conditionalPanel(
      sprintf("input['%s'] == 'top'", ns("gene_mode")),
      shiny::numericInput(ns("top"), i18n("Top N", "Top N"), value = 20, min = 3,
                          max = 100, step = 1)),
    shiny::conditionalPanel(
      sprintf("input['%s'] == 'list'", ns("gene_mode")),
      shiny::textAreaInput(ns("genes"), i18n("Genes (comma separated)", "基因（逗号分隔）"),
                           placeholder = "TP53, FLT3, DNMT3A, NPM1", rows = 3)),
    shiny::uiOutput(ns("clin_ui")),
    shiny::checkboxInput(ns("sort_anno"),
                         i18n("Sort samples by annotation", "按注释排序样本"), value = FALSE),
    shiny::checkboxInput(ns("draw_titv"),
                         i18n("Add TiTv panel", "附加 TiTv 面板"), value = FALSE),
    shiny::checkboxInput(ns("show_pct"),
                         i18n("Show mutation percentages", "显示突变百分比"), value = TRUE),
    run_button(ns("run"), "Draw oncoplot", "绘制 Oncoplot")
  )
  step_container(
    title     = list(en = "Oncoplot", zh = "Oncoplot"),
    explainer = explainer,
    controls  = controls,
    summary   = shiny::uiOutput(ns("summary")),
    preview   = preview_plot_ui(ns("plot"), download = TRUE)
  )
}

#' @rdname mod_wes_onco
#' @keywords internal
mod_wes_onco_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    ns  <- session$ns
    cfg <- shiny::reactiveVal(NULL)

    output$clin_ui <- shiny::renderUI({
      shiny::req(rv$maf)
      cols <- wes_clinical_cols(rv$maf)
      if (!length(cols)) {
        return(shiny::div(class = "omicone-status-empty",
                          i18n("No clinical columns attached to this MAF.",
                               "该 MAF 未附带临床列。")))
      }
      shiny::tagList(
        label_with_help("Clinical annotations",
                        "Drawn as coloured bars above the oncoplot.",
                        label_zh = "临床注释", tip_zh = "以彩色条带绘制在 Oncoplot 上方。"),
        shiny::selectizeInput(ns("clin_feats"), NULL, choices = cols, multiple = TRUE)
      )
    })

    parse_genes <- function(txt) {
      if (is.null(txt) || !nzchar(trimws(txt))) return(character(0))
      g <- trimws(unlist(strsplit(txt, "[,\n\t ]+")))
      g[nzchar(g)]
    }

    shiny::observeEvent(input$run, {
      shiny::req(rv$maf)
      if (!require_pkgs("maftools", "Oncoplot")) return(NULL)
      genes <- if (identical(input$gene_mode, "list")) parse_genes(input$genes) else NULL
      if (identical(input$gene_mode, "list")) {
        known <- intersect(genes, wes_genes(rv$maf))
        if (!length(known)) {
          shiny::showNotification("None of those genes are mutated in this cohort.",
                                  type = "error", duration = 10)
          return(NULL)
        }
        if (length(known) < length(genes)) {
          shiny::showNotification(
            sprintf("%d gene(s) not mutated in this cohort were dropped.",
                    length(genes) - length(known)),
            type = "warning", duration = 8)
        }
        genes <- known
      }
      cfg(list(genes = genes, top = max(3, as.integer(input$top %||% 20)),
               clin = input$clin_feats,
               sort_anno = isTRUE(input$sort_anno),
               titv = isTRUE(input$draw_titv),
               pct = isTRUE(input$show_pct)))
      mark_done(rv, "wes_onco")
      log_step(log_rv, "WES oncoplot",
               params = list(genes = if (is.null(genes)) paste("top", input$top)
                                     else paste(genes, collapse = ", "),
                             annotations = paste(input$clin_feats, collapse = ", ")),
               code = if (is.null(genes))
                 sprintf('maftools::oncoplot(maf, top = %s)', input$top)
               else
                 sprintf('maftools::oncoplot(maf, genes = c("%s"))',
                         paste(genes, collapse = '", "')))
    })

    output$summary <- shiny::renderUI({
      if (is.null(rv$maf)) return(wes_no_maf())
      c0 <- cfg()
      if (is.null(c0)) {
        return(wes_prompt("Choose genes and click <b>Draw oncoplot</b>.",
                          "选择基因后点击<b>绘制 Oncoplot</b>。"))
      }
      ov <- wes_overview(rv$maf)
      n_genes <- if (is.null(c0$genes)) c0$top else length(c0$genes)
      bslib::layout_columns(
        col_widths = c(4, 4, 4),
        stat_tile(i18n("Genes shown", "显示基因数"), n_genes),
        stat_tile(i18n("Samples", "样本数"), format(ov$samples, big.mark = ",")),
        stat_tile(i18n("Annotations", "注释列"),
                  if (length(c0$clin)) length(c0$clin) else 0)
      )
    })

    draw_onco <- function() {
      c0 <- cfg(); shiny::req(rv$maf, c0)
      n_genes <- if (is.null(c0$genes)) c0$top else length(c0$genes)
      args <- list(maf = rv$maf, draw_titv = c0$titv, showTumorSampleBarcodes = FALSE,
                   removeNonMutated = TRUE,
                   fontSize = adaptive_cex(n_genes, base = 1.0, n_ref = 20,
                                           lo = 0.7, hi = 1.15))
      if (is.null(c0$genes)) args$top <- c0$top else args$genes <- c0$genes
      if (length(c0$clin)) {
        args$clinicalFeatures <- c0$clin
        if (c0$sort_anno) args$sortByAnnotation <- TRUE
      }
      if (!c0$pct) args$showPct <- FALSE
      do.call(maftools::oncoplot, args)
    }
    output$plot <- render_base_plot(draw_onco)
    register_figure_download(
      output, input, "plot", draw_onco, "wes_oncoplot",
      width = function() {
        n <- tryCatch(nrow(maftools::getSampleSummary(rv$maf)), error = function(e) 50)
        max(10, min(30, 4 + 0.09 * n))
      },
      height = function() {
        c0 <- cfg()
        n <- if (is.null(c0)) 20 else if (is.null(c0$genes)) c0$top else length(c0$genes)
        max(7, min(20, 3 + 0.28 * n))
      })
  })
}
