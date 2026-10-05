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
            mutation. Empty tile = that sample has no non-synonymous mutation in
            that gene.",
      zh = "基因为行、样本为列，每个格子按突变类型着色。空白格表示该样本在该基因上没有非同义突变。"),
    why  = list(
      en = "It shows at a glance which genes define the cohort, how often they are
            hit, and — from the staircase pattern — which genes tend not to be
            mutated in the same patient.",
      zh = "一眼就能看出哪些基因定义了这个队列、被击中的频率，以及从阶梯状排布中看出哪些基因倾向于不在同一个患者中共同突变。"),
    how  = list(
      en = "Start with the <b>top N</b> genes. Once you have a hypothesis, paste a
            specific <b>gene list</b> instead (matched without regard to case).
            Add <b>clinical annotations</b> to colour the bar above the plot by
            subtype, sex, treatment, and tick <b>sort by annotation</b> to group
            samples by it.",
      zh = "先用 <b>Top N</b> 基因。有了假设之后，改为粘贴具体的<b>基因列表</b>（不区分大小写）。加上<b>临床注释</b>可以按亚型、性别、治疗给顶部的注释条着色，勾选<b>按注释排序</b>可让样本按注释分组。"),
    read = list(
      en = "Each row is a gene, each column a sample; a coloured tile is a
            mutation, its colour the consequence (multi-hit samples are marked
            separately). Bars on top count mutations per sample; bars on the
            right count mutated samples per gene, with the percentage of all
            samples beside. Samples are sorted so co-mutation and mutual
            exclusivity show up as a staircase; the visual pattern is a hint —
            the <b>Drivers &amp; interactions</b> step tests it.",
      zh = "每行一个基因、每列一个样本；有色格子代表突变，颜色标记后果（多重打击的样本单独标记）。顶部柱条是每样本突变数；右侧柱条是每基因的突变样本数，旁边标注占全部样本的百分比。样本经过排序，使共突变与互斥呈现为阶梯状；这种视觉模式只是线索——<b>驱动基因与互作</b>步骤会做统计检验。"),
    example = list(
      en = "Top 20 genes in TCGA LAML: <code>FLT3</code> (27%),
               <code>DNMT3A</code> (25%) and <code>NPM1</code> (17%) lead.
               Annotate by <code>FAB_classification</code> to see how each is
               spread across the FAB subtypes.",
      zh = "TCGA LAML 的 Top 20 基因：<code>FLT3</code>（27%）、<code>DNMT3A</code>（25%）、<code>NPM1</code>（17%）居前。加上 <code>FAB_classification</code> 注释，可以看到它们在各 FAB 亚型中的分布。")
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
  step_container(id = id, 
    title     = list(en = "Oncoplot", zh = "Oncoplot"),
    subtitle  = list(en = "Gene × sample mutation matrix — the waterfall figure.",
                     zh = "基因 × 样本突变矩阵——瀑布图。"),
    explainer = explainer,
    controls  = controls,
    summary   = shiny::uiOutput(ns("summary")),
    preview   = shiny::tagList(
      shiny::uiOutput(ns("insight")),
      preview_plot_ui(ns("plot"), download = TRUE,
                      guide = list(en = "The cohort's oncoplot (waterfall) will be drawn here.",
                                   zh = "运行后，这里将绘制该队列的 Oncoplot（瀑布图）。"),
                      caption = list(en = "Rows = genes, columns = samples with at least one shown gene mutated; tile colour = mutation consequence. Top bars: per-sample burden; right bars and %: mutated samples per gene, of all samples.",
                                     zh = "行＝基因，列＝至少有一个所示基因突变的样本；格子颜色＝突变后果。顶部柱条：每样本突变数；右侧柱条与百分比：每基因的突变样本数（占全部样本）。")))
  )
}

#' @rdname mod_wes_onco
#' @keywords internal
mod_wes_onco_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    ns  <- session$ns
    cfg <- step_result(rv, "wes")

    output$clin_ui <- shiny::renderUI({
      shiny::req(rv$maf)
      cols <- wes_clinical_cols(rv$maf)
      if (!length(cols)) {
        return(shiny::div(class = "omicone-status-empty",
                          i18n("No clinical columns attached to this MAF.",
                               "该 MAF 未附带临床列。")))
      }
      wes_col_select(ns, "clin_feats",
                     label = list(en = "Clinical annotations", zh = "临床注释"),
                     tip = list(en = "Drawn as coloured bars above the oncoplot.",
                                zh = "以彩色条带绘制在 Oncoplot 上方。"),
                     choices = cols, multiple = TRUE)
    })

    draw_with <- function(maf, c0) {
      wes_oncoplot_gg(maf, c0$shown, clin = c0$clin, sort_anno = c0$sort_anno,
                      show_pct = c0$pct, titv = c0$titv)
    }

    shiny::observeEvent(input$run, {
      shiny::req(rv$maf)
      if (!require_pkgs("maftools", "Oncoplot")) return(NULL)
      maf <- rv$maf
      known <- wes_genes(maf)
      genes <- NULL
      if (identical(input$gene_mode, "list")) {
        typed <- parse_genes(input$genes)
        genes <- wes_match_genes(typed, known)
        if (!length(genes)) {
          wes_notify("None of those genes are mutated in this cohort.",
                     "所输入的基因在该队列中均无突变。")
          return(NULL)
        }
        if (length(genes) < length(typed)) {
          wes_notify(sprintf("%d gene(s) not mutated in this cohort were dropped.",
                             length(typed) - length(genes)),
                     sprintf("已去除 %d 个在该队列中无突变的基因。", length(typed) - length(genes)),
                     type = "warning", duration = 8)
        }
      }
      top <- wes_int(input$top, 20, lo = 3, hi = 100)
      # oncoplot(top =) ranks genes by altered samples, as getGeneSummary does
      shown <- if (is.null(genes)) utils::head(known, top) else genes
      c0 <- list(genes = genes, top = top, shown = shown,
                 clin = input$clin_feats,
                 sort_anno = isTRUE(input$sort_anno),
                 titv = isTRUE(input$draw_titv),
                 pct = isTRUE(input$show_pct))
      ok <- with_progress_notify(wes_dry_run(function() draw_plot_object(draw_with(maf, c0))),
                                 message = "Drawing the oncoplot...")
      if (is.null(ok)) return(NULL)
      cfg(c0)
      mark_done(rv, "wes_onco")
      args <- list(maf = quote(maf), top = if (is.null(genes)) top else NULL, genes = genes,
                   clinicalFeatures = if (length(c0$clin)) c0$clin else NULL,
                   sortByAnnotation = if (length(c0$clin) && c0$sort_anno) TRUE else NULL,
                   draw_titv = c0$titv, showPct = c0$pct, removeNonMutated = TRUE)
      log_step(log_rv, "WES oncoplot",
               params = list(genes = if (is.null(genes)) paste("top", top)
                                     else paste(genes, collapse = ", "),
                             annotations = paste(c0$clin, collapse = ", ")),
               code = wes_code("maftools::oncoplot", args))
    })

    stats <- shiny::reactive({
      c0 <- cfg()
      shiny::req(rv$maf, c0)
      fr <- wes_gene_freq(rv$maf, c0$shown)
      list(freq = fr, n_genes = length(c0$shown),
           samples = wes_overview(rv$maf)$samples,
           n_clin = length(c0$clin))
    })

    output$summary <- shiny::renderUI({
      if (is.null(rv$maf)) return(wes_no_maf())
      if (is.null(cfg())) {
        return(wes_prompt("Choose genes and click <b>Draw oncoplot</b>.",
                          "选择基因后点击<b>绘制 Oncoplot</b>。"))
      }
      s <- stats()
      shiny::tagList(
        stat_tile(i18n("Genes shown", "显示基因数"), s$n_genes),
        stat_tile(i18n("Samples", "样本数"), wes_fmt(s$samples)),
        stat_tile(i18n("Annotations", "注释列"), s$n_clin)
      )
    })

    output$insight <- shiny::renderUI({
      if (is.null(rv$maf) || is.null(cfg())) return(NULL)
      s <- stats()
      fr <- s$freq
      if (!nrow(fr)) return(NULL)
      k <- min(3, nrow(fr))
      items <- sprintf("<b>%s</b> %.0f%%", fr$gene[seq_len(k)], fr$pct[seq_len(k)])
      insight_bar(
        sprintf("Most mutated of the %d genes shown: %s of all %s samples. In the matrix, check whether they tile different columns (mutually exclusive) or stack together (co-occurring) — a visual hint only; the interactions test quantifies it.",
                s$n_genes, paste(items, collapse = ", "), wes_fmt(s$samples)),
        sprintf("所示 %d 个基因中突变最多的：%s（占全部 %s 个样本）。在矩阵中看它们是铺满不同列（互斥）还是堆在一起（共现）——这只是视觉线索，互作检验才给出定量结论。",
                s$n_genes, paste(items, collapse = "、"), wes_fmt(s$samples)))
    })

    draw_onco <- function() {
      c0 <- cfg()
      shiny::req(rv$maf, c0)
      draw_with(rv$maf, c0)
    }
    render_step_plot(output, input, "plot", draw_onco, name = "wes_oncoplot",
      width = function() {
        n <- tryCatch(nrow(maftools::getSampleSummary(rv$maf)), error = function(e) 50)
        max(10, min(30, 4 + 0.09 * n))
      },
      height = function() {
        c0 <- cfg()
        n <- if (is.null(c0)) 20 else length(c0$shown)
        max(7, min(20, 3.5 + 0.28 * n + 0.3 * length(c0$clin %||% character(0))))
      })
  })
}
