#' Module: Enrichment & GSEA
#'
#' Turn per-group differential expression into interpretable biology: which
#' annotated gene sets (GO, KEGG, Reactome, WikiPathways, MSigDB Hallmark) are
#' over-represented among a group's up-regulated genes (ORA) or shifted along
#' its full ranked gene list (GSEA). The DE table both read comes from
#' scop::RunDEtest(); see sc_enrichment() / sc_gsea() in fct_scop.R.
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_enrichment
NULL

#' @rdname mod_enrichment
#' @keywords internal
mod_enrichment_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = list(en = "Enrichment & GSEA", zh = "富集与 GSEA"),
    what = list(
      en = "Summarise each group's differentially expressed genes into the pathways
            and functions they represent.",
      zh = "把每个分组的差异表达基因归纳为它们所代表的通路与功能。"),
    why  = list(
      en = "A list of gene names is hard to interpret; enrichment says, for example,
            that a cluster is enriched for <i>antigen presentation</i>. Without a
            proper background, ORA flags almost any expressed pathway as enriched.",
      zh = "一串基因名很难解读；富集分析能告诉你某个簇富集于例如<i>抗原呈递</i>。如果背景基因集不合适，ORA 几乎会把任何表达的通路都判为富集。"),
    how  = list(
      en = "The step first runs a one-vs-rest Wilcoxon test per group
            (scop::RunDEtest, both directions, BH). <b>ORA</b> takes each group's
            up-regulated genes (log2FC above the cut-off, q &lt; 0.05) and tests
            them against the genes that were tested, not the whole database.
            <b>GSEA</b> ranks every tested gene by log2FC.",
      zh = "这一步先对每个分组做“一对其余”的 Wilcoxon 检验（scop::RunDEtest，双向，BH 校正）。<b>ORA</b> 取每组上调基因（log2FC 高于阈值且 q &lt; 0.05），以“参与检验的基因”而非整个数据库为背景进行检验。<b>GSEA</b> 按 log2FC 对全部参与检验的基因排序。"),
    read = list(
      en = "Rows are gene sets, columns are groups; only sets with BH q &lt; 0.05
            are drawn. For ORA dot size is the gene ratio and colour the q-value;
            for GSEA colour is the normalised enrichment score (NES).",
      zh = "行为基因集，列为分组；只绘制 BH q &lt; 0.05 的基因集。ORA 中点大小为基因比例、颜色为 q 值；GSEA 中颜色为标准化富集分数（NES）。"),
    example = list(
      en = "For a cytotoxic T-cell cluster, GO terms such as
            <code>T cell mediated cytotoxicity</code> should rise to the top; a hit
            resting on two genes is fragile.",
      zh = "对于细胞毒性 T 细胞簇，<code>T cell mediated cytotoxicity</code> 之类的 GO 条目应排在前列；仅由两个基因支撑的结果不可靠。")
  )
  controls <- shiny::tagList(
    label_with_help("Analysis type",
                    "ORA tests each group's up-regulated genes; GSEA uses every tested gene ranked by log2FC.",
                    label_zh = "分析类型",
                    tip_zh = "ORA 检验每组的上调基因；GSEA 使用按 log2FC 排序的全部参与检验的基因。"),
    shiny::radioButtons(ns("analysis"), NULL,
                        c("Over-representation (ORA)" = "ora", "GSEA" = "gsea"),
                        selected = "ora"),
    label_with_help("Group by (metadata column)",
                    "Groups compared one-vs-rest in the DE test (e.g. seurat_clusters, celltype).",
                    label_zh = "分组依据（元数据列）",
                    tip_zh = "差异检验中“一对其余”比较的分组（例如 seurat_clusters、celltype）。"),
    shiny::selectInput(ns("group_by"), NULL, choices = NULL),
    label_with_help("Database",
                    "Gene-set collection, as named by scop::PrepareDB().",
                    label_zh = "数据库",
                    tip_zh = "基因集集合，名称与 scop::PrepareDB() 一致。"),
    shiny::selectInput(ns("db"), NULL, choices = enrich_db_choices(), selected = "GO_BP"),
    label_with_help("Species",
                    "Must match your data; guessed from gene-symbol casing on import.",
                    label_zh = "物种",
                    tip_zh = "必须与数据一致；导入时根据基因符号大小写自动推断。"),
    shiny::selectInput(ns("species"), NULL,
                       choices = c("Homo sapiens" = "Homo_sapiens",
                                   "Mus musculus" = "Mus_musculus"),
                       selected = "Homo_sapiens"),
    shiny::conditionalPanel(
      sprintf("input['%s'] == 'ora'", ns("analysis")),
      label_with_help("Minimum log2 fold-change (ORA)",
                      "A gene enters a group's list when avg log2FC exceeds this and BH q < 0.05.",
                      label_zh = "最小 log2 倍数变化（ORA）",
                      tip_zh = "当 avg log2FC 超过该值且 BH q < 0.05 时，基因进入该组的基因列表。"),
      shiny::numericInput(ns("lfc"), NULL, value = 0.25, min = 0, max = 5, step = 0.05)),
    run_button(ns("run"), "Run enrichment", "运行富集分析")
  )
  step_container(id = id, 
    title     = list(en = "Enrichment & GSEA", zh = "富集与 GSEA"),
    subtitle  = list(en = "Pathways and gene sets behind each group's DE genes.",
                     zh = "每个分组差异基因背后的通路与基因集。"),
    explainer = explainer,
    controls  = controls,
    summary   = shiny::uiOutput(ns("summary")),
    preview   = shiny::tagList(
      shiny::uiOutput(ns("insight")),
      preview_plot_ui(ns("preview"), download = TRUE,
        guide = list(en = "Significant gene sets per group will be drawn here.",
                     zh = "运行后，这里将绘制每个分组的显著基因集。"),
        caption = list(
          en = "Top gene sets (rows, BH q < 0.05) per group (columns). ORA: size = gene ratio, colour = q. GSEA: colour = NES, size = set size.",
          zh = "每个分组（列）的头部基因集（行，BH q < 0.05）。ORA：点大小＝基因比例，颜色＝q 值。GSEA：颜色＝NES，点大小＝基因集大小。")))
  )
}

#' @rdname mod_enrichment
#' @keywords internal
mod_enrichment_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    res <- step_results(rv, "sc", done = FALSE, analysis = NULL, group_by = NULL,
                        db = NULL, species = NULL, lfc = NA, summary = NULL,
                        universe = NA)

    shiny::observe({
      obj <- rv$obj
      cols <- categorical_cols(obj_meta(obj))
      def <- default_group_col(cols, obj_misc(obj, "omicone_cluster_col"))
      shiny::updateSelectInput(session, "group_by", choices = cols,
                               selected = keep_selected(shiny::isolate(input$group_by),
                                                        cols, def))
    })

    # A new dataset: start from the species its gene symbols suggest.
    shiny::observeEvent(rv$epoch_sc, {
      shiny::req(rv$obj)
      shiny::updateSelectInput(session, "species", selected = scop_species(rv$obj))
    })

    shiny::observeEvent(input$run, {
      shiny::req(rv$obj)
      if (!require_pkgs("scop", "Enrichment / GSEA")) return(NULL)
      analysis <- input$analysis
      group_by <- input$group_by
      db <- input$db
      species <- input$species
      lfc <- num_input(input$lfc, lo = 0, hi = 5)
      shiny::req(group_by, db, species, !is.na(lfc))
      # The stored DE table is reused only while nothing upstream has re-run.
      reuse <- identical(step_state(rv, "enrichment"), "done")
      obj <- with_progress_notify({
        if (identical(analysis, "ora")) {
          sc_enrichment(rv$obj, group_by = group_by, db = db, species = species,
                        lfc = lfc, reuse_de = reuse)
        } else {
          sc_gsea(rv$obj, group_by = group_by, db = db, species = species,
                  reuse_de = reuse)
        }
      }, message = if (identical(analysis, "ora"))
        "Running DE test and over-representation analysis..." else "Running DE test and GSEA...")
      if (is.null(obj)) return(NULL)
      rv$obj <- obj
      res$done     <- TRUE
      res$analysis <- analysis
      res$group_by <- group_by
      res$db       <- db
      res$species  <- species
      res$lfc      <- lfc
      res$summary  <- enrichment_summary(enrichment_table(obj, analysis, group_by))
      res$universe <- if (identical(analysis, "ora")) {
        obj@tools[[paste("Enrichment", group_by, "wilcox", sep = "_")]]$universe_size %||% NA
      } else {
        NA
      }
      mark_done(rv, "enrichment")
      log_step(log_rv, if (identical(analysis, "ora")) "Enrichment (ORA)" else "GSEA",
               params = list(analysis = analysis, group_by = group_by, db = db,
                             species = species,
                             lfc = if (identical(analysis, "ora")) lfc else NULL),
               code = enrichment_log_code(analysis, group_by, db, species, lfc = lfc))
      shiny::showNotification(
        i18n(sprintf("%s finished on '%s' (%s).",
                     if (identical(analysis, "ora")) "Enrichment" else "GSEA", group_by, db),
             sprintf("已完成 %s（分组：%s，数据库：%s）。",
                     if (identical(analysis, "ora")) "富集分析" else "GSEA", group_by, db)),
        type = "message")
    })

    output$summary <- shiny::renderUI({
      if (!isTRUE(res$done)) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("Choose a grouping and a database, then click <b>Run enrichment</b>.",
                               "选择分组与数据库，然后点击<b>运行富集分析</b>。")))
      }
      s <- res$summary
      shiny::tagList(
        stat_tile(i18n("Analysis", "分析类型"),
                  if (identical(res$analysis, "ora")) "ORA" else "GSEA"),
        stat_tile(i18n("Database", "数据库"), res$db),
        stat_tile(i18n("Group by", "分组依据"), res$group_by),
        stat_tile(i18n("Sets with q < 0.05", "q < 0.05 的基因集"),
                  format(s$n_sig, big.mark = ","))
      )
    })

    output$insight <- shiny::renderUI({
      if (!isTRUE(res$done)) return(NULL)
      s <- res$summary
      bg_en <- if (identical(res$analysis, "ora")) {
        sprintf(" Background: %s annotated genes that were tested in the DE step.",
                format(res$universe, big.mark = ","))
      } else {
        " Each group's ranking uses every tested gene (avg log2FC)."
      }
      bg_zh <- if (identical(res$analysis, "ora")) {
        sprintf("背景：DE 检验中参与检验且有注释的 %s 个基因。", format(res$universe, big.mark = ","))
      } else {
        "每组的排序使用全部参与检验的基因（avg log2FC）。"
      }
      insight_bar(
        sprintf(paste0("<b>%s</b> of %s gene-set tests reach BH q &lt; 0.05, in <b>%d</b> of %d ",
                       "groups.%s q-values are adjusted within each group &times; database, ",
                       "not across groups; the DE is cell-level, so read hits as marker ",
                       "biology, not as a between-condition result."),
                format(s$n_sig, big.mark = ","), format(s$n_tested, big.mark = ","),
                s$groups_hit, s$groups, bg_en),
        sprintf(paste0("共 %s 次基因集检验中有 <b>%s</b> 次达到 BH q &lt; 0.05，涉及 %d 个分组中的 ",
                       "<b>%d</b> 个。%sq 值只在每个“分组 × 数据库”内校正，未跨分组校正；",
                       "差异检验以细胞为单位，结果应理解为标志基因层面的生物学，而非条件间比较的结论。"),
                format(s$n_tested, big.mark = ","), format(s$n_sig, big.mark = ","),
                s$groups, s$groups_hit, bg_zh))
    })

    render_step_plot(output, input, "preview", function() {
      shiny::req(res$done, res$group_by)
      sc_enrichplot(rv$obj, analysis = res$analysis, group_by = res$group_by, db = res$db)
    }, name = "enrichment")
  })
}
