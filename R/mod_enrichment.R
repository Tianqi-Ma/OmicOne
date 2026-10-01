#' Module: Enrichment & GSEA
#'
#' Turn a list of differentially-expressed genes into interpretable biology by
#' testing which annotated gene sets (GO terms, KEGG/Reactome/WikiPathway
#' pathways, MSigDB signatures) are over-represented (ORA) or coordinately
#' shifted (GSEA) in each group. Wraps scop::RunEnrichment / scop::RunGSEA.
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
      en = "Summarise the marker/DE genes of each group into the biological
            pathways and functions they represent.",
      zh = "把每个分组的标志基因 / 差异表达基因归纳为它们所代表的生物学通路与功能。"),
    why  = list(
      en = "A long list of gene names is hard to interpret. Enrichment tells you,
            for example, that a cluster is enriched for <i>antigen presentation</i>
            or <i>oxidative phosphorylation</i> -- turning genes into biology.",
      zh = "一长串基因名很难解读。富集分析能告诉你，例如某个簇富集于<i>抗原呈递</i>或<i>氧化磷酸化</i>——把基因转化为生物学意义。"),
    how  = list(
      en = "<b>Over-representation (ORA)</b> tests whether a group's marker genes
            hit a pathway more than by chance. <b>GSEA</b> instead uses the full
            ranked gene list, so it can detect coordinated weak shifts. Pick a
            database and matching species. <b>Run the Markers step first</b> so
            the DE results are available.",
      zh = "<b>过表达分析（ORA）</b>检验某分组的标志基因命中某通路是否超出随机预期。<b>GSEA</b> 则使用完整的排序基因列表，因此能检测到协同的微弱变化。请选择数据库并匹配物种。<b>请先运行“标志基因”步骤</b>，以便获得差异表达结果。"),
    read = list(
      en = "Each bar or point is a gene set, ranked by significance. Read the
            top hits together with the genes that drove them — a hit resting on
            a single gene is fragile.",
      zh = "每根条/每个点是一个基因集，按显著性排序。解读头部结果时要结合驱动它的基因——仅由单个基因支撑的结果不可靠。"),
    example = list(
      en = "For a cytotoxic T-cell cluster, GO terms like
               <code>T cell mediated cytotoxicity</code> should rise to the top.",
      zh = "对于一个细胞毒性 T 细胞簇，<code>T cell mediated cytotoxicity</code> 之类的 GO 条目应排在最前列。")
  )
  controls <- shiny::tagList(
    label_with_help("Analysis type",
                    "ORA tests marker-gene overlap with pathways; GSEA uses the full ranked gene list.",
                    label_zh = "分析类型",
                    tip_zh = "ORA 检验标志基因与通路的重叠；GSEA 使用完整的排序基因列表。"),
    shiny::radioButtons(ns("analysis"), NULL,
                        c("Over-representation (ORA)" = "ora", "GSEA" = "gsea"),
                        selected = "ora"),
    label_with_help("Group by (metadata column)",
                    "Which grouping the DE/enrichment is computed per (e.g. seurat_clusters, celltype).",
                    label_zh = "分组依据（元数据列）",
                    tip_zh = "按哪一列分组来计算差异 / 富集（例如 seurat_clusters、celltype）。"),
    shiny::selectInput(ns("group_by"), NULL, choices = NULL),
    label_with_help("Database",
                    "Source of gene sets: GO terms or curated pathway/signature collections.",
                    label_zh = "数据库",
                    tip_zh = "基因集来源：GO 条目或经过整理的通路 / 特征集合。"),
    shiny::selectInput(ns("db"), NULL,
                       choices = c("GO" = "GO", "KEGG" = "KEGG",
                                   "Reactome" = "Reactome",
                                   "WikiPathway" = "WikiPathway",
                                   "MSigDB" = "MSigDB"),
                       selected = "GO"),
    label_with_help("Species",
                    "Must match your data; sets the annotation database used.",
                    label_zh = "物种",
                    tip_zh = "必须与你的数据一致；决定所使用的注释数据库。"),
    shiny::selectInput(ns("species"), NULL,
                       choices = c("Homo sapiens" = "Homo_sapiens",
                                   "Mus musculus" = "Mus_musculus"),
                       selected = "Homo_sapiens"),
    run_button(ns("run"), "Run enrichment", "运行富集分析")
  )
  step_container(
    title     = list(en = "Enrichment & GSEA", zh = "富集与 GSEA"),
    subtitle  = list(en = "Pathways and gene sets behind your gene list.",
                     zh = "基因列表背后的通路与基因集。"),
    explainer = explainer,
    controls  = controls,
    summary   = shiny::uiOutput(ns("summary")),
    preview   = preview_plot_ui(ns("preview"),
      guide = list(en = "Enrichment results will be drawn here.",
                   zh = "运行后，这里将绘制富集分析结果。"),
      caption = list(en = "Enriched gene sets ranked by significance.",
                     zh = "按显著性排序的富集基因集。"))
  )
}

#' @rdname mod_enrichment
#' @keywords internal
mod_enrichment_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    res <- shiny::reactiveValues(done = FALSE, analysis = NULL, group_by = NULL,
                                 db = NULL, species = NULL)

    # Keep the group-by selector in sync with the current object.
    shiny::observe({
      obj <- rv$obj
      cols <- if (is.null(obj)) character(0) else obj_meta_cols(obj)
      sel <- if ("celltype" %in% cols) "celltype"
             else if ("seurat_clusters" %in% cols) "seurat_clusters"
             else if (length(cols)) cols[1] else NULL
      shiny::updateSelectInput(session, "group_by", choices = cols, selected = sel)
    })

    shiny::observeEvent(input$run, {
      shiny::req(rv$obj)
      # All enrichment engines live in scop.
      if (!require_pkgs("scop", "Enrichment / GSEA")) return(NULL)
      group_by <- input$group_by
      shiny::req(group_by)
      obj <- with_progress_notify({
        if (input$analysis == "ora") {
          sc_enrichment(rv$obj, group_by = group_by, db = input$db,
                        species = input$species)
        } else {
          sc_gsea(rv$obj, group_by = group_by, db = input$db,
                  species = input$species)
        }
      }, message = if (input$analysis == "ora")
        "Running over-representation analysis..." else "Running GSEA...")
      if (is.null(obj)) return(NULL)
      rv$obj <- obj
      res$done     <- TRUE
      res$analysis <- input$analysis
      res$group_by <- group_by
      res$db       <- input$db
      res$species  <- input$species
      mark_done(rv, "enrichment")
      log_step(log_rv, "Enrichment",
               params = list(analysis = input$analysis, group_by = group_by,
                             db = input$db, species = input$species),
               code = sprintf(
                 'obj <- scop::%s(obj, group_by="%s", db="%s", species="%s")',
                 if (input$analysis == "ora") "RunEnrichment" else "RunGSEA",
                 group_by, input$db, input$species))
      shiny::showNotification(
        i18n(sprintf("%s finished on '%s' (%s).",
                     if (input$analysis == "ora") "Enrichment" else "GSEA",
                     group_by, input$db),
             sprintf("已完成 %s（分组：%s，数据库：%s）。",
                     if (input$analysis == "ora") "富集分析" else "GSEA",
                     group_by, input$db)),
        type = "message")
    })

    output$summary <- shiny::renderUI({
      if (!isTRUE(res$done)) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("Run the Markers step first, then choose a database and click <b>Run enrichment</b>.",
                               "请先运行“标志基因”步骤，然后选择数据库并点击<b>运行富集分析</b>。")))
      }
      bslib::layout_columns(
        col_widths = c(4, 4, 4),
        stat_tile(i18n("Analysis", "分析类型"),
                  if (res$analysis == "ora") "ORA" else "GSEA"),
        stat_tile(i18n("Database", "数据库"), res$db),
        stat_tile(i18n("Group by", "分组依据"), res$group_by)
      )
    })

    output$preview <- render_scop_plot(function() {
      shiny::req(res$done, res$group_by)
      # ORA -> bar plot of enriched terms; GSEA -> running-score plot.
      if (identical(res$analysis, "ora")) {
        sc_enrichplot(rv$obj, group_by = res$group_by, plot_type = "bar")
      } else {
        sc_gseaplot(rv$obj, group_by = res$group_by)
      }
    })
  })
}
