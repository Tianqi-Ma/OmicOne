#' Module: Marker genes
#'
#' Find genes that are differentially expressed in each cluster relative to the
#' rest of the cells. These marker genes are what you use to give a cluster a
#' biological identity in the next (annotation) step.
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_markers
NULL

#' @rdname mod_markers
#' @keywords internal
mod_markers_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = list(en = "Marker genes", zh = "标志基因"),
    what = list(
      en = "Find the genes that are specifically up- (or down-) regulated in each
            cluster compared with all other cells.",
      zh = "找出与其他所有细胞相比，在每个簇中特异性上调（或下调）的基因。"),
    why  = list(
      en = "Clusters are just groups of similar cells until you know what makes
            them different. Marker genes are the evidence you use to call a
            cluster a cell type (e.g. <code>CD3D</code> for T cells).",
      zh = "在你弄清各簇之间的区别之前，簇只是相似细胞的分组而已。标志基因是你把某个簇判定为某种细胞类型的依据（例如 <code>CD3D</code> 对应 T 细胞）。"),
    how  = list(
      en = "<b>Wilcoxon</b> is the fast, robust default. Raise the log fold-change
            or min.pct to keep only stronger, more specific markers. Keep
            <b>only positive</b> markers if you only care about what a cluster
            expresses <i>more</i> than others.",
      zh = "<b>Wilcoxon</b> 是快速、稳健的默认差异检验。提高 log fold-change 或 min.pct 可只保留更强、更特异的标志基因。若只关心某个簇比其他簇<i>更高</i>表达的基因，可只保留<b>正向</b>标志基因。"),
    read = list(
      en = "The plot shows the top markers per cluster: dot size = fraction of
            cells expressing the gene, colour = mean expression. A good marker
            is strong in one cluster and absent everywhere else.",
      zh = "图中展示每个簇的头部标志基因：点大小＝表达该基因的细胞比例，颜色＝平均表达量。好的标志基因在一个簇中强表达、在其他簇中缺失。"),
    example = list(
      en = "For a T-cell cluster you would expect markers like <code>CD3D</code>,
               <code>CD3E</code> and <code>TRAC</code> at the top of the list.",
      zh = "对于一个 T 细胞簇，你会期望 <code>CD3D</code>、<code>CD3E</code> 和 <code>TRAC</code> 这类标志基因排在列表前列。")
  )
  controls <- shiny::tagList(
    label_with_help("Statistical test",
                    "Wilcoxon = fast rank test (default). ROC = ranks genes by classification power. MAST = models dropout (needs the MAST package).",
                    "差异检验方法",
                    "Wilcoxon = 快速的秩检验（默认）。ROC = 按分类能力对基因排序。MAST = 对 dropout 建模（需要 MAST 包）。"),
    shiny::selectInput(ns("test"), NULL,
                       choices = c("Wilcoxon" = "wilcox", "ROC" = "roc", "MAST" = "MAST"),
                       selected = "wilcox"),
    label_with_help("Log fold-change threshold",
                    "Minimum log2 fold-change to test a gene. Higher = fewer, stronger markers.",
                    "Log fold-change 阈值",
                    "检验某个基因所需的最小 log2 fold-change。越高 = 标志基因越少、越强。"),
    shiny::numericInput(ns("logfc"), NULL, value = 0.25, min = 0, step = 0.05),
    label_with_help("Min fraction expressing (min.pct)",
                    "A gene must be detected in at least this fraction of cells in one of the two groups.",
                    "最小表达比例（min.pct）",
                    "某个基因必须在两组之一中至少这一比例的细胞里被检测到。"),
    shiny::numericInput(ns("min_pct"), NULL, value = 0.1, min = 0, max = 1, step = 0.05),
    shiny::checkboxInput(ns("only_pos"),
                         i18n("Only positive markers", "仅保留正向标志基因"),
                         value = TRUE),
    label_with_help("Top N per cluster",
                    "How many top markers per cluster to show in the preview chart.",
                    "每簇 Top N",
                    "在预览图中每个簇显示多少个排名靠前的标志基因。"),
    shiny::numericInput(ns("top_n"), NULL, value = 10, min = 1, max = 50, step = 1),
    run_button(ns("run"), "Find markers", "查找标志基因")
  )
  step_container(
    title     = list(en = "Marker genes", zh = "标志基因"),
    subtitle  = list(en = "Find the genes that define each cluster.",
                     zh = "找出定义每个簇的基因。"),
    explainer = explainer,
    controls  = controls,
    summary   = shiny::uiOutput(ns("summary")),
    preview   = shiny::tagList(
      preview_plot_ui(ns("preview"),
        guide = list(en = "Top markers per cluster will be drawn here.",
                     zh = "运行后，这里将绘制每个簇的头部标志基因。"),
        caption = list(en = "Top markers per cluster (dot size = % expressing, colour = mean level).",
                       zh = "各簇头部标志基因（点大小＝表达比例，颜色＝平均表达量）。")),
      shiny::div(class = "omicone-table", shiny::uiOutput(ns("table")))
    )
  )
}

#' @rdname mod_markers
#' @keywords internal
mod_markers_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    markers <- shiny::reactiveVal(NULL)

    shiny::observeEvent(input$run, {
      shiny::req(rv$obj)
      if (!require_pkgs("Seurat", "Marker genes")) return(NULL)
      df <- with_progress_notify({
        markers_obj(rv$obj, test = input$test, logfc = input$logfc,
                    min_pct = input$min_pct, only_pos = input$only_pos)
      }, message = "Finding marker genes...")
      if (is.null(df)) return(NULL)
      markers(df)
      rv$markers <- df
      mark_done(rv, "markers")
      log_step(log_rv, "Markers",
               params = list(test = input$test, logfc = input$logfc,
                             min_pct = input$min_pct, only_pos = input$only_pos),
               code = sprintf(
                 'markers <- Seurat::FindAllMarkers(obj, test.use="%s", logfc.threshold=%s, min.pct=%s, only.pos=%s)',
                 input$test, input$logfc, input$min_pct, input$only_pos))
      shiny::showNotification(sprintf("Found %d markers across %d clusters.",
                                      nrow(df), length(unique(df$cluster))),
                              type = "message")
    })

    output$summary <- shiny::renderUI({
      df <- markers()
      if (is.null(df)) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("Set options and click <b>Find markers</b>.",
                               "设置选项后点击<b>查找标志基因</b>。")))
      }
      n_clusters <- length(unique(df$cluster))
      genes <- if ("gene" %in% colnames(df)) length(unique(df$gene)) else NA_integer_
      bslib::layout_columns(
        col_widths = c(4, 4, 4),
        stat_tile(i18n("Markers", "标志基因数"), format(nrow(df), big.mark = ",")),
        stat_tile(i18n("Clusters", "簇数"), format(n_clusters, big.mark = ",")),
        stat_tile(i18n("Unique genes", "去重基因数"), format(genes, big.mark = ","))
      )
    })

    # Top-N markers per cluster as a tidy data.frame (shared by chart + table)
    top_markers <- shiny::reactive({
      df <- markers()
      shiny::req(df)
      n <- max(1, as.integer(input$top_n))
      parts <- split(df, df$cluster)
      picked <- lapply(parts, function(d) {
        ord <- if ("avg_log2FC" %in% colnames(d)) order(-d$avg_log2FC) else seq_len(nrow(d))
        utils::head(d[ord, , drop = FALSE], n)
      })
      do.call(rbind, picked)
    })

    output$preview <- render_preview_plot(function() {
      d <- top_markers()
      shiny::req(d)
      d$cluster <- factor(d$cluster)
      lfc <- if ("avg_log2FC" %in% colnames(d)) d$avg_log2FC else rep(NA_real_, nrow(d))
      padj <- if ("p_val_adj" %in% colnames(d)) d$p_val_adj else rep(NA_real_, nrow(d))
      gene <- if ("gene" %in% colnames(d)) d$gene else rownames(d)
      d$gene_lab <- gene
      d$avg_log2FC <- lfc
      d$p_val_adj  <- padj
      d$text <- sprintf("gene: %s\ncluster: %s\navg_log2FC: %.2f\np_val_adj: %.2g",
                        gene, as.character(d$cluster), lfc, padj)
      ggplot2::ggplot(d, ggplot2::aes(x = stats::reorder(gene_lab, avg_log2FC),
                                      y = avg_log2FC, fill = cluster)) +
        ggplot2::geom_col() +
        ggplot2::coord_flip() +
        ggplot2::facet_wrap(~cluster, scales = "free_y") +
        ggplot2::scale_fill_manual(values = sc_palette(nlevels(d$cluster)),
                                   guide = "none") +
        ggplot2::labs(x = NULL, y = "avg log2 fold-change",
                      title = "Top marker genes per cluster") +
        omicone_theme()
    })

    output$table <- shiny::renderUI({
      shiny::req(markers())
      ns <- session$ns
      if (has_pkg("DT")) {
        DT::dataTableOutput(ns("dt"))
      } else {
        shiny::verbatimTextOutput(ns("txt"))
      }
    })

    if (has_pkg("DT")) {
      output$dt <- DT::renderDataTable({
        df <- markers()
        shiny::req(df)
        DT::datatable(df, filter = "top", rownames = FALSE,
                      options = list(pageLength = 10, scrollX = TRUE))
      })
    } else {
      output$txt <- shiny::renderPrint({
        df <- markers()
        shiny::req(df)
        utils::head(df, 20)
      })
    }
  })
}
