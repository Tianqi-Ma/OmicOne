#' Module: Marker genes
#'
#' Find genes that are differentially expressed in each cluster (or any other
#' grouping) relative to the rest of the cells. These marker genes are what you
#' use to give a cluster a biological identity in the next (annotation) step.
#' This is cell-level marker discovery, not a between-condition test.
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
            cluster a cell type (e.g. <code>CD3D</code> for T cells). Cells are
            treated as replicates, so these p-values support marker discovery,
            not claims that conditions or patients differ.",
      zh = "在你弄清各簇之间的区别之前，簇只是相似细胞的分组而已。标志基因是你把某个簇判定为某种细胞类型的依据（例如 <code>CD3D</code> 对应 T 细胞）。这里把细胞当作重复，因此这些 p 值只用于寻找标志基因，不能用于声称不同条件或患者之间存在差异。"),
    how  = list(
      en = "<b>Wilcoxon</b> is the fast, robust default. Raise the log fold-change
            or min.pct to keep only stronger, more specific markers. Keep
            <b>only positive</b> markers if you only care about what a group
            expresses <i>more</i> than others. <b>Group by</b> defaults to the
            active clustering.",
      zh = "<b>Wilcoxon</b> 是快速、稳健的默认差异检验。提高 log fold-change 或 min.pct 可只保留更强、更特异的标志基因。若只关心某组比其他组<i>更高</i>表达的基因，可只保留<b>正向</b>标志基因。<b>分组</b>默认是当前聚类。"),
    read = list(
      en = "Dot plot of the top markers per group (BH-adjusted p < 0.05, then
            largest fold change): dot size = % of cells in the group expressing
            the gene, colour = scaled mean expression. A good marker is a big,
            dark dot in one group and small or pale elsewhere.",
      zh = "每组头部标志基因的点图（先筛选 BH 校正 p < 0.05，再按 fold change 排序）：点大小 = 组内表达该基因的细胞比例，颜色 = 标准化后的平均表达。好的标志基因在一个组中是又大又深的点，在其他组中小或浅。"),
    example = list(
      en = "For a T-cell cluster you would expect markers like <code>CD3D</code>,
               <code>CD3E</code> and <code>TRAC</code> at the top of the list.",
      zh = "对于一个 T 细胞簇，你会期望 <code>CD3D</code>、<code>CD3E</code> 和 <code>TRAC</code> 这类标志基因排在列表前列。")
  )
  tbl_out <- if (has_pkg("DT")) DT::dataTableOutput(ns("tbl")) else shiny::verbatimTextOutput(ns("tbl"))
  controls <- shiny::tagList(
    label_with_help("Group by",
                    "Groups compared one against the rest. Defaults to the active clustering; any categorical column (e.g. celltype) works.",
                    "分组",
                    "每组与其余细胞比较。默认使用当前聚类；任意分类列（例如 celltype）均可。"),
    shiny::uiOutput(ns("group_ui")),
    label_with_help("Statistical test",
                    "Wilcoxon = fast rank test (default). ROC = ranks genes by classification power (no p-values). MAST = hurdle model (needs the MAST package).",
                    "差异检验方法",
                    "Wilcoxon = 快速的秩检验（默认）。ROC = 按分类能力对基因排序（不给 p 值）。MAST = hurdle 模型（需要 MAST 包）。"),
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
    label_with_help("Top N per group",
                    "How many significant markers per group to show in the figures (display only).",
                    "每组 Top N",
                    "图中每组显示多少个显著的标志基因（仅影响显示）。"),
    shiny::numericInput(ns("top_n"), NULL, value = 5, min = 1, max = 50, step = 1),
    run_button(ns("run"), "Find markers", "查找标志基因")
  )
  step_container(id = id, 
    title     = list(en = "Marker genes", zh = "标志基因"),
    subtitle  = list(en = "Find the genes that define each cluster.",
                     zh = "找出定义每个簇的基因。"),
    explainer = explainer,
    controls  = controls,
    summary   = shiny::uiOutput(ns("summary")),
    preview   = shiny::tagList(
      shiny::uiOutput(ns("insight")),
      bslib::navset_card_tab(
        bslib::nav_panel(i18n("Dot plot", "点图"),
          preview_plot_ui(ns("preview"), download = TRUE,
            guide = list(en = "A dot plot of the top markers per group will be drawn here.",
                         zh = "运行后，这里将绘制每组头部标志基因的点图。"),
            caption = list(en = "Dot size = % of cells expressing; colour = scaled mean expression (top N per group, adjusted p < 0.05).",
                           zh = "点大小＝表达比例；颜色＝标准化平均表达（每组 Top N，校正 p < 0.05）。"))),
        bslib::nav_panel(i18n("Fold changes", "倍数变化"),
          preview_plot_ui(ns("bars"), download = TRUE,
            guide = list(en = "Bars of the top markers' fold changes will be drawn here.",
                         zh = "运行后，这里将绘制头部标志基因的倍数变化条形图。"),
            caption = list(en = "One bar = one marker's average log2 fold change in its group versus all other cells.",
                           zh = "每根条为一个标志基因在其所在组相对其他所有细胞的平均 log2 倍数变化。"))),
        bslib::nav_panel(i18n("Table", "表格"),
          shiny::div(class = "omicone-table",
                     shiny::downloadButton(ns("csv"), i18n("Download full table (CSV)",
                                                           "下载完整表格（CSV）"),
                                           class = "btn-sm"),
                     tbl_out))
      )
    )
  )
}

#' @rdname mod_markers
#' @keywords internal
mod_markers_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    res <- step_results(rv, "sc", df = NULL, group_by = NULL, test = NULL)

    output$group_ui <- shiny::renderUI({
      md <- obj_meta(rv$obj)
      cols <- sc_group_cols(md)
      if (!length(cols)) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("No grouping yet — run clustering first.",
                               "尚无分组——请先运行聚类。")))
      }
      active <- obj_misc(rv$obj, "omicone_cluster_col")
      default <- intersect(c(active, "seurat_clusters"), cols)
      default <- if (length(default)) default[1] else cols[1]
      shiny::selectInput(session$ns("group_by"), NULL, choices = cols,
                         selected = keep_selected(shiny::isolate(input$group_by), cols, default))
    })

    shiny::observeEvent(input$run, {
      shiny::req(rv$obj)
      test <- input$test
      group_by <- input$group_by
      logfc <- num_input(input$logfc, 0)
      min_pct <- num_input(input$min_pct, 0, 1)
      only_pos <- isTRUE(input$only_pos)
      if (is.na(logfc) || is.na(min_pct)) {
        shiny::showNotification(i18n("Enter valid numbers for the fold-change and min.pct thresholds.",
                                     "请为 fold change 和 min.pct 阈值输入有效数值。"), type = "error")
        return(NULL)
      }
      if (is.null(group_by) || !group_by %in% obj_meta_cols(rv$obj)) {
        shiny::showNotification(i18n("Choose a grouping column (run clustering first).",
                                     "请选择分组列（请先运行聚类）。"), type = "error")
        return(NULL)
      }
      pkgs <- c("Seurat", if (test == "MAST") "MAST")
      if (!require_pkgs(pkgs, "Marker genes")) return(NULL)
      df <- with_progress_notify({
        markers_obj(rv$obj, group_by = group_by, test = test, logfc = logfc,
                    min_pct = min_pct, only_pos = only_pos)
      }, message = "Finding marker genes...")
      if (is.null(df)) return(NULL)
      res$df <- df
      res$group_by <- group_by
      res$test <- test
      rv$markers <- df
      mark_done(rv, "markers")
      log_step(log_rv, "Markers",
               params = list(group_by = group_by, test = test, logfc = logfc,
                             min_pct = min_pct, only_pos = only_pos),
               code = markers_log_code(group_by, test, logfc, min_pct, only_pos))
      shiny::showNotification(
        i18n(sprintf("Found %s markers across %d groups.", format(nrow(df), big.mark = ","),
                     length(unique(df$cluster))),
             sprintf("在 %d 个组中找到 %s 个标志基因。", length(unique(df$cluster)),
                     format(nrow(df), big.mark = ","))),
        type = "message")
    })

    top_markers <- shiny::reactive({
      df <- res$df
      shiny::req(df)
      n <- int_input(input$top_n, 1, 200)
      shiny::req(!is.na(n))
      marker_top_n(df, n)
    })

    output$summary <- shiny::renderUI({
      df <- res$df
      if (is.null(df)) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("Set options and click <b>Find markers</b>.",
                               "设置选项后点击<b>查找标志基因</b>。")))
      }
      sig <- if (!is.null(df$p_val_adj)) sum(df$p_val_adj < 0.05, na.rm = TRUE) else NA
      shiny::tagList(
        stat_tile(i18n("Genes tested, kept", "检验后保留"), format(nrow(df), big.mark = ",")),
        stat_tile(i18n("Adjusted p < 0.05", "校正 p < 0.05"),
                  if (is.na(sig)) "n/a" else format(sig, big.mark = ",")),
        stat_tile(i18n("Groups", "组数"), format(length(unique(df$cluster)))),
        stat_tile(i18n("Grouped by", "分组"), res$group_by)
      )
    })

    output$insight <- shiny::renderUI({
      df <- res$df
      if (is.null(df)) return(NULL)
      n_grp <- length(unique(df$cluster))
      if (is.null(df$p_val_adj)) {
        return(insight_bar(
          sprintf("%s markers across %d groups ranked by %s; the ROC test gives no p-values. Cell-level marker discovery, not a comparison of conditions.",
                  format(nrow(df), big.mark = ","), n_grp, if (!is.null(df$power)) "power" else "AUC"),
          sprintf("%d 个组共 %s 个标志基因，按%s排序；ROC 检验不给 p 值。这是细胞层面的标志基因发现，不是条件间比较。",
                  n_grp, format(nrow(df), big.mark = ","), if (!is.null(df$power)) " power " else " AUC ")))
      }
      sig <- df[!is.na(df$p_val_adj) & df$p_val_adj < 0.05, , drop = FALSE]
      none <- setdiff(unique(as.character(df$cluster)), unique(as.character(sig$cluster)))
      insight_bar(
        sprintf("%s of %s markers have BH-adjusted p < 0.05 across %d groups%s. Cells are the unit of this test (marker discovery, not between-condition inference).",
                format(nrow(sig), big.mark = ","), format(nrow(df), big.mark = ","), n_grp,
                if (length(none)) sprintf("; no significant marker for %s", paste(none, collapse = ", ")) else ""),
        sprintf("%d 个组共 %s 个标志基因中，%s 个 BH 校正 p < 0.05%s。该检验以细胞为单位（用于发现标志基因，不能用于条件间推断）。",
                n_grp, format(nrow(df), big.mark = ","), format(nrow(sig), big.mark = ","),
                if (length(none)) sprintf("；%s 没有显著标志基因", paste(none, collapse = "、")) else ""))
    })

    render_step_plot(output, input, "preview", function() {
      top <- top_markers()
      shiny::req(nrow(top) > 0, rv$obj, res$group_by %in% obj_meta_cols(rv$obj))
      markers_dot_plot(rv$obj, top, res$group_by)
    }, name = "markers_dotplot", width = 12, height = 7)

    render_step_plot(output, input, "bars", function() {
      top <- top_markers()
      shiny::req(nrow(top) > 0)
      markers_bar_plot(top, colors = if (res$group_by %in% obj_meta_cols(rv$obj)) group_colors(rv$obj, res$group_by))
    }, name = "markers_bars", width = 12, height = 9)

    output$tbl <- render_tbl_wrap(function() {
      df <- res$df
      shiny::req(df)
      utils::head(df, 5000)   # cap rows for the browser; the CSV has them all
    })

    output$csv <- shiny::downloadHandler(
      filename = function() sprintf("markers_%s.csv", format(Sys.time(), "%Y%m%d_%H%M%S")),
      content = function(file) {
        df <- res$df
        if (is.null(df)) df <- data.frame()
        utils::write.csv(df, file, row.names = FALSE)
      }
    )
  })
}
