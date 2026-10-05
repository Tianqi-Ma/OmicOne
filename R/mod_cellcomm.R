#' Module: Cell-cell communication
#'
#' Infer ligand-receptor signalling between cell groups with LIANA (consensus
#' of several scoring methods, aggregated with liana_aggregate()) or CellChat
#' (full standard workflow). CellPhoneDB and NicheNet need external setup and
#' are not run in-app. The result is stored in `rv$cellcomm`, not in the
#' Seurat object.
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_cellcomm
NULL

#' @rdname mod_cellcomm
#' @keywords internal
mod_cellcomm_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = list(en = "Cell-cell communication", zh = "细胞间通讯"),
    what = list(
      en = "Infer which cell groups signal to which through ligand-receptor (LR)
            pairs.",
      zh = "推断哪些细胞群体通过配体-受体（LR）对向哪些群体发出信号。"),
    why  = list(
      en = "Tissue function emerges from crosstalk between cell types; LR inference
            proposes the wiring behind niches, immune responses and the tumour
            microenvironment. It predicts potential signalling from expression,
            not measured binding.",
      zh = "组织功能源于细胞类型之间的相互作用；LR 推断给出微环境、免疫应答和肿瘤微环境背后可能的连接。它是基于表达对潜在信号的预测，而非测得的结合。"),
    how  = list(
      en = "Pick the column that labels the cell groups, then a method.
            <b>LIANA</b> runs several scoring methods and ranks LR pairs by their
            consensus (liana_aggregate); pairs with aggregate rank &lt; 0.05 are
            counted. <b>CellChat</b> runs its full workflow and keeps pairs with
            permutation p &lt; 0.05. Mouse data use the mouse databases.",
      zh = "选择标注细胞群体的列，再选择方法。<b>LIANA</b> 运行多种打分方法，并按共识对 LR 对排序（liana_aggregate）；统计 aggregate rank &lt; 0.05 的配对。<b>CellChat</b> 运行其完整流程，保留置换检验 p &lt; 0.05 的配对。小鼠数据使用小鼠数据库。"),
    read = list(
      en = "The heatmap counts significant LR pairs from each sender (row) to each
            receiver (column). A count is not a strength: one strong pair can
            matter more than many weak ones, and groups with more cells reach
            significance more easily.",
      zh = "热图统计每个发送方（行）到每个接收方（列）的显著 LR 对数量。数量不等于强度：一个强配对可能比许多弱配对更重要，细胞数更多的群体也更容易达到显著。"),
    example = list(
      en = "Macrophages signalling to T cells through a checkpoint ligand appear as
            a high count in the macrophage row, T-cell column.",
      zh = "巨噬细胞通过检查点配体向 T 细胞发出信号，会表现为巨噬细胞行、T 细胞列中较高的计数。")
  )
  controls <- shiny::tagList(
    label_with_help("Group-by column",
                    "Metadata column labelling the cell groups tested for communication (e.g. cell type, cluster).",
                    label_zh = "分组列",
                    tip_zh = "用于检验通讯的细胞群体标注列（如细胞类型、簇）。"),
    shiny::selectInput(ns("group"), NULL, choices = NULL),
    label_with_help("Method",
                    "LIANA / CellChat run in R. * = extra setup (often Python) required, not run in-app.",
                    label_zh = "方法",
                    tip_zh = "LIANA / CellChat 在 R 中运行。* = 需要额外设置（通常是 Python），不在应用内运行。"),
    shiny::selectInput(ns("method"), NULL,
                       choices = c("LIANA" = "liana", "CellChat" = "cellchat",
                                   "CellPhoneDB *" = "cellphonedb",
                                   "NicheNet *" = "nichenet")),
    shiny::div(class = "omicone-note",
               i18n("* CellPhoneDB and NicheNet require extra setup (often a Python environment) and are not run in-app.",
                    "* CellPhoneDB 和 NicheNet 需要额外设置（通常是 Python 环境），不在应用内运行。")),
    run_button(ns("run"), "Infer communication", "推断通讯")
  )
  step_container(id = id, title = list(en = "Cell-cell communication", zh = "细胞间通讯"),
                 subtitle = list(en = "Ligand-receptor signalling between cell types.",
                                 zh = "细胞类型之间的配体-受体信号。"),
                 explainer = explainer, controls = controls,
                 summary = shiny::uiOutput(ns("summary")),
                 preview = shiny::tagList(
                   shiny::uiOutput(ns("insight")),
                   preview_plot_ui(ns("preview"), download = TRUE,
                     guide = list(en = "A sender x receiver count heatmap will be drawn here.",
                                  zh = "运行后，这里将绘制发送方 × 接收方的计数热图。"),
                     caption = list(en = "Number of significant ligand-receptor pairs from each sender group (row) to each receiver group (column).",
                                    zh = "每个发送方分组（行）到每个接收方分组（列）的显著配体-受体对数量。"))))
}

#' @rdname mod_cellcomm
#' @keywords internal
mod_cellcomm_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    res <- step_results(rv, "sc", done = FALSE, method = NULL, group = NULL,
                        species = NULL, table = NULL, n_groups = NA)

    shiny::observe({
      obj <- rv$obj
      cols <- categorical_cols(obj_meta(obj))
      def <- default_group_col(cols, obj_misc(obj, "omicone_cluster_col"))
      shiny::updateSelectInput(session, "group", choices = cols,
                               selected = keep_selected(shiny::isolate(input$group), cols, def))
    })

    shiny::observeEvent(input$run, {
      shiny::req(rv$obj)
      method <- input$method
      group <- input$group
      shiny::req(method, group)
      if (method %in% c("cellphonedb", "nichenet")) {
        shiny::showNotification(
          i18n(sprintf("'%s' needs extra external setup (often Python) and is not run in-app.", method),
               sprintf("'%s' 需要额外的外部设置（通常是 Python），不在应用内运行。", method)),
          type = "warning", duration = 10)
        return(NULL)
      }
      pkgs <- if (identical(method, "liana")) c("liana", "Seurat") else c("CellChat", "Seurat")
      if (!require_pkgs(pkgs, "Cell-cell communication")) return(NULL)
      species <- guess_species(rv$obj)
      labels <- obj_meta(rv$obj)[[group]]
      prefixed <- identical(method, "cellchat") && any(as.character(labels) == "0", na.rm = TRUE)
      do_fast <- has_pkg("presto")
      out <- with_progress_notify({
        result <- if (identical(method, "liana")) {
          sc_liana(rv$obj, group_by = group, species = species)
        } else {
          sc_cellchat(rv$obj, group_by = group, species = species, do_fast = do_fast)
        }
        list(result = result, table = cellcomm_sig_table(result, method))
      }, message = "Inferring cell-cell communication...")
      if (is.null(out)) return(NULL)
      rv$cellcomm <- list(method = method, group = group, result = out$result,
                          table = out$table)
      res$done     <- TRUE
      res$method   <- method
      res$group    <- group
      res$species  <- species
      res$table    <- out$table
      res$n_groups <- length(unique(stats::na.omit(as.character(labels))))
      mark_done(rv, "cellcomm")
      log_step(log_rv, "Cell-cell communication",
               params = list(method = method, group_by = group, species = species,
                             threshold = if (identical(method, "liana"))
                               "aggregate_rank < 0.05" else "permutation p < 0.05"),
               code = if (identical(method, "liana")) {
                 liana_log_code(group, species)
               } else {
                 cellchat_log_code(group, species, prefixed = prefixed, do_fast = do_fast)
               })
      shiny::showNotification(
        i18n(sprintf("Communication inference done (%s).", method),
             sprintf("通讯推断完成（%s）。", method)),
        type = "message")
    })

    output$summary <- shiny::renderUI({
      if (!isTRUE(res$done)) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("Choose a group column and method, then click <b>Infer communication</b>.",
                               "选择分组列和方法，然后点击<b>推断通讯</b>。")))
      }
      shiny::tagList(
        stat_tile(i18n("Method", "方法"), res$method),
        stat_tile(i18n("Group column", "分组列"), res$group),
        stat_tile(i18n("Significant LR pairs", "显著 LR 对"),
                  format(nrow(res$table), big.mark = ","))
      )
    })

    output$insight <- shiny::renderUI({
      if (!isTRUE(res$done)) return(NULL)
      tab <- res$table
      n_pairs <- nrow(unique(tab[, c("ligand", "receptor"), drop = FALSE]))
      rule_en <- if (identical(res$method, "liana")) {
        "LIANA aggregate rank &lt; 0.05 (a rank-aggregation score, not FDR-adjusted)"
      } else {
        "CellChat permutation p &lt; 0.05 per pair (not FDR-adjusted)"
      }
      rule_zh <- if (identical(res$method, "liana")) {
        "LIANA aggregate rank &lt; 0.05（秩聚合分数，未做 FDR 校正）"
      } else {
        "CellChat 每对置换检验 p &lt; 0.05（未做 FDR 校正）"
      }
      insight_bar(
        sprintf("<b>%s</b> significant sender-receiver-LR rows (%s distinct LR pairs) among %d groups, by %s. Cells from all samples are pooled; treat these as hypotheses for validation.",
                format(nrow(tab), big.mark = ","), format(n_pairs, big.mark = ","),
                res$n_groups, rule_en),
        sprintf("%d 个分组之间共 <b>%s</b> 条显著的 发送方-接收方-LR 记录（%s 个不同的 LR 对），判定标准为 %s。所有样本的细胞被合并分析；请将结果视为有待验证的假设。",
                res$n_groups, format(nrow(tab), big.mark = ","), format(n_pairs, big.mark = ","),
                rule_zh))
    })

    render_step_plot(output, input, "preview", function() {
      shiny::req(res$done, res$table)
      cellcomm_count_plot(res$table)
    }, name = "cellcomm")
  })
}
