#' Module: Cell-cell communication
#'
#' Infer ligand-receptor signalling between cell groups (cell types / clusters)
#' using one of several backends. LIANA and CellChat run in R; CellPhoneDB and
#' NicheNet need extra (often Python) setup and are marked accordingly.
#' The result is stored in `rv$cellcomm` (not the Seurat object).
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
      en = "Infer which cell groups are talking to which, via ligand-receptor
            pairs, and how strong each interaction is.",
      zh = "推断哪些细胞群体在通过配体-受体对相互通讯，以及每种相互作用的强度。"),
    why  = list(
      en = "Tissue function emerges from crosstalk between cell types. Mapping
            ligand-receptor signalling reveals the wiring behind niches, immune
            responses and the tumour microenvironment.",
      zh = "组织功能源于细胞类型之间的相互作用。绘制配体-受体信号可揭示微环境、免疫应答
            和肿瘤微环境背后的连接关系。"),
    how  = list(
      en = "Pick the metadata column that labels your cell groups (usually cell
            type or cluster), then a method. <b>LIANA</b> aggregates several
            scoring methods and is a robust R default; <b>CellChat</b> adds
            pathway-level views. Methods marked <b>*</b> (CellPhoneDB, NicheNet)
            need extra setup (often Python).",
      zh = "选择标注细胞群体的元数据列（通常是细胞类型或簇），再选择方法。
            <b>LIANA</b> 聚合多种打分方法，是稳健的 R 默认；<b>CellChat</b> 提供通路级视图。
            标 <b>*</b> 的方法（CellPhoneDB、NicheNet）需要额外设置（通常是 Python）。"),
    read = list(
      en = "Circles are cell types, edges are significant ligand–receptor
            pairs; a thicker edge means stronger communication. Look for the
            sender→receiver pairs that could explain a behaviour you see.",
      zh = "圆圈为细胞类型，连线为显著的配体–受体对；线越粗通讯越强。寻找能解释你所观察现象的发送→接收对。"),
    example = list(
      en = "Macrophages signalling to T cells via a checkpoint ligand-receptor
               pair would appear as a strong edge between those two groups.",
      zh = "巨噬细胞通过某个免疫检查点配体-受体对向 T 细胞发出信号，会表现为两群之间的一条强边。")
  )
  controls <- shiny::tagList(
    label_with_help("Group-by column",
                    "The metadata column labelling the cell groups to test for communication (e.g. cell type, cluster).",
                    label_zh = "分组列",
                    tip_zh = "用于检验通讯的细胞群体标注元数据列（如细胞类型、簇）。"),
    shiny::uiOutput(ns("group_ui")),
    label_with_help("Method",
                    "LIANA / CellChat run in R. * = extra setup (often Python) required.",
                    label_zh = "方法",
                    tip_zh = "LIANA / CellChat 在 R 中运行。* = 需要额外设置（通常是 Python）。"),
    shiny::selectInput(ns("method"), NULL,
                       choices = c("LIANA" = "liana", "CellChat" = "cellchat",
                                   "CellPhoneDB *" = "cellphonedb",
                                   "NicheNet *" = "nichenet")),
    shiny::div(class = "omicone-note",
               i18n("* CellPhoneDB and NicheNet require extra setup (often a Python environment) and are not run in-app.",
                    "* CellPhoneDB 和 NicheNet 需要额外设置（通常是 Python 环境），不在应用内运行。")),
    run_button(ns("run"), "Infer communication", "推断通讯")
  )
  step_container(title = list(en = "Cell-cell communication", zh = "细胞间通讯"),
                 subtitle = list(en = "Ligand–receptor conversations between cell types.",
                                 zh = "细胞类型之间的配体–受体对话。"),
                 explainer = explainer, controls = controls,
                 summary = shiny::uiOutput(ns("summary")),
                 preview = preview_plot_ui(ns("preview"),
                   guide = list(en = "The inferred communication network will be drawn here.",
                                zh = "运行后，这里将绘制推断的通讯网络。"),
                   caption = list(en = "Cell-cell communication network (edge width = strength).",
                                  zh = "细胞间通讯网络（连线粗细＝强度）。")))
}

#' @rdname mod_cellcomm
#' @keywords internal
mod_cellcomm_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    res <- shiny::reactiveValues(done = FALSE, method = NULL, group = NULL,
                                 n_interactions = NA)

    # Populate the group-by selector from the object metadata.
    output$group_ui <- shiny::renderUI({
      cols <- obj_meta_cols(rv$obj)
      if (length(cols) == 0) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("Load and annotate a dataset to choose a group column.",
                               "加载并注释数据集以选择分组列。")))
      }
      # Prefer an annotation-like column if present.
      pref <- cols[cols %in% c("cell_type", "celltype", "CellType", "SingleR",
                               "seurat_clusters")]
      default <- if (length(pref)) pref[1] else cols[1]
      shiny::selectInput(session$ns("group"), NULL, choices = cols, selected = default)
    })

    shiny::observeEvent(input$run, {
      shiny::req(rv$obj)
      shiny::req(input$group)
      method <- input$method
      group  <- input$group
      # CellPhoneDB / NicheNet are not runnable in-app.
      if (method %in% c("cellphonedb", "nichenet")) {
        shiny::showNotification(
          i18n(sprintf("'%s' needs extra external setup (often Python) and is not run in-app.", method),
               sprintf("'%s' 需要额外的外部设置（通常是 Python），不在应用内运行。", method)),
          type = "warning", duration = 10)
        return(NULL)
      }
      pkg <- if (method == "liana") "liana" else "CellChat"
      if (!require_pkgs(pkg, "Cell-cell communication")) return(NULL)
      result <- with_progress_notify({
        sc_cellcomm(rv$obj, group_by = group, method = method)
      }, message = "Inferring cell-cell communication...")
      if (is.null(result)) return(NULL)
      rv$cellcomm <- list(method = method, group = group, result = result)
      res$done   <- TRUE
      res$method <- method
      res$group  <- group
      res$n_interactions <- tryCatch({
        if (is.data.frame(result)) nrow(result) else NA_integer_
      }, error = function(e) NA_integer_)
      mark_done(rv, "cellcomm")
      log_step(log_rv, "Cell-cell communication",
               params = list(method = method, group_by = group),
               code = sprintf("cellcomm <- sc_cellcomm(obj, group_by = '%s', method = '%s')",
                              group, method))
      shiny::showNotification(
        i18n(sprintf("Communication inference done (%s).", method),
             sprintf("通讯推断完成（%s）。", method)),
        type = "message")
    })

    output$summary <- shiny::renderUI({
      if (!isTRUE(res$done)) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("Choose a group column and method, then click Infer communication.",
                               "选择分组列和方法，然后点击推断通讯。")))
      }
      n_lab <- if (is.na(res$n_interactions)) "-" else
        format(res$n_interactions, big.mark = ",")
      bslib::layout_columns(
        col_widths = c(4, 4, 4),
        stat_tile(i18n("Method", "方法"), res$method),
        stat_tile(i18n("Group column", "分组列"), res$group),
        stat_tile(i18n("Interactions", "相互作用数"), n_lab)
      )
    })

    output$preview <- render_scop_plot(function() {
      shiny::req(res$done)
      cc <- rv$cellcomm
      shiny::req(cc)
      tryCatch(
        cellcomm_summary_plot(cc$result, cc$method),
        error = function(e) {
          shiny::showNotification(
            i18n(paste("Communication preview unavailable:", conditionMessage(e)),
                 paste("通讯预览不可用：", conditionMessage(e))),
            type = "error", duration = 10)
          NULL
        })
    })
  })
}

#' Build a robust summary plot from a cell-communication result
#'
#' Tries to render a source->target interaction-count heatmap. For LIANA the
#' result is a data.frame (or a named list of them) with `source`/`target`
#' columns; for CellChat we fall back to its own netVisual heatmap when
#' available. Any failure surfaces via the caller's tryCatch.
#' @keywords internal
cellcomm_summary_plot <- function(result, method) {
  method <- tolower(method)

  # LIANA: aggregate to a source x target interaction-count matrix.
  if (method == "liana") {
    df <- result
    # liana_wrap can return a named list of per-method data.frames.
    if (is.list(df) && !is.data.frame(df)) {
      hit <- Filter(function(x) is.data.frame(x) &&
                      all(c("source", "target") %in% colnames(x)), df)
      if (length(hit)) df <- hit[[1]]
    }
    if (!is.data.frame(df) || !all(c("source", "target") %in% colnames(df))) {
      stop("LIANA result has no source/target columns to summarise.")
    }
    counts <- as.data.frame(table(source = df$source, target = df$target),
                            stringsAsFactors = FALSE)
    counts$Freq <- as.numeric(counts$Freq)
    return(
      ggplot2::ggplot(counts, ggplot2::aes(x = target, y = source, fill = Freq)) +
        ggplot2::geom_tile(colour = "white") +
        ggplot2::scale_fill_gradient(low = "#eef3f8", high = "#2f81c7",
                                     name = "interactions") +
        ggplot2::labs(x = "Target (receiver) / 靶细胞（接收方）",
                      y = "Source (sender) / 源细胞（发送方）",
                      title = "Ligand-receptor interactions / 配体-受体相互作用") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
        omicone_theme()
    )
  }

  # CellChat: use its own network heatmap if the pipeline has been run.
  if (method == "cellchat") {
    if (has_pkg("CellChat")) {
      p <- tryCatch(CellChat::netVisual_heatmap(result), error = function(e) NULL)
      if (!is.null(p)) return(p)
    }
    stop("CellChat object needs the full inference pipeline before plotting.")
  }

  stop("No preview available for method '", method, "'.")
}
