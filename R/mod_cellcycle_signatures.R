#' Module: Cell cycle & signature scoring
#'
#' Two related sub-actions on the working object:
#'  (a) Cell-cycle scoring — adds S.Score, G2M.Score and a Phase call.
#'  (b) Signature scoring — scores custom or built-in gene sets per cell
#'      (UCell or Seurat AddModuleScore).
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_cellcycle_signatures
NULL

#' Built-in example signature gene sets (small, illustrative)
#' @keywords internal
cellcycle_example_sets <- function() {
  list(
    proliferation = c("MKI67", "TOP2A", "PCNA", "CCNB1", "CCNB2", "CDK1",
                      "BIRC5", "AURKB", "CENPF"),
    EMT           = c("VIM", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1",
                      "FN1", "CDH2", "MMP2", "MMP9"),
    hypoxia       = c("HIF1A", "VEGFA", "CA9", "SLC2A1", "LDHA", "PGK1",
                      "ALDOA", "ENO1", "BNIP3"),
    inflammation  = c("IL6", "IL1B", "TNF", "CXCL8", "CCL2", "NFKB1",
                      "STAT1", "IRF1", "PTGS2")
  )
}

#' @rdname mod_cellcycle_signatures
#' @keywords internal
mod_cellcycle_signatures_ui <- function(id) {
  ns <- shiny::NS(id)
  sets <- cellcycle_example_sets()
  explainer <- explainer_card(
    title = list(en = "Cell cycle & signatures", zh = "细胞周期与信号评分"),
    what = list(
      en = "Score each cell for its cell-cycle phase and for the activity of gene
            signatures (e.g. proliferation, EMT, hypoxia, inflammation).",
      zh = "为每个细胞评估其细胞周期时相，以及基因信号（如增殖、EMT、缺氧、炎症）的活性。"),
    why  = list(
      en = "Cell-cycle differences can dominate clustering and be mistaken for
            biology; signature scores turn a curated gene list into a single,
            comparable per-cell activity value.",
      zh = "细胞周期差异可能主导聚类并被误认为生物学差异；信号评分将一份精选基因列表
            转化为单一、可比较的每细胞活性值。"),
    how  = list(
      en = "Click <b>Score cell cycle</b> to add S/G2M scores and a Phase call.
            For signatures, tick built-in sets and/or type your own as
            <code>SetName: GENE1, GENE2, ...</code> (one per line), pick a method,
            then <b>Score signatures</b>. <b>UCell</b> is rank-based and robust;
            <b>AddModuleScore</b> is the Seurat default.",
      zh = "点击<b>细胞周期评分</b>以添加 S/G2M 分数和 Phase 判定。
            对于信号评分，勾选内置基因集，和/或按
            <code>集合名: GENE1, GENE2, ...</code>（每行一个）自定义，选择方法后点击
            <b>信号评分</b>。<b>UCell</b> 基于排名且稳健；<b>AddModuleScore</b> 为 Seurat 默认。"),
    read = list(
      en = "Each cell gets one score per gene set; high scores on the embedding
            show where a programme is active. A cluster dominated by cell-cycle
            scores is a cycling state, not a separate cell type.",
      zh = "每个细胞对每个基因集得到一个分数；嵌入图上的高分区域即该程序活跃的位置。被周期分数主导的簇是一种增殖状态，而非独立的细胞类型。"),
    example = list(
      en = "A tumour cluster scoring high on proliferation and hypoxia while
               cycling in G2M points to an actively growing, oxygen-starved niche.",
      zh = "某肿瘤簇在增殖和缺氧上得分高且处于 G2M 周期，提示一个活跃增殖、缺氧的微环境。")
  )
  controls <- shiny::tagList(
    # (a) Cell-cycle scoring
    shiny::tags$h6(i18n("1. Cell cycle", "1. 细胞周期")),
    label_with_help("Cell-cycle scoring",
                    "Uses Seurat's updated 2019 S/G2M gene lists to assign each cell a phase.",
                    label_zh = "细胞周期评分",
                    tip_zh = "使用 Seurat 2019 更新的 S/G2M 基因列表为每个细胞分配时相。"),
    run_button(ns("run_cc"), "Score cell cycle", "细胞周期评分"),
    shiny::tags$hr(),
    # (b) Signature scoring
    shiny::tags$h6(i18n("2. Signature scoring", "2. 信号评分")),
    label_with_help("Built-in example sets",
                    "Curated illustrative gene sets. Tick any to include them.",
                    label_zh = "内置示例基因集",
                    tip_zh = "精选的示例基因集。勾选以纳入评分。"),
    shiny::checkboxGroupInput(ns("builtin"), NULL,
                              choices = stats::setNames(names(sets), names(sets))),
    label_with_help("Custom gene sets",
                    "One set per line, format: SetName: GENE1, GENE2, GENE3",
                    label_zh = "自定义基因集",
                    tip_zh = "每行一个，格式：集合名: GENE1, GENE2, GENE3"),
    shiny::textAreaInput(ns("custom"), NULL, value = "", rows = 4,
                         placeholder = "MySet: CD3D, CD3E, TRAC"),
    label_with_help("Method",
                    "UCell = rank-based, robust to depth (needs the UCell package). AddModuleScore = Seurat default.",
                    label_zh = "方法",
                    tip_zh = "UCell = 基于排名、对测序深度稳健（需要 UCell 包）。AddModuleScore = Seurat 默认。"),
    shiny::selectInput(ns("method"), NULL,
                       c("UCell" = "UCell", "AddModuleScore" = "AddModuleScore")),
    run_button(ns("run_sig"), "Score signatures", "信号评分"),
    shiny::tags$hr(),
    # Preview control
    label_with_help("Preview feature",
                    "Choose a computed score, or Phase, to display on the embedding.",
                    label_zh = "预览特征",
                    tip_zh = "选择要在嵌入上显示的评分或 Phase。"),
    shiny::uiOutput(ns("preview_ui"))
  )
  step_container(title = list(en = "Cell cycle & signatures", zh = "细胞周期与信号评分"),
                 subtitle = list(en = "Score cell-cycle phase and gene-set activity per cell.",
                                 zh = "为每个细胞评定细胞周期时相与基因集活性。"),
                 explainer = explainer, controls = controls,
                 summary = shiny::uiOutput(ns("summary")),
                 preview = preview_plot_ui(ns("preview"),
                   guide = list(en = "Per-cell signature scores will be drawn here.",
                                zh = "运行后，这里将绘制每细胞的信号评分。"),
                   caption = list(en = "Signature / cell-cycle scores over the embedding.",
                                  zh = "叠加在嵌入图上的信号 / 细胞周期评分。")))
}

#' @rdname mod_cellcycle_signatures
#' @keywords internal
mod_cellcycle_signatures_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    res <- shiny::reactiveValues(cc_done = FALSE, sig_done = FALSE,
                                 sets = NULL, method = NULL)

    # Parse the custom textarea into a named list of gene-set vectors.
    parse_custom_sets <- function(txt) {
      if (is.null(txt) || !nzchar(trimws(txt))) return(list())
      lines <- strsplit(txt, "\n", fixed = TRUE)[[1]]
      out <- list()
      for (ln in lines) {
        ln <- trimws(ln)
        if (!nzchar(ln) || !grepl(":", ln, fixed = TRUE)) next
        nm <- trimws(sub(":.*$", "", ln))
        genes <- trimws(strsplit(sub("^[^:]*:", "", ln), ",", fixed = TRUE)[[1]])
        genes <- genes[nzchar(genes)]
        if (nzchar(nm) && length(genes) > 0) out[[nm]] <- genes
      }
      out
    }

    # Combine ticked built-in sets with parsed custom sets into one named list.
    collect_sets <- shiny::reactive({
      builtin <- cellcycle_example_sets()
      chosen  <- if (is.null(input$builtin)) list() else builtin[input$builtin]
      utils::modifyList(chosen, parse_custom_sets(input$custom))
    })

    # (a) Cell-cycle scoring ---------------------------------------------------
    shiny::observeEvent(input$run_cc, {
      shiny::req(rv$obj)
      if (!require_pkgs("Seurat", "Cell-cycle scoring")) return(NULL)
      obj <- with_progress_notify({
        sc_cellcycle(rv$obj)
      }, message = "Scoring cell cycle (S / G2M / Phase)...")
      if (is.null(obj)) return(NULL)
      rv$obj <- obj
      res$cc_done <- TRUE
      mark_done(rv, "cellcycle")
      log_step(log_rv, "Cell-cycle scoring",
               params = list(),
               code = "obj <- Seurat::CellCycleScoring(obj, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)")
      shiny::showNotification(
        i18n("Cell-cycle scoring done: S.Score, G2M.Score and Phase added.",
             "细胞周期评分完成：已添加 S.Score、G2M.Score 和 Phase。"),
        type = "message")
    })

    # (b) Signature scoring ----------------------------------------------------
    shiny::observeEvent(input$run_sig, {
      shiny::req(rv$obj)
      sets <- collect_sets()
      if (length(sets) == 0) {
        shiny::showNotification(
          i18n("Select a built-in set or enter a custom gene set first.",
               "请先选择内置基因集或输入自定义基因集。"),
          type = "warning")
        return(NULL)
      }
      pkgs <- if (identical(input$method, "UCell")) c("Seurat", "UCell") else "Seurat"
      if (!require_pkgs(pkgs, "Signature scoring")) return(NULL)
      method <- input$method
      obj <- with_progress_notify({
        sc_modulescore(rv$obj, features = sets, method = method)
      }, message = "Scoring gene signatures...")
      if (is.null(obj)) return(NULL)
      rv$obj <- obj
      res$sig_done <- TRUE
      res$sets     <- names(sets)
      res$method   <- method
      mark_done(rv, "cellcycle")
      log_step(log_rv, "Signature scoring",
               params = list(method = method, sets = names(sets)),
               code = sprintf(
                 "obj <- sc_modulescore(obj, features = list(%s), method = '%s')",
                 paste(names(sets), collapse = ", "), method))
      shiny::showNotification(
        i18n(sprintf("Scored %d signature set(s).", length(sets)),
             sprintf("已评分 %d 个信号基因集。", length(sets))),
        type = "message")
    })

    # Preview feature choices: Phase (if present) + numeric score columns.
    output$preview_ui <- shiny::renderUI({
      cols <- obj_meta_cols(rv$obj)
      md   <- obj_meta(rv$obj)
      score_cols <- cols[grepl("Score|score|_UCell$|Phase", cols)]
      score_cols <- unique(score_cols)
      if (length(score_cols) == 0) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("Run scoring to enable the preview.",
                               "运行评分以启用预览。")))
      }
      shiny::selectInput(session$ns("feature"), NULL, choices = score_cols)
    })

    output$summary <- shiny::renderUI({
      if (!isTRUE(res$cc_done) && !isTRUE(res$sig_done)) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("Score the cell cycle and/or gene signatures.",
                               "对细胞周期和/或基因信号进行评分。")))
      }
      bslib::layout_columns(
        col_widths = c(4, 4, 4),
        stat_tile(i18n("Cell cycle", "细胞周期"),
                  if (isTRUE(res$cc_done)) i18n("Scored", "已评分") else i18n("No", "否")),
        stat_tile(i18n("Signatures", "信号集"),
                  if (isTRUE(res$sig_done)) length(res$sets) else 0L),
        stat_tile(i18n("Method", "方法"), res$method %||% "-")
      )
    })

    output$preview <- render_scop_plot(function() {
      obj <- rv$obj
      shiny::req(obj)
      feature <- input$feature
      shiny::req(feature)
      if (!require_pkgs("scop", "Score preview")) return(NULL)
      tryCatch({
        if (identical(feature, "Phase")) {
          # Categorical: colour the embedding by cell-cycle phase.
          sc_dimplot(obj, group_by = "Phase")
        } else {
          # Continuous: paint the score onto the embedding.
          sc_featureplot(obj, features = feature)
        }
      }, error = function(e) {
        shiny::showNotification(
          i18n(paste("Preview unavailable:", conditionMessage(e)),
               paste("预览不可用：", conditionMessage(e))),
          type = "error", duration = 10)
        NULL
      })
    })
  })
}
