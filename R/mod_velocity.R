#' Module: RNA velocity
#'
#' Estimate the future transcriptional state of each cell from the ratio of
#' spliced to unspliced mRNA, giving a directional "velocity" field on the
#' embedding. Wraps scop::RunSCVELO (scVelo, Python) via sc_velocity().
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_velocity
NULL

#' @rdname mod_velocity
#' @keywords internal
mod_velocity_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = list(en = "RNA velocity", zh = "RNA 速率"),
    what = list(
      en = "Predict where each cell is heading by comparing its unspliced
            (nascent) and spliced (mature) mRNA.",
      zh = "通过比较每个细胞的未剪接（新生）与已剪接（成熟）mRNA，预测细胞的走向。"),
    why  = list(
      en = "Pseudotime gives an order but not a direction. Velocity adds an arrow
            per cell, showing the likely direction of differentiation on the
            embedding.",
      zh = "拟时序给出顺序但不给出方向。速率为每个细胞添加一个箭头，在降维图上展示可能的分化方向。"),
    how  = list(
      en = "You need <b>spliced / unspliced layers</b> (from velocyto or
            kallisto|bustools) in the object, and a <b>Python conda environment</b>
            (run scop::PrepareEnv() once). The <b>dynamical</b> mode is most
            accurate; <b>stochastic</b> / <b>deterministic</b> are faster.",
      zh = "对象中需要包含<b>剪接 / 未剪接图层</b>（来自 velocyto 或 kallisto|bustools），并需要 <b>Python conda 环境</b>（首次使用请运行 scop::PrepareEnv()）。<b>dynamical</b> 模式最准确；<b>stochastic</b> / <b>deterministic</b> 更快。"),
    read = list(
      en = "Arrows point from a cell's current state toward its likely future
            state. Long, coherent arrows mean a strong directional flow; short,
            random ones mean the signal is weak in that region.",
      zh = "箭头从细胞当前状态指向其可能的未来状态。长而一致的箭头＝强的定向流；短而杂乱的箭头＝该区域信号弱。"),
    example = list(
      en = "In a differentiating system, velocity arrows should flow from
               progenitor cells outward toward the mature cell types.",
      zh = "在一个分化系统中，速率箭头应从祖细胞向外流向成熟细胞类型。")
  )
  controls <- shiny::tagList(
    shiny::div(class = "omicone-note",
               i18n("Requires spliced/unspliced layers and a Python conda env (scop::PrepareEnv()).",
                    "需要剪接/未剪接图层以及 Python conda 环境（scop::PrepareEnv()）。")),
    label_with_help("Mode",
                    "dynamical = most accurate (recommended); stochastic/deterministic are faster approximations.",
                    label_zh = "模式",
                    tip_zh = "dynamical = 最准确（推荐）；stochastic/deterministic 为更快的近似方法。"),
    shiny::selectInput(ns("mode"), NULL,
                       choices = c("Dynamical" = "dynamical",
                                   "Stochastic" = "stochastic",
                                   "Deterministic" = "deterministic"),
                       selected = "dynamical"),
    label_with_help("Group by (metadata column)",
                    "Cell grouping used to summarise and colour the velocity field (e.g. seurat_clusters, celltype).",
                    label_zh = "分组依据（元数据列）",
                    tip_zh = "用于汇总并为速率场着色的细胞分组（例如 seurat_clusters、celltype）。"),
    shiny::selectInput(ns("group_by"), NULL, choices = NULL),
    run_button(ns("run"), "Run velocity", "运行 RNA 速率")
  )
  step_container(
    title     = list(en = "RNA velocity", zh = "RNA 速率"),
    subtitle  = list(en = "Each cell's future direction from splicing kinetics.",
                     zh = "基于剪接动力学推断每个细胞的未来方向。"),
    explainer = explainer,
    controls  = controls,
    summary   = shiny::uiOutput(ns("summary")),
    preview   = preview_plot_ui(ns("preview"),
      guide = list(en = "The velocity field will be drawn here.",
                   zh = "运行后，这里将绘制速率场。"),
      caption = list(en = "RNA velocity field over the embedding.",
                     zh = "叠加在嵌入图上的 RNA 速率场。"))
  )
}

#' @rdname mod_velocity
#' @keywords internal
mod_velocity_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    res <- shiny::reactiveValues(done = FALSE, mode = NULL, group_by = NULL)

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
      # scop::RunSCVELO drives scVelo through a Python conda environment.
      if (!require_pkgs("scop", "RNA velocity")) return(NULL)
      shiny::showNotification(
        i18n("RNA velocity runs in Python and needs spliced/unspliced layers. If it fails, run scop::PrepareEnv().",
             "RNA 速率在 Python 中运行，且需要剪接/未剪接图层。如果失败，请运行 scop::PrepareEnv()。"),
        type = "warning", duration = 8)
      group_by <- input$group_by
      shiny::req(group_by)
      obj <- with_progress_notify({
        sc_velocity(rv$obj, group_by = group_by, mode = input$mode)
      }, message = sprintf("Estimating RNA velocity (%s)...", input$mode))
      if (is.null(obj)) return(NULL)
      rv$obj <- obj
      res$done     <- TRUE
      res$mode     <- input$mode
      res$group_by <- group_by
      mark_done(rv, "velocity")
      log_step(log_rv, "RNA velocity",
               params = list(mode = input$mode, group_by = group_by),
               code = sprintf(
                 'obj <- sc_velocity(obj, group_by="%s", mode="%s")',
                 group_by, input$mode))
      shiny::showNotification(
        i18n(sprintf("RNA velocity (%s) finished on '%s'.", input$mode, group_by),
             sprintf("RNA 速率（%s）已完成，分组：%s。", input$mode, group_by)),
        type = "message")
    })

    output$summary <- shiny::renderUI({
      if (!isTRUE(res$done)) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("Ensure spliced/unspliced layers exist, then click <b>Run velocity</b>.",
                               "确认存在剪接/未剪接图层后，点击<b>运行 RNA 速率</b>。")))
      }
      bslib::layout_columns(
        col_widths = c(6, 6),
        stat_tile(i18n("Mode", "模式"), res$mode),
        stat_tile(i18n("Group by", "分组依据"), res$group_by)
      )
    })

    output$preview <- render_scop_plot(function() {
      shiny::req(res$done)
      # Velocity stream/grid overlaid on the embedding.
      sc_velocityplot(rv$obj)
    })
  })
}
