#' Module: RNA velocity
#'
#' Estimate each cell's likely next state from the ratio of unspliced to
#' spliced mRNA, giving a directional field on the UMAP. Wraps
#' scop::RunSCVELO() (scVelo, Python) via sc_velocity().
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
      en = "Predict where each cell is heading by comparing its unspliced (nascent)
            and spliced (mature) mRNA.",
      zh = "通过比较每个细胞的未剪接（新生）与已剪接（成熟）mRNA，预测细胞的走向。"),
    why  = list(
      en = "Pseudotime gives an order but not a direction. Velocity adds an arrow per
            cell, an independent check on which way differentiation runs.",
      zh = "拟时序给出顺序但不给出方向。速率为每个细胞加上一个箭头，可独立检验分化的走向。"),
    how  = list(
      en = "The object needs <b>spliced</b> and <b>unspliced</b> assays (velocyto or
            kallisto|bustools), a PCA and a UMAP, plus a Python environment
            (scop::PrepareEnv()). <b>Stochastic</b> (default) is fast and stable;
            <b>dynamical</b> fits full splicing kinetics and is much slower.",
      zh = "对象需包含 <b>spliced</b> 与 <b>unspliced</b> 两个 assay（来自 velocyto 或 kallisto|bustools）、PCA 与 UMAP，以及 Python 环境（scop::PrepareEnv()）。<b>stochastic</b>（默认）快速稳定；<b>dynamical</b> 拟合完整剪接动力学，速度慢得多。"),
    read = list(
      en = "Streamlines on the UMAP follow the projected velocity; cells are
            coloured by group. Long, coherent streams mean a strong directional
            flow; tangled or absent streams mean the signal is weak there.",
      zh = "UMAP 上的流线沿投影后的速率方向；细胞按分组着色。长而一致的流线＝强的定向流；杂乱或缺失的流线＝该区域信号弱。"),
    example = list(
      en = "In pancreatic endocrinogenesis, streams run from Ngn3+ progenitors toward
            alpha, beta and delta cells.",
      zh = "在胰腺内分泌发育中，流线从 Ngn3+ 祖细胞流向 alpha、beta 与 delta 细胞。")
  )
  controls <- shiny::tagList(
    shiny::div(class = "omicone-note",
               i18n("Requires spliced/unspliced assays, PCA + UMAP, and a Python conda env (scop::PrepareEnv()).",
                    "需要 spliced/unspliced assay、PCA 与 UMAP，以及 Python conda 环境（scop::PrepareEnv()）。")),
    label_with_help("Mode",
                    "stochastic = fast, stable default; deterministic = steady-state model; dynamical = full kinetics, slow.",
                    label_zh = "模式",
                    tip_zh = "stochastic = 快速稳定的默认；deterministic = 稳态模型；dynamical = 完整动力学，较慢。"),
    shiny::selectInput(ns("mode"), NULL,
                       choices = c("Stochastic" = "stochastic",
                                   "Deterministic" = "deterministic",
                                   "Dynamical" = "dynamical"),
                       selected = "stochastic"),
    label_with_help("Group by (metadata column)",
                    "Cell grouping used to colour the velocity plot (e.g. seurat_clusters, celltype).",
                    label_zh = "分组依据（元数据列）",
                    tip_zh = "用于为速率图着色的细胞分组（例如 seurat_clusters、celltype）。"),
    shiny::selectInput(ns("group_by"), NULL, choices = NULL),
    run_button(ns("run"), "Run velocity", "运行 RNA 速率")
  )
  step_container(id = id, 
    title     = list(en = "RNA velocity", zh = "RNA 速率"),
    subtitle  = list(en = "Each cell's future direction from splicing kinetics.",
                     zh = "基于剪接动力学推断每个细胞的未来方向。"),
    explainer = explainer,
    controls  = controls,
    summary   = shiny::uiOutput(ns("summary")),
    preview   = shiny::tagList(
      shiny::uiOutput(ns("insight")),
      preview_plot_ui(ns("preview"), download = TRUE,
        guide = list(en = "The velocity stream plot will be drawn here.",
                     zh = "运行后，这里将绘制速率流线图。"),
        caption = list(en = "RNA-velocity streamlines on the UMAP; cells coloured by group.",
                       zh = "UMAP 上的 RNA 速率流线；细胞按分组着色。")))
  )
}

#' @rdname mod_velocity
#' @keywords internal
mod_velocity_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    res <- step_results(rv, "sc", done = FALSE, mode = NULL, group_by = NULL)

    shiny::observe({
      obj <- rv$obj
      cols <- categorical_cols(obj_meta(obj))
      def <- default_group_col(cols, obj_misc(obj, "omicone_cluster_col"))
      shiny::updateSelectInput(session, "group_by", choices = cols,
                               selected = keep_selected(shiny::isolate(input$group_by),
                                                        cols, def))
    })

    shiny::observeEvent(input$run, {
      shiny::req(rv$obj)
      mode <- input$mode
      group_by <- input$group_by
      shiny::req(mode, group_by)
      miss <- velocity_missing(obj_assays(rv$obj), obj_reductions(rv$obj))
      if (length(miss)) {
        shiny::showNotification(
          i18n(sprintf("RNA velocity needs: %s. Spliced/unspliced counts come from velocyto or kallisto|bustools.",
                       paste(miss, collapse = ", ")),
               sprintf("RNA 速率缺少：%s。剪接/未剪接计数来自 velocyto 或 kallisto|bustools。",
                       paste(miss, collapse = ", "))),
          type = "warning", duration = 12)
        return(NULL)
      }
      if (!require_pkgs("scop", "RNA velocity")) return(NULL)
      shiny::showNotification(
        i18n("RNA velocity runs in Python. If it fails, run scop::PrepareEnv().",
             "RNA 速率在 Python 中运行。如果失败，请运行 scop::PrepareEnv()。"),
        type = "warning", duration = 8)
      obj <- with_progress_notify({
        sc_velocity(rv$obj, group_by = group_by, mode = mode)
      }, message = sprintf("Estimating RNA velocity (%s)...", mode))
      if (is.null(obj)) return(NULL)
      rv$obj <- obj
      res$done     <- TRUE
      res$mode     <- mode
      res$group_by <- group_by
      mark_done(rv, "velocity")
      log_step(log_rv, "RNA velocity",
               params = list(mode = mode, group_by = group_by,
                             linear_reduction = "pca", nonlinear_reduction = "umap"),
               code = velocity_log_code(group_by, mode))
      shiny::showNotification(
        i18n(sprintf("RNA velocity (%s) finished on '%s'.", mode, group_by),
             sprintf("RNA 速率（%s）已完成，分组：%s。", mode, group_by)),
        type = "message")
    })

    output$summary <- shiny::renderUI({
      if (!isTRUE(res$done)) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("Check the spliced/unspliced assays exist, then click <b>Run velocity</b>.",
                               "确认存在 spliced/unspliced assay 后，点击<b>运行 RNA 速率</b>。")))
      }
      shiny::tagList(
        stat_tile(i18n("Mode", "模式"), res$mode),
        stat_tile(i18n("Group by", "分组依据"), res$group_by),
        stat_tile(i18n("Embedding", "嵌入"), "umap")
      )
    })

    output$insight <- shiny::renderUI({
      if (!isTRUE(res$done)) return(NULL)
      insight_bar(
        sprintf("scVelo <b>%s</b> model, neighbours on the PCA, velocity projected onto the UMAP. Projection onto 2-D can bend or invent directions; confirm a flow with marker genes or pseudotime before relying on it.",
                res$mode),
        sprintf("scVelo <b>%s</b> 模型，近邻基于 PCA，速率投影到 UMAP。二维投影可能扭曲甚至凭空产生方向；依赖某一流向之前，请用标志基因或拟时序加以确认。",
                res$mode))
    })

    render_step_plot(output, input, "preview", function() {
      shiny::req(res$done)
      sc_velocityplot(rv$obj, mode = res$mode, group_by = res$group_by)
    }, name = "velocity")
  })
}
