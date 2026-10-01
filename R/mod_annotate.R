#' Module: Cell-type annotation
#'
#' Give each cluster a biological identity. You can label clusters by hand using
#' the marker genes from the previous step, or predict labels automatically with
#' a reference (SingleR / Azimuth). The result is stored in a `celltype` column.
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_annotate
NULL

#' @rdname mod_annotate
#' @keywords internal
mod_annotate_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = list(en = "Cell-type annotation", zh = "细胞类型注释"),
    what = list(
      en = "Turn anonymous clusters into named cell types (e.g. 'CD8 T cell').",
      zh = "把匿名的簇转变为有名字的细胞类型（例如“CD8 T 细胞”）。"),
    why  = list(
      en = "Numbered clusters mean nothing biologically. Annotation is what lets
            you talk about your data in terms of real cell populations.",
      zh = "带编号的簇在生物学上没有意义。注释能让你用真实的细胞群体来描述数据。"),
    how  = list(
      en = "<b>Manual</b> uses your marker genes plus prior knowledge -- most
            reliable but needs expertise. <b>SingleR</b> and <b>Azimuth</b> compare
            each cell to a labelled reference; they are fast but need internet to
            download the reference and can be wrong for unusual tissues.",
      zh = "<b>手动</b>方式结合你的标志基因和先验知识——最可靠但需要专业经验。<b>SingleR</b> 和 <b>Azimuth</b> 会将每个细胞与带标签的参考数据集比对；它们速度快，但需要联网下载参考数据集，且对于不常见的组织可能出错。"),
    read = list(
      en = "Each cluster gets the reference label with the highest score. A
            low-confidence label means no good match exists — check that
            cluster's markers by hand before trusting the name.",
      zh = "每个簇获得得分最高的参考标签。低置信度标签意味着没有良好匹配——采信该名称前，请手动核对这个簇的标志基因。"),
    example = list(
      en = "A cluster whose top markers are <code>CD3D</code>/<code>CD8A</code>
               would be labelled a 'CD8 T cell'. Hover a point in the plot for a
               plain-language explanation of each cell type.",
      zh = "如果某个簇的主要标志基因是 <code>CD3D</code>/<code>CD8A</code>，就会被标注为“CD8 T 细胞”。将鼠标悬停在图中的点上，可查看每种细胞类型的通俗说明。")
  )
  controls <- shiny::tagList(
    label_with_help("Method",
                    "Manual = you name each cluster. SingleR/Azimuth = automatic reference-based prediction (needs internet).",
                    "方法",
                    "手动 = 你为每个簇命名。SingleR/Azimuth = 基于参考数据集的自动预测（需要联网）。"),
    shiny::selectInput(ns("method"), NULL,
                       choices = c("Manual (marker-based)" = "manual",
                                   "SingleR (reference)"   = "singler",
                                   "Azimuth (reference)"   = "azimuth"),
                       selected = "manual"),
    shiny::conditionalPanel(
      sprintf("input['%s'] == 'manual'", ns("method")),
      label_with_help("Label each cluster",
                      "Type a cell-type name for every cluster, then click Apply labels.",
                      "为每个簇打标签",
                      "为每个簇输入一个细胞类型名称，然后点击“应用标签”。"),
      shiny::uiOutput(ns("manual_inputs")),
      shiny::actionButton(ns("apply_manual"),
                          i18n("Apply labels", "应用标签"),
                          icon = shiny::icon("check"),
                          class = "btn-primary omicone-run w-100")
    ),
    shiny::conditionalPanel(
      sprintf("input['%s'] == 'singler'", ns("method")),
      shiny::div(class = "omicone-note",
                 i18n("SingleR downloads a celldex reference (needs internet) the first time.",
                      "首次运行时 SingleR 会下载一个 celldex 参考数据集（需要联网）。")),
      label_with_help("Reference",
                      "A labelled dataset to compare your cells against. Pick one that matches your tissue.",
                      "参考数据集",
                      "用于与你的细胞比对的带标签数据集。请选择与你的组织相匹配的一个。"),
      shiny::selectInput(ns("ref"), NULL,
                         choices = c("Human Primary Cell Atlas" = "HumanPrimaryCellAtlasData",
                                     "Blueprint/ENCODE"         = "BlueprintEncodeData",
                                     "Monaco immune"            = "MonacoImmuneData",
                                     "Mouse RNA-seq (ImmGen)"   = "ImmGenData"),
                         selected = "HumanPrimaryCellAtlasData"),
      run_button(ns("run_singler"), "Run SingleR", "运行 SingleR")
    ),
    shiny::conditionalPanel(
      sprintf("input['%s'] == 'azimuth'", ns("method")),
      shiny::div(class = "omicone-note",
                 i18n("Azimuth maps to a curated reference (needs internet and the Azimuth package).",
                      "Azimuth 会映射到一个精选的参考数据集（需要联网和 Azimuth 包）。")),
      run_button(ns("run_azimuth"), "Run Azimuth", "运行 Azimuth")
    )
  )
  step_container(
    title     = list(en = "Cell-type annotation", zh = "细胞类型注释"),
    subtitle  = list(en = "Name the clusters with reference databases.",
                     zh = "借助参考数据库为簇命名。"),
    explainer = explainer,
    controls  = controls,
    summary   = shiny::uiOutput(ns("summary")),
    preview   = preview_plot_ui(ns("preview"),
      guide = list(en = "Annotation scores and the labelled embedding will be drawn here.",
                   zh = "运行后，这里将绘制注释得分与标注后的嵌入图。"),
      caption = list(en = "Clusters labelled by their best reference match.",
                     zh = "按最佳参考匹配标注的簇。"))
  )
}

#' @rdname mod_annotate
#' @keywords internal
mod_annotate_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Cluster levels from Idents (or seurat_clusters) of the working object.
    cluster_levels <- shiny::reactive({
      obj <- rv$obj
      shiny::req(obj)
      lv <- tryCatch(levels(Seurat::Idents(obj)), error = function(e) NULL)
      if (is.null(lv) || !length(lv)) {
        md <- obj_meta(obj)
        if ("seurat_clusters" %in% colnames(md)) lv <- levels(factor(md$seurat_clusters))
      }
      lv
    })

    # One text box per cluster for manual labelling.
    output$manual_inputs <- shiny::renderUI({
      lv <- cluster_levels()
      if (is.null(lv) || !length(lv)) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("No clusters found. Run clustering first.",
                               "未找到簇。请先运行聚类。")))
      }
      shiny::tagList(lapply(lv, function(cl) {
        shiny::textInput(ns(paste0("lab_", cl)),
                         label = i18n(paste0("Cluster ", cl), paste0("簇 ", cl)),
                         value = "")
      }))
    })

    # ---- Manual ----
    shiny::observeEvent(input$apply_manual, {
      shiny::req(rv$obj)
      lv <- cluster_levels()
      shiny::req(lv)
      labels <- vapply(lv, function(cl) {
        v <- input[[paste0("lab_", cl)]]
        if (is.null(v) || !nzchar(trimws(v))) as.character(cl) else trimws(v)
      }, character(1))
      names(labels) <- as.character(lv)
      obj <- rv$obj
      idents <- tryCatch(as.character(Seurat::Idents(obj)), error = function(e) NULL)
      if (is.null(idents)) {
        md <- obj_meta(obj)
        idents <- as.character(md$seurat_clusters)
      }
      obj$celltype <- unname(labels[idents])
      rv$obj <- obj
      mark_done(rv, "annotate")
      log_step(log_rv, "Annotate (manual)",
               params = as.list(labels),
               code = c("labels <- c(  # cluster -> cell type",
                        paste0("  '", names(labels), "' = '", labels, "'",
                               c(rep(",", length(labels) - 1), ""), collapse = "\n"),
                        ")",
                        "obj$celltype <- labels[as.character(Seurat::Idents(obj))]"))
      shiny::showNotification("Applied manual cell-type labels.", type = "message")
    })

    # ---- SingleR ----
    shiny::observeEvent(input$run_singler, {
      shiny::req(rv$obj)
      if (!require_pkgs(c("Seurat", "SingleR", "celldex"), "SingleR annotation")) return(NULL)
      obj <- with_progress_notify({
        ref_se <- switch(input$ref,
          HumanPrimaryCellAtlasData = celldex::HumanPrimaryCellAtlasData(),
          BlueprintEncodeData       = celldex::BlueprintEncodeData(),
          MonacoImmuneData          = celldex::MonacoImmuneData(),
          ImmGenData                = celldex::ImmGenData())
        ref <- list(data = ref_se, labels = ref_se$label.main)
        o <- annotate_singler(rv$obj, ref)
        o$celltype <- o[["SingleR"]][, 1]
        o
      }, message = "Running SingleR (may download a reference)...")
      if (is.null(obj)) return(NULL)
      rv$obj <- obj
      mark_done(rv, "annotate")
      log_step(log_rv, "Annotate (SingleR)",
               params = list(reference = input$ref),
               code = c(sprintf('ref <- celldex::%s()', input$ref),
                        "sce <- Seurat::as.SingleCellExperiment(obj)",
                        "pred <- SingleR::SingleR(test=sce, ref=ref, labels=ref$label.main)",
                        "obj$celltype <- pred$labels"))
      shiny::showNotification("SingleR annotation done.", type = "message")
    })

    # ---- Azimuth ----
    shiny::observeEvent(input$run_azimuth, {
      shiny::req(rv$obj)
      if (!require_pkgs(c("Seurat", "Azimuth"), "Azimuth annotation")) return(NULL)
      obj <- with_progress_notify({
        o <- Azimuth::RunAzimuth(rv$obj, reference = "pbmcref")
        md <- obj_meta(o)
        pc <- grep("^predicted.celltype", colnames(md), value = TRUE)
        if (length(pc)) o$celltype <- md[[pc[1]]]
        o
      }, message = "Running Azimuth (needs internet)...")
      if (is.null(obj)) return(NULL)
      rv$obj <- obj
      mark_done(rv, "annotate")
      log_step(log_rv, "Annotate (Azimuth)",
               params = list(reference = "pbmcref"),
               code = c('obj <- Azimuth::RunAzimuth(obj, reference = "pbmcref")',
                        'obj$celltype <- obj$predicted.celltype.l2'))
      shiny::showNotification("Azimuth annotation done.", type = "message")
    })

    output$summary <- shiny::renderUI({
      md <- obj_meta(rv$obj)
      if (is.null(md) || !("celltype" %in% colnames(md))) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("No annotation yet. Pick a method and run it.",
                               "还没有注释结果。请选择一种方法并运行。")))
      }
      tab <- sort(table(md$celltype), decreasing = TRUE)
      df <- data.frame(celltype = names(tab), n = as.integer(tab),
                       stringsAsFactors = FALSE)
      shiny::tagList(
        bslib::layout_columns(
          col_widths = c(6, 6),
          stat_tile(i18n("Cell types", "细胞类型数"), format(nrow(df), big.mark = ",")),
          stat_tile(i18n("Cells annotated", "已注释细胞数"), format(sum(df$n), big.mark = ","))
        ),
        shiny::tags$table(
          class = "table table-sm omicone-dist",
          shiny::tags$thead(shiny::tags$tr(
            shiny::tags$th(i18n("Cell type", "细胞类型")),
            shiny::tags$th(i18n("Cells", "细胞数")))),
          shiny::tags$tbody(lapply(seq_len(nrow(df)), function(i) {
            shiny::tags$tr(
              shiny::tags$td(df$celltype[i]),
              shiny::tags$td(format(df$n[i], big.mark = ",")))
          }))
        )
      )
    })

    output$preview <- render_preview_plot(function() {
      obj <- rv$obj
      shiny::req(obj)
      md <- obj_meta(obj)
      shiny::req("celltype" %in% colnames(md))
      red <- if (has_reduction(obj, "umap")) "umap" else obj_reductions(obj)[1]
      shiny::req(!is.na(red), length(red) > 0)
      df <- embedding_df(obj, red, color_by = "celltype")
      counts <- table(df$color)
      df$n <- as.integer(counts[as.character(df$color)])
      df$explanation <- explain_celltype(df$color)
      df$text <- sprintf("cell type: %s\ncells of this type: %s\n%s",
                         df$color, format(df$n, big.mark = ","),
                         ifelse(nzchar(df$explanation), df$explanation,
                                "(no plain-language description available)"))
      cats <- length(unique(df$color))
      ggplot2::ggplot(df, ggplot2::aes(x = dim1, y = dim2, colour = color)) +
        ggplot2::geom_point(size = 0.5, alpha = 0.7) +
        ggplot2::scale_colour_manual(values = sc_palette(cats), name = "Cell type") +
        ggplot2::labs(x = paste0(red, " 1"), y = paste0(red, " 2"),
                      title = "Cells coloured by annotated cell type") +
        omicone_theme()
    })
  })
}
