#' Module: Cell-type annotation
#'
#' Give each cluster a biological identity. You can label clusters by hand using
#' the marker genes from the previous step, or predict labels automatically with
#' a reference (SingleR per cluster / Azimuth per cell). The result is stored in
#' a `celltype` column.
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
            reliable but needs expertise. <b>SingleR</b> labels each cluster from
            its pooled profile against a celldex reference (pick one matching
            your species and tissue). <b>Azimuth</b> maps each cell onto a
            curated human reference (PBMC, bone marrow, lung, ...). Both
            references are downloaded on first use and can be wrong for
            tissues they do not cover.",
      zh = "<b>手动</b>方式结合你的标志基因和先验知识——最可靠但需要专业经验。<b>SingleR</b> 用每个簇的汇总表达谱与 celldex 参考比对来为簇命名（请选择与物种和组织匹配的参考）。<b>Azimuth</b> 把每个细胞映射到精选的人类参考图谱（PBMC、骨髓、肺……）。两种参考首次使用时需下载，且对参考未覆盖的组织可能出错。"),
    read = list(
      en = "Left: the map coloured by cell type. Right: one bar per cluster.
            SingleR shows the delta (how far the best label beats the next);
            red = pruned as ambiguous, labelled Unassigned. Azimuth shows the
            mean per-cell prediction score; red = below 0.5. Manual labels
            show cluster sizes.",
      zh = "左：按细胞类型着色的嵌入图。右：每个簇一根条。SingleR 显示 delta（最佳标签领先次佳标签的幅度）；红色 = 因不明确被剔除，标为 Unassigned。Azimuth 显示每细胞预测分数的均值；红色 = 低于 0.5。手动标签显示各簇大小。"),
    example = list(
      en = "A cluster whose top markers are <code>CD3D</code>/<code>CD8A</code>
               would be labelled a 'CD8 T cell'. The Cell types tab adds a
               plain-language note for common labels.",
      zh = "如果某个簇的主要标志基因是 <code>CD3D</code>/<code>CD8A</code>，就会被标注为“CD8 T 细胞”。“细胞类型”页签会为常见标签附上通俗说明。")
  )
  controls <- shiny::tagList(
    label_with_help("Method",
                    "Manual = you name each cluster. SingleR/Azimuth = automatic reference-based prediction (needs internet the first time).",
                    "方法",
                    "手动 = 你为每个簇命名。SingleR/Azimuth = 基于参考数据集的自动预测（首次需要联网）。"),
    shiny::selectInput(ns("method"), NULL,
                       choices = c("Manual (marker-based)" = "manual",
                                   "SingleR (reference)"   = "singler",
                                   "Azimuth (reference)"   = "azimuth"),
                       selected = "manual"),
    shiny::conditionalPanel(
      sprintf("input['%s'] == 'manual'", ns("method")),
      label_with_help("Label each cluster",
                      "Type a cell-type name for every cluster of the active clustering, then click Apply labels. Blank = keep the cluster number.",
                      "为每个簇打标签",
                      "为当前聚类的每个簇输入一个细胞类型名称，然后点击“应用标签”。留空 = 保留簇编号。"),
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
                      "A labelled dataset to compare your cells against. Pick one that matches your species and tissue.",
                      "参考数据集",
                      "用于与你的细胞比对的带标签数据集。请选择与你的物种和组织相匹配的一个。"),
      shiny::selectInput(ns("ref"), NULL,
                         choices = c("Human Primary Cell Atlas (human)" = "HumanPrimaryCellAtlasData",
                                     "Blueprint/ENCODE (human)"         = "BlueprintEncodeData",
                                     "Monaco immune (human)"            = "MonacoImmuneData",
                                     "DICE immune (human)"              = "DatabaseImmuneCellExpressionData",
                                     "Novershtern haematopoietic (human)" = "NovershternHematopoieticData",
                                     "ImmGen (mouse)"                   = "ImmGenData",
                                     "Mouse RNA-seq (mouse)"            = "MouseRNAseqData"),
                         selected = "HumanPrimaryCellAtlasData"),
      label_with_help("Label level", "main = broad types (e.g. T cells); fine = subtypes (e.g. CD8 effector memory).",
                      "标签层级", "main = 大类（如 T 细胞）；fine = 亚型（如 CD8 效应记忆）。"),
      shiny::radioButtons(ns("level"), NULL, c("main" = "main", "fine" = "fine"), inline = TRUE),
      run_button(ns("run_singler"), "Run SingleR", "运行 SingleR")
    ),
    shiny::conditionalPanel(
      sprintf("input['%s'] == 'azimuth'", ns("method")),
      shiny::div(class = "omicone-note",
                 i18n("Azimuth maps human cells onto a tissue-specific reference (needs internet and the Azimuth package). Pick the reference for your tissue.",
                      "Azimuth 将人类细胞映射到组织特异的参考图谱（需要联网和 Azimuth 包）。请选择与你的组织匹配的参考。")),
      label_with_help("Reference", "Azimuth reference atlas; results are only meaningful for the tissue it covers.",
                      "参考图谱", "Azimuth 参考图谱；只有对其覆盖的组织结果才有意义。"),
      shiny::selectInput(ns("az_ref"), NULL,
                         choices = c("PBMC" = "pbmcref", "Bone marrow" = "bonemarrowref",
                                     "Tonsil" = "tonsilref", "Lung" = "lungref",
                                     "Kidney" = "kidneyref", "Heart" = "heartref",
                                     "Pancreas" = "pancreasref", "Adipose" = "adiposeref",
                                     "Fetal development" = "fetusref",
                                     "Human motor cortex" = "humancortexref"),
                         selected = "pbmcref"),
      label_with_help("Annotation level", "1 = coarse, 2 = intermediate, 3 = fine. References without numbered levels use their first annotation.",
                      "注释层级", "1 = 粗，2 = 中，3 = 细。没有编号层级的参考使用其第一个注释列。"),
      shiny::radioButtons(ns("az_level"), NULL, c("1" = "1", "2" = "2", "3" = "3"),
                          selected = "2", inline = TRUE),
      run_button(ns("run_azimuth"), "Run Azimuth", "运行 Azimuth")
    )
  )
  step_container(
    title     = list(en = "Cell-type annotation", zh = "细胞类型注释"),
    subtitle  = list(en = "Name each cluster by hand from its markers, or predict labels with SingleR / Azimuth.",
                     zh = "根据标志基因手动为簇命名，或用 SingleR / Azimuth 预测标签。"),
    explainer = explainer,
    controls  = controls,
    summary   = shiny::uiOutput(ns("summary")),
    preview   = shiny::tagList(
      shiny::uiOutput(ns("insight")),
      bslib::navset_card_tab(
        bslib::nav_panel(i18n("Annotation", "注释"),
          preview_plot_ui(ns("preview"), download = TRUE,
            guide = list(en = "The labelled map and per-cluster confidence will be drawn here.",
                         zh = "运行后，这里将绘制标注后的嵌入图与每个簇的置信度。"),
            caption = list(en = "Left: cells coloured by label. Right: one bar per cluster (SingleR delta, Azimuth mean score, or cluster size); red = low confidence.",
                           zh = "左：按标签着色的细胞。右：每个簇一根条（SingleR delta、Azimuth 平均分数或簇大小）；红色 = 低置信度。"))),
        bslib::nav_panel(i18n("Cell types", "细胞类型"),
          shiny::div(class = "omicone-table",
                     if (has_pkg("DT")) DT::dataTableOutput(ns("tbl"))
                     else shiny::verbatimTextOutput(ns("tbl"))))
      )
    )
  )
}

#' @rdname mod_annotate
#' @keywords internal
mod_annotate_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    res <- step_results(rv, "sc", tab = NULL, method = NULL, labels = NULL,
                        labels_sig = NULL, species_warn = NULL)

    # Cluster levels of the active clustering (Idents).
    cluster_levels <- shiny::reactive({
      obj <- rv$obj
      shiny::req(obj)
      lv <- tryCatch(levels(SeuratObject::Idents(obj)), error = function(e) NULL)
      if (is.null(lv) || !length(lv)) {
        md <- obj_meta(obj)
        if ("seurat_clusters" %in% colnames(md)) lv <- levels(factor(md$seurat_clusters))
      }
      lv
    })

    current_idents <- function(obj) {
      id <- tryCatch(SeuratObject::Idents(obj), error = function(e) NULL)
      if (is.null(id)) id <- factor(obj_meta(obj)$seurat_clusters)
      id
    }

    # Which clustering a set of manual labels belongs to: stored labels are
    # only filled back in while the clustering is the same one.
    clustering_sig <- function(obj) {
      id <- current_idents(obj)
      paste(obj_misc(obj, "omicone_cluster_col") %||% "", length(id),
            paste(as.integer(table(id)), collapse = ","))
    }

    # One text box per cluster; ids use the position, so cluster names with
    # spaces or symbols still make valid input ids. Values typed earlier (only
    # while the clusters are unchanged, since ids are positional) or the
    # labels last applied (matched by cluster name) are filled back in.
    last_sig <- NULL
    output$manual_inputs <- shiny::renderUI({
      lv <- cluster_levels()
      if (is.null(lv) || !length(lv)) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("No clusters found. Run clustering first.",
                               "未找到簇。请先运行聚类。")))
      }
      sig <- clustering_sig(shiny::isolate(rv$obj))
      prev <- shiny::isolate(res$labels)
      if (!identical(shiny::isolate(res$labels_sig), sig)) prev <- NULL
      same <- identical(sig, last_sig)
      last_sig <<- sig
      shiny::tagList(lapply(seq_along(lv), function(i) {
        cl <- lv[i]
        typed <- if (same) shiny::isolate(input[[paste0("lab_", i)]]) else NULL
        val <- if (!is.null(typed) && nzchar(typed)) typed
               else if (!is.null(prev) && cl %in% names(prev) && prev[[cl]] != cl) prev[[cl]]
               else ""
        shiny::textInput(ns(paste0("lab_", i)),
                         label = i18n(paste0("Cluster ", cl), paste0("簇 ", cl)),
                         value = val)
      }))
    })

    species_check <- function(ref_species) {
      sp <- guess_species(rv$obj)
      if (identical(sp, ref_species)) return(NULL)
      msg <- i18n(sprintf("The data look %s but the reference is %s: labels will be unreliable.",
                          sp, ref_species),
                  sprintf("数据看起来是%s，但参考是%s：标签将不可靠。",
                          if (sp == "human") "人类" else "小鼠",
                          if (ref_species == "human") "人类" else "小鼠"))
      shiny::showNotification(msg, type = "warning", duration = 15)
      sprintf("%s data, %s reference", sp, ref_species)
    }

    # ---- Manual ----
    shiny::observeEvent(input$apply_manual, {
      shiny::req(rv$obj)
      lv <- cluster_levels()
      shiny::req(lv)
      labels <- vapply(seq_along(lv), function(i) {
        v <- input[[paste0("lab_", i)]]
        if (is.null(v) || !nzchar(trimws(v))) as.character(lv[i]) else trimws(v)
      }, character(1))
      names(labels) <- as.character(lv)
      obj <- rv$obj
      idents <- as.character(current_idents(obj))
      obj$celltype <- unname(labels[idents])
      rv$obj <- obj
      res$labels <- labels
      res$labels_sig <- clustering_sig(obj)
      res$method <- "manual"
      res$species_warn <- NULL
      res$tab <- data.frame(cluster = names(labels), label = unname(labels),
                            raw_label = unname(labels), confidence = NA_real_,
                            low_conf = FALSE,
                            n = as.integer(table(factor(idents, levels = lv))),
                            stringsAsFactors = FALSE)
      mark_done(rv, "annotate")
      log_step(log_rv, "Annotate", params = list(method = "manual", labels = as.list(labels)),
               code = annotate_manual_log_code(labels))
      shiny::showNotification(i18n("Applied manual cell-type labels.", "已应用手动细胞类型标签。"),
                              type = "message")
    })

    # ---- SingleR ----
    shiny::observeEvent(input$run_singler, {
      shiny::req(rv$obj)
      if (!require_pkgs(c("Seurat", "SingleR", "celldex"), "SingleR annotation")) return(NULL)
      ref_name <- input$ref
      level <- input$level
      warn <- species_check(singler_ref_species(ref_name))
      out <- with_progress_notify({
        ref <- getExportedValue("celldex", ref_name)()
        labels <- SummarizedExperiment::colData(ref)[[paste0("label.", level)]]
        annotate_singler(rv$obj, ref, labels, current_idents(rv$obj))
      }, message = "Running SingleR (may download a reference)...")
      if (is.null(out)) return(NULL)
      rv$obj <- out$obj
      res$tab <- out$table
      res$method <- "singler"
      res$labels <- NULL
      res$species_warn <- warn
      mark_done(rv, "annotate")
      log_step(log_rv, "Annotate",
               params = list(method = "SingleR", reference = ref_name, level = level),
               code = annotate_singler_log_code(ref_name, level))
      shiny::showNotification(i18n("SingleR annotation done.", "SingleR 注释完成。"),
                              type = "message")
    })

    # ---- Azimuth ----
    shiny::observeEvent(input$run_azimuth, {
      shiny::req(rv$obj)
      if (!require_pkgs(c("Seurat", "Azimuth"), "Azimuth annotation")) return(NULL)
      reference <- input$az_ref
      level <- input$az_level
      warn <- species_check("human")
      assay0 <- obj_default_assay(rv$obj) %||% "RNA"
      out <- with_progress_notify({
        annotate_azimuth(rv$obj, reference = reference, level = level,
                         clusters = current_idents(rv$obj))
      }, message = "Running Azimuth (needs internet)...")
      if (is.null(out)) return(NULL)
      rv$obj <- out$obj
      res$tab <- out$table
      res$method <- "azimuth"
      res$labels <- NULL
      res$species_warn <- warn
      mark_done(rv, "annotate")
      log_step(log_rv, "Annotate",
               params = list(method = "Azimuth", reference = reference, level = level,
                             column = out$col),
               code = annotate_azimuth_log_code(reference, out$col, assay0))
      if (!grepl(paste0("(\\.l|_level_)", level, "$"), out$col)) {
        shiny::showNotification(
          i18n(sprintf("This reference has no level %s; used %s.", level, out$col),
               sprintf("该参考没有第 %s 层级；已使用 %s。", level, out$col)),
          type = "warning", duration = 12)
      }
      shiny::showNotification(i18n("Azimuth annotation done.", "Azimuth 注释完成。"),
                              type = "message")
    })

    # The per-cluster table: from the last run, or derived from an existing
    # celltype column (e.g. an imported, already annotated object).
    cluster_tab <- shiny::reactive({
      tab <- res$tab
      if (!is.null(tab)) return(tab)
      obj <- rv$obj
      md <- obj_meta(obj)
      shiny::req("celltype" %in% colnames(md))
      annotation_cluster_table(current_idents(obj), md$celltype)
    })

    output$summary <- shiny::renderUI({
      md <- obj_meta(rv$obj)
      if (!("celltype" %in% colnames(md))) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("No annotation yet. Pick a method and run it.",
                               "还没有注释结果。请选择一种方法并运行。")))
      }
      tab <- res$tab
      shiny::tagList(
        stat_tile(i18n("Cell types", "细胞类型数"),
                  format(length(unique(md$celltype)), big.mark = ",")),
        stat_tile(i18n("Cells annotated", "已注释细胞数"),
                  format(sum(!is.na(md$celltype) & md$celltype != "Unassigned"), big.mark = ",")),
        stat_tile(i18n("Method", "方法"), res$method %||% i18n("existing column", "已有列")),
        if (!is.null(tab)) {
          stat_tile(i18n("Low-confidence clusters", "低置信度簇"),
                    sprintf("%d of %d", sum(tab$low_conf), nrow(tab)))
        }
      )
    })

    output$insight <- shiny::renderUI({
      tab <- res$tab
      if (is.null(tab)) return(NULL)
      low <- sum(tab$low_conf)
      warn_en <- if (!is.null(res$species_warn)) sprintf(" Species mismatch (%s).", res$species_warn) else ""
      warn_zh <- if (!is.null(res$species_warn)) sprintf("物种不匹配（%s）。", res$species_warn) else ""
      if (identical(res$method, "manual")) {
        return(insight_bar(
          sprintf("%d clusters labelled by hand into %d cell types.", nrow(tab),
                  length(unique(tab$label))),
          sprintf("已手动将 %d 个簇标注为 %d 种细胞类型。", nrow(tab), length(unique(tab$label)))))
      }
      what_en <- if (identical(res$method, "singler")) "pruned as ambiguous (Unassigned)"
                 else "with mean prediction score below 0.5"
      what_zh <- if (identical(res$method, "singler")) "因不明确被剔除（Unassigned）"
                 else "平均预测分数低于 0.5"
      insight_bar(
        sprintf("%d clusters -> %d labels; %d cluster(s) %s. Check those clusters' markers before trusting their names.%s",
                nrow(tab), length(unique(tab$label)), low, what_en, warn_en),
        sprintf("%d 个簇 -> %d 种标签；%d 个簇%s。采信这些簇的名称前请核对其标志基因。%s",
                nrow(tab), length(unique(tab$label)), low, what_zh, warn_zh))
    })

    render_step_plot(output, input, "preview", function() {
      obj <- rv$obj
      shiny::req(obj, "celltype" %in% obj_meta_cols(obj))
      annotate_plot(obj, cluster_tab(), res$method %||% "manual")
    }, name = "annotate", width = 12, height = 7)

    output$tbl <- render_tbl_wrap(function() {
      tab <- cluster_tab()
      shiny::req(tab)
      out <- tab[, intersect(c("cluster", "label", "raw_label", "confidence", "low_conf", "n",
                               "share"), names(tab)), drop = FALSE]
      out$about <- explain_celltype(out$label)
      out
    })
  })
}
