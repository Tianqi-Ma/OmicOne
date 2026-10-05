#' Shared controls of the sample-level steps (pseudobulk DE, abundance)
#'
#' Both steps need the same three columns -- which cells belong to which
#' sample, the sample's condition, and the cell type -- and the same way out
#' when the condition is not in the metadata yet: a small sample sheet (one row
#' per sample) joined onto the cells. These helpers are called from inside the
#' two modules' UI / server functions (same namespace), not as modules.
#'
#' @name mod_sample_design
#' @keywords internal
NULL

#' @rdname mod_sample_design
#' @param ns The calling module's namespace function.
#' @param levels Also offer the reference / comparison levels.
#' @param covariate Also offer an optional sample-level covariate.
#' @keywords internal
sample_design_controls <- function(ns, levels = TRUE, covariate = FALSE) {
  shiny::tagList(
    label_with_help("Sample column",
                    "Which sample (patient, donor, library) each cell came from. Samples are the replicates of the test.",
                    "样本列",
                    "每个细胞来自哪个样本（患者、供体、文库）。样本是检验中的重复单位。"),
    shiny::uiOutput(ns("sample_ui")),
    label_with_help("Condition column",
                    "The sample-level variable to compare (e.g. response, treatment, timepoint). It must be the same for every cell of a sample.",
                    "条件列",
                    "要比较的样本层面变量（如疗效、处理、时间点）。同一样本的所有细胞必须相同。"),
    shiny::uiOutput(ns("cond_ui")),
    if (levels) shiny::uiOutput(ns("levels_ui")),
    if (covariate) shiny::tagList(
      label_with_help("Covariate (optional)",
                      "A second sample-level variable added to the model: a patient id for paired pre/post samples, or a batch. Leave at none if unsure.",
                      "协变量（可选）",
                      "加入模型的第二个样本层面变量：配对的治疗前后样本用患者 ID，或者批次。不确定就选无。"),
      shiny::uiOutput(ns("cov_ui"))),
    label_with_help("Cell-type column",
                    "The cell populations to analyse one by one. Defaults to the annotation (celltype) when there is one, else the active clustering.",
                    "细胞类型列",
                    "逐个分析的细胞群。有注释（celltype）时默认用注释，否则用当前聚类。"),
    shiny::uiOutput(ns("group_ui")),
    shiny::tags$details(
      class = "omicone-sheet",
      shiny::tags$summary(i18n("No condition column? Add a sample sheet",
                               "没有条件列？添加样本表")),
      shiny::div(class = "omicone-help-text",
                 i18n("A CSV / TSV with one row per sample: the first column holds the sample ids (as in the sample column), the other columns become new metadata columns (e.g. <code>response</code>).",
                      "CSV / TSV，每个样本一行：第一列是样本 ID（与样本列一致），其余各列会成为新的元数据列（例如 <code>response</code>）。")),
      shiny::fileInput(ns("sheet"), NULL, accept = c(".csv", ".tsv", ".txt"),
                       buttonLabel = i18n("Browse…", "选择文件…")))
  )
}

#' @rdname mod_sample_design
#' @param input,output,session The calling module's.
#' @param rv Shared hub.
#' @return A reactive returning list(sample, condition, group, ref, alt, covariate).
#' @keywords internal
sample_design_server <- function(input, output, session, rv, levels = TRUE, covariate = FALSE) {
  cols <- shiny::reactive({
    md <- obj_meta(rv$obj)
    if (is.null(md)) character(0) else sc_group_cols(md)
  })

  output$sample_ui <- shiny::renderUI({
    cc <- cols()
    if (!length(cc)) {
      return(shiny::div(class = "omicone-placeholder",
                        i18n("Load a dataset with several samples first.", "请先载入包含多个样本的数据。")))
    }
    default <- guess_batch_col(obj_meta(rv$obj)) %||% cc[1]
    if (!default %in% cc) default <- cc[1]
    shiny::selectInput(session$ns("sample"), NULL, choices = cc,
                       selected = keep_selected(shiny::isolate(input$sample), cc, default))
  })

  output$cond_ui <- shiny::renderUI({
    cc <- setdiff(cols(), input$sample)
    if (!length(cc)) {
      return(shiny::div(class = "omicone-placeholder",
                        i18n("No other categorical column: add a sample sheet below.",
                             "没有其他分类列：请在下方添加样本表。")))
    }
    hit <- cc[grepl("condition|group|response|treat|status|arm|timepoint|disease", cc, ignore.case = TRUE)]
    shiny::selectInput(session$ns("condition"), NULL, choices = cc,
                       selected = keep_selected(shiny::isolate(input$condition), cc,
                                                if (length(hit)) hit[1] else cc[1]))
  })

  cond_levels <- shiny::reactive({
    md <- obj_meta(rv$obj)
    shiny::req(md, input$condition %in% names(md))
    lv <- md[[input$condition]]
    if (is.factor(lv)) levels(droplevels(lv)) else sort(unique(as.character(lv[!is.na(lv)])))
  })

  if (levels) {
    output$levels_ui <- shiny::renderUI({
      lv <- cond_levels()
      if (length(lv) < 2) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("This column has a single value.", "该列只有一个取值。")))
      }
      ref <- keep_selected(shiny::isolate(input$ref), lv, lv[1])
      alt <- keep_selected(shiny::isolate(input$alt), setdiff(lv, ref), setdiff(lv, ref)[1])
      shiny::div(class = "omicone-row2",
        shiny::div(label_with_help("Reference", "The baseline group (logFC < 0 = lower in the comparison group).",
                                   "参照组", "基线组（logFC < 0 = 比较组中更低）。"),
                   shiny::selectInput(session$ns("ref"), NULL, choices = lv, selected = ref)),
        shiny::div(label_with_help("Compared group", "Fold changes are this group versus the reference.",
                                   "比较组", "倍数变化 = 该组相对参照组。"),
                   shiny::selectInput(session$ns("alt"), NULL, choices = lv, selected = alt)))
    })
  }

  if (covariate) {
    output$cov_ui <- shiny::renderUI({
      cc <- setdiff(cols(), c(input$sample, input$condition))
      ch <- c(stats::setNames("", "none"), stats::setNames(cc, cc))
      shiny::selectInput(session$ns("covariate"), NULL, choices = ch,
                         selected = keep_selected(shiny::isolate(input$covariate), ch, ""))
    })
  }

  output$group_ui <- shiny::renderUI({
    cc <- cols()
    shiny::req(length(cc))
    active <- obj_misc(rv$obj, "omicone_cluster_col")
    default <- intersect(c("celltype", active, "seurat_clusters"), cc)
    default <- if (length(default)) default[1] else cc[1]
    shiny::selectInput(session$ns("group"), NULL, choices = cc,
                       selected = keep_selected(shiny::isolate(input$group), cc, default))
  })

  shiny::observeEvent(input$sheet, {
    shiny::req(rv$obj, input$sample)
    f <- input$sheet
    sheet <- tryCatch(utils::read.delim(f$datapath, sep = guess_table_sep(f$name), check.names = FALSE,
                                        stringsAsFactors = FALSE, colClasses = "character"),
                      error = function(e) NULL)
    if (is.null(sheet) || ncol(sheet) < 2) {
      shiny::showNotification(i18n("Could not read the sheet: it needs a sample column and at least one more column.",
                                   "无法读取样本表：需要样本列和至少一个其他列。"), type = "error")
      return(NULL)
    }
    md <- obj_meta(rv$obj)
    j <- tryCatch(apply_sample_sheet(md, sheet, input$sample), error = function(e) conditionMessage(e))
    if (is.character(j)) {
      shiny::showNotification(j, type = "error")
      return(NULL)
    }
    if (!length(j$matched)) {
      shiny::showNotification(i18n(sprintf("No sample id in the sheet matches the column '%s'.", input$sample),
                                   sprintf("样本表中没有与列「%s」匹配的样本 ID。", input$sample)),
                              type = "error")
      return(NULL)
    }
    obj <- rv$obj
    for (nm in names(j$cols)) obj[[nm]] <- j$cols[[nm]]
    rv$obj <- obj
    shiny::showNotification(
      i18n(sprintf("Added %s for %d sample(s)%s.", paste(names(j$cols), collapse = ", "), length(j$matched),
                   if (length(j$missing)) sprintf("; %d sample(s) not in the sheet are NA (%s)", length(j$missing),
                                                  paste(utils::head(j$missing, 5), collapse = ", ")) else ""),
           sprintf("已为 %d 个样本添加 %s%s。", length(j$matched), paste(names(j$cols), collapse = "、"),
                   if (length(j$missing)) sprintf("；%d 个不在样本表中的样本为 NA（%s）", length(j$missing),
                                                  paste(utils::head(j$missing, 5), collapse = "、")) else "")),
      type = if (length(j$missing)) "warning" else "message", duration = 8)
  })

  shiny::reactive(list(sample = input$sample, condition = input$condition, group = input$group,
                       ref = input$ref, alt = input$alt,
                       covariate = if (isTRUE(nzchar(input$covariate %||% ""))) input$covariate else NULL))
}
