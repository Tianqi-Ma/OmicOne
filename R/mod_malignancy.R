#' Module: Malignant cells / CNV
#'
#' Estimate large-scale copy-number variation (CNV) to separate malignant
#' (aneuploid) tumour cells from normal (diploid) cells, and score a stemness
#' signature. This step is OPTIONAL and only meaningful for cancer datasets.
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_malignancy
NULL

#' @rdname mod_malignancy
#' @keywords internal
mod_malignancy_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = list(en = "Malignant cells / CNV", zh = "恶性细胞 / CNV"),
    what = list(
      en = "Infer chromosome-scale copy-number changes to tell malignant
            (aneuploid) tumour cells apart from normal (diploid) cells, and
            score a stemness signature. <b>Optional — only for cancer data.</b>",
      zh = "推断染色体尺度的拷贝数变化，以区分恶性（非整倍体）肿瘤细胞与正常（二倍体）细胞，并对干性特征打分。<b>可选步骤 —— 仅适用于癌症数据。</b>"),
    why  = list(
      en = "Tumours carry broad copy-number gains/losses that normal cells lack.
            Detecting them separates cancer cells from the microenvironment and
            reveals malignant clones.",
      zh = "肿瘤携带正常细胞所没有的大范围拷贝数增益/缺失。检测它们可将癌细胞与微环境区分开，并揭示恶性克隆。"),
    how  = list(
      en = "Pick a method and, if you can, a set of known-normal reference cells
            (e.g. immune/stromal cell types). CopyKAT works without a reference;
            inferCNV and Numbat need extra setup. Stemness scores an mRNAsi-style
            gene signature you provide.",
      zh = "选择一种方法；如果可能，指定一组已知的正常参考细胞（例如免疫/基质细胞类型）。CopyKAT 无需参考；inferCNV 与 Numbat 需要额外配置。干性打分基于你提供的 mRNAsi 风格基因特征。"),
    read = list(
      en = "The heatmap shows chromosome-scale gains (warm colours) and losses
            (cool) per cell; the score ranks each cell's malignancy. The
            reference cells anchor what 'normal' looks like — tumour cells are
            the ones that depart from it.",
      zh = "热图展示每个细胞染色体尺度的扩增（暖色）与缺失（冷色）；评分给出每个细胞的恶性程度。参考细胞锚定了「正常」的基线——偏离它的即肿瘤细胞。"),
    example = list(
      en = "Epithelial cells flagged <b>aneuploid</b> are likely the tumour;
               matched <b>diploid</b> T/B cells act as the normal reference.",
      zh = "被标记为<b>非整倍体</b>的上皮细胞很可能是肿瘤；配对的<b>二倍体</b> T/B 细胞则充当正常参考。")
  )
  controls <- shiny::tagList(
    label_with_help("CNV method",
                    "CopyKAT works without a reference. inferCNV and Numbat need extra setup (* = extra install).",
                    label_zh = "CNV 方法",
                    tip_zh = "CopyKAT 无需参考。inferCNV 与 Numbat 需要额外配置（* = 需额外安装）。"),
    shiny::selectInput(ns("method"), NULL,
                       choices = c("CopyKAT"    = "copykat",
                                   "inferCNV *" = "infercnv",
                                   "Numbat *"   = "numbat"),
                       selected = "copykat"),
    label_with_help("Reference (normal) cell-type column",
                    "Metadata column holding cell-type / cluster labels used to pick the normal reference cells.",
                    label_zh = "参考（正常）细胞类型列",
                    tip_zh = "包含细胞类型/聚类标签的元数据列，用于选取正常参考细胞。"),
    shiny::uiOutput(ns("ref_col_ui")),
    label_with_help("Normal (reference) groups",
                    "Which groups of the chosen column are known-normal cells. Leave empty to run without a reference.",
                    label_zh = "正常（参考）分组",
                    tip_zh = "所选列中哪些分组是已知的正常细胞。留空则在无参考的情况下运行。"),
    shiny::uiOutput(ns("ref_groups_ui")),
    run_button(ns("run_cnv"), "Run CNV", "运行 CNV"),
    shiny::tags$hr(),
    label_with_help("Stemness gene signature",
                    "Comma / newline separated gene symbols scored as a stemness (mRNAsi-style) module.",
                    label_zh = "干性基因特征",
                    tip_zh = "以逗号/换行分隔的基因符号，作为干性（mRNAsi 风格）模块进行打分。"),
    shiny::textAreaInput(ns("stem_genes"), NULL, rows = 4,
                         placeholder = "e.g. SOX2, POU5F1, NANOG, LIN28A, MYC"),
    run_button(ns("run_stem"), "Score stemness", "计算干性评分")
  )
  step_container(title = list(en = "Malignant cells / CNV", zh = "恶性细胞 / CNV"),
                 subtitle = list(en = "Separate malignant cells by CNV-like signal.",
                                 zh = "按 CNV 样信号区分恶性细胞。"),
                 explainer = explainer, controls = controls,
                 summary = shiny::uiOutput(ns("summary")),
                 preview = preview_plot_ui(ns("preview"),
                   guide = list(en = "The CNV-style heatmap and malignancy scores will be drawn here.",
                                zh = "运行后，这里将绘制 CNV 热图与恶性评分。"),
                   caption = list(en = "Chromosome-scale gain/loss heatmap per cell.",
                                  zh = "每细胞染色体尺度扩增/缺失热图。")))
}

#' @rdname mod_malignancy
#' @keywords internal
mod_malignancy_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    res <- shiny::reactiveValues(cnv_done = FALSE, method = NULL, n_malignant = NA,
                                 n_normal = NA, stem_done = FALSE, n_stem_genes = NA)

    # Parse a free-text gene list (commas, whitespace or newlines).
    parse_genes <- function(txt) {
      if (is.null(txt) || !nzchar(trimws(txt))) return(character(0))
      g <- trimws(unlist(strsplit(txt, "[,\\s]+", perl = TRUE)))
      unique(g[nzchar(g)])
    }

    # Reference-column selector: any metadata column (cell-type / clusters).
    output$ref_col_ui <- shiny::renderUI({
      cols <- if (is.null(rv$obj)) character(0) else obj_meta_cols(rv$obj)
      if (length(cols) == 0) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("No metadata yet — cluster or annotate first.",
                               "还没有元数据 —— 请先聚类或注释。")))
      }
      sel <- if ("celltype" %in% cols) "celltype"
             else if ("seurat_clusters" %in% cols) "seurat_clusters"
             else cols[1]
      shiny::selectInput(ns("ref_col"), NULL, choices = cols, selected = sel)
    })

    # Normal-group selector: levels of the chosen reference column.
    output$ref_groups_ui <- shiny::renderUI({
      shiny::req(rv$obj, input$ref_col)
      md <- obj_meta(rv$obj)
      if (!input$ref_col %in% colnames(md)) return(NULL)
      levs <- sort(unique(as.character(md[[input$ref_col]])))
      shiny::selectInput(ns("ref_groups"), NULL, choices = levs,
                         selected = NULL, multiple = TRUE)
    })

    # ---- Run CNV / malignant-cell calling -----------------------------------
    shiny::observeEvent(input$run_cnv, {
      shiny::req(rv$obj)
      if (!require_pkgs("copykat", "CNV / malignant cells")) return(NULL)
      method <- input$method

      # Cell names of the chosen normal reference groups (if any).
      ref_cells <- NULL
      if (!is.null(input$ref_col) && length(input$ref_groups)) {
        md <- obj_meta(rv$obj)
        if (input$ref_col %in% colnames(md)) {
          ref_cells <- rownames(md)[as.character(md[[input$ref_col]]) %in%
                                      input$ref_groups]
        }
      }

      out <- with_progress_notify({
        sc_cnv(rv$obj, method = method, ref_cells = ref_cells)
      }, message = sprintf("Running %s (CNV)...", method))
      if (is.null(out)) return(NULL)
      rv$cnv <- out

      # If the tool returns per-cell malignant calls, write them to the object.
      # CopyKAT returns a list with a `prediction` data.frame (cell.names,
      # copykat.pred = "aneuploid"/"diploid").
      calls <- NULL
      pred <- tryCatch(out$prediction, error = function(e) NULL)
      if (!is.null(pred) && all(c("cell.names", "copykat.pred") %in% colnames(pred))) {
        malignant <- ifelse(pred$copykat.pred == "aneuploid", "malignant", "normal")
        calls <- stats::setNames(malignant, pred$cell.names)
      }
      if (!is.null(calls)) {
        md <- obj_meta(rv$obj)
        vec <- rep(NA_character_, nrow(md))
        names(vec) <- rownames(md)
        common <- intersect(names(calls), names(vec))
        vec[common] <- calls[common]
        rv$obj$malignant <- vec
        res$n_malignant <- sum(vec == "malignant", na.rm = TRUE)
        res$n_normal    <- sum(vec == "normal", na.rm = TRUE)
      } else {
        res$n_malignant <- NA
        res$n_normal    <- NA
      }

      res$cnv_done <- TRUE
      res$method   <- method
      mark_done(rv, "malignancy")
      log_step(log_rv, "Malignant cells / CNV",
               params = list(method = method,
                             ref_col = input$ref_col %||% "",
                             ref_groups = input$ref_groups %||% character(0)),
               code = sprintf(
                 "cnv <- sc_cnv(obj, method = '%s', ref_cells = ref_cells)  # writes obj$malignant",
                 method))
      shiny::showNotification(
        i18n(sprintf("CNV (%s) done.", method), sprintf("CNV（%s）完成。", method)),
        type = "message")
    })

    # ---- Stemness scoring ----------------------------------------------------
    shiny::observeEvent(input$run_stem, {
      shiny::req(rv$obj)
      genes <- parse_genes(input$stem_genes)
      if (length(genes) == 0) {
        shiny::showNotification(
          i18n("Enter at least one gene for the stemness signature.",
               "请为干性特征输入至少一个基因。"),
          type = "warning")
        return(NULL)
      }
      if (!require_pkgs("scop", "Stemness scoring")) return(NULL)
      obj <- with_progress_notify({
        sc_stemness(rv$obj, features = genes)
      }, message = "Scoring stemness signature...")
      if (is.null(obj)) return(NULL)
      rv$obj <- obj
      res$stem_done <- TRUE
      res$n_stem_genes <- length(genes)
      mark_done(rv, "malignancy")
      log_step(log_rv, "Stemness score",
               params = list(n_genes = length(genes),
                             genes = paste(genes, collapse = ", ")),
               code = "obj <- sc_stemness(obj, features = stemness_genes)")
      shiny::showNotification(
        i18n(sprintf("Stemness scored on %d gene(s).", length(genes)),
             sprintf("已基于 %d 个基因计算干性评分。", length(genes))),
        type = "message")
    })

    output$summary <- shiny::renderUI({
      if (!isTRUE(res$cnv_done) && !isTRUE(res$stem_done)) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("Optional (cancer data). Run CNV or score stemness.",
                               "可选（癌症数据）。运行 CNV 或计算干性评分。")))
      }
      tiles <- list(stat_tile(i18n("Method", "方法"),
                              if (is.null(res$method)) i18n("-", "-") else res$method))
      if (!is.na(res$n_malignant)) {
        tiles <- c(tiles, list(
          stat_tile(i18n("Malignant", "恶性"), format(res$n_malignant, big.mark = ",")),
          stat_tile(i18n("Normal", "正常"), format(res$n_normal, big.mark = ","))))
      }
      if (isTRUE(res$stem_done)) {
        tiles <- c(tiles, list(
          stat_tile(i18n("Stemness genes", "干性基因"), res$n_stem_genes)))
      }
      n <- length(tiles)
      w <- rep(floor(12 / max(n, 1)), n)
      do.call(bslib::layout_columns, c(list(col_widths = w), tiles))
    })

    output$preview <- render_preview_plot(function() {
      obj <- rv$obj
      shiny::req(obj)
      # Prefer the malignant-call dim plot when we have per-cell calls.
      if ("malignant" %in% obj_meta_cols(obj)) {
        p <- tryCatch(sc_dimplot(obj, group_by = "malignant"),
                      error = function(e) NULL)
        if (!is.null(p)) return(p)
      }
      # Otherwise, if the CNV tool produced its own plot, try to show it.
      if (!is.null(rv$cnv)) {
        p <- tryCatch({
          cand <- rv$cnv$plot %||% NULL
          if (methods::is(cand, "ggplot") || methods::is(cand, "Heatmap") ||
              methods::is(cand, "HeatmapList")) cand else NULL
        }, error = function(e) NULL)
        if (!is.null(p)) return(p)
      }
      # Fallback message.
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0, y = 0, size = 5,
                          label = paste0(
                            "Run CNV to call malignant cells.\n",
                            "The CNV heatmap opens in the tool's own window.")) +
        ggplot2::theme_void()
    })
  })
}
