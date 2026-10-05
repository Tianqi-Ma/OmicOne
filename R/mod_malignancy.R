#' Module: Malignant cells / CNV
#'
#' Infer chromosome-scale copy-number changes with copykat (one run per
#' sample) to separate aneuploid (likely malignant) from diploid cells, and
#' score a user-defined stemness gene set. Optional, cancer data only.
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
      en = "Infer chromosome-scale copy-number changes to tell aneuploid (likely
            malignant) cells from diploid ones, and score a stemness gene set you
            supply. <b>Optional, cancer data only.</b>",
      zh = "推断染色体尺度的拷贝数变化，以区分非整倍体（可能为恶性）与二倍体细胞，并对你提供的干性基因集打分。<b>可选步骤，仅适用于癌症数据。</b>"),
    why  = list(
      en = "Tumours carry broad gains and losses that normal cells lack; detecting
            them separates cancer cells from the microenvironment without relying
            on marker genes.",
      zh = "肿瘤携带正常细胞所没有的大范围扩增与缺失；检测它们可以不依赖标志基因，把癌细胞与微环境区分开。"),
    how  = list(
      en = "copykat runs once per sample (pick the sample column). Choosing known
            normal groups (immune / stromal) as reference anchors the diploid
            baseline; without one copykat infers it. Mouse data use the mm10
            genome. The stemness score is a UCell score of your genes, not
            mRNAsi or CytoTRACE.",
      zh = "copykat 按样本逐一运行（请选择样本列）。选择已知的正常分组（免疫 / 基质）作为参考，可锚定二倍体基线；不选时由 copykat 自行推断。小鼠数据使用 mm10 基因组。干性分数是你所给基因的 UCell 分数，并非 mRNAsi 或 CytoTRACE。"),
    read = list(
      en = "The map colours cells by the copykat call: malignant (aneuploid,
            including low-confidence calls), normal (diploid) or no call. Check
            that the reference groups come out normal; if many are called
            malignant the reference or the calls are suspect.",
      zh = "图中细胞按 copykat 判定着色：恶性（非整倍体，含低置信判定）、正常（二倍体）或无判定。请确认参考分组大多被判为正常；若其中很多被判为恶性，说明参考选择或判定结果可疑。"),
    example = list(
      en = "Epithelial cells called <b>malignant</b> with the T and B cells of the
            same sample called <b>normal</b> is the expected pattern.",
      zh = "上皮细胞被判为<b>恶性</b>，而同一样本的 T、B 细胞被判为<b>正常</b>，是预期的模式。")
  )
  controls <- shiny::tagList(
    label_with_help("CNV method",
                    "copykat runs in-app. inferCNV and Numbat (*) need their own setup (gene positions, allele counts) and are not run here.",
                    label_zh = "CNV 方法",
                    tip_zh = "copykat 在应用内运行。inferCNV 与 Numbat（*）需要各自的额外准备（基因位置、等位基因计数），不在此运行。"),
    shiny::selectInput(ns("method"), NULL,
                       choices = c("copykat"    = "copykat",
                                   "inferCNV *" = "infercnv",
                                   "Numbat *"   = "numbat"),
                       selected = "copykat"),
    label_with_help("Sample column",
                    "copykat runs once per sample so baselines never mix samples.",
                    label_zh = "样本列",
                    tip_zh = "copykat 按样本逐一运行，使基线不会跨样本混合。"),
    shiny::selectInput(ns("sample_col"), NULL, choices = c("(one run, all cells)" = "")),
    label_with_help("Reference (normal) cell-type column",
                    "Metadata column with cell-type / cluster labels used to pick the normal reference cells.",
                    label_zh = "参考（正常）细胞类型列",
                    tip_zh = "包含细胞类型/聚类标签的元数据列，用于选取正常参考细胞。"),
    shiny::selectInput(ns("ref_col"), NULL, choices = NULL),
    label_with_help("Normal (reference) groups",
                    "Groups of that column that are known normal cells. Leave empty to let copykat infer the baseline.",
                    label_zh = "正常（参考）分组",
                    tip_zh = "该列中已知为正常细胞的分组。留空则由 copykat 自行推断基线。"),
    shiny::selectInput(ns("ref_groups"), NULL, choices = NULL, multiple = TRUE),
    run_button(ns("run_cnv"), "Run CNV", "运行 CNV"),
    shiny::tags$hr(),
    label_with_help("Stemness gene set (user-defined)",
                    "Gene symbols scored with UCell as a user-defined stemness gene-set score (at least 3 found).",
                    label_zh = "干性基因集（自定义）",
                    tip_zh = "以 UCell 计算的自定义干性基因集分数（至少需在数据中找到 3 个基因）。"),
    shiny::textAreaInput(ns("stem_genes"), NULL, rows = 4,
                         placeholder = "e.g. SOX2, POU5F1, NANOG, LIN28A, MYC"),
    run_button(ns("run_stem"), "Score stemness", "计算干性评分"),
    shiny::tags$hr(),
    label_with_help("Show", "Which result the figure shows.",
                    label_zh = "显示", tip_zh = "图中显示哪一项结果。"),
    shiny::radioButtons(ns("show"), NULL,
                        c("CNV calls" = "calls", "Stemness score" = "stem"),
                        selected = "calls", inline = TRUE)
  )
  step_container(title = list(en = "Malignant cells / CNV", zh = "恶性细胞 / CNV"),
                 subtitle = list(en = "Separate malignant cells by inferred copy-number changes.",
                                 zh = "按推断的拷贝数变化区分恶性细胞。"),
                 explainer = explainer, controls = controls,
                 summary = shiny::uiOutput(ns("summary")),
                 preview = shiny::tagList(
                   shiny::uiOutput(ns("insight")),
                   preview_plot_ui(ns("preview"), download = TRUE,
                     guide = list(en = "Per-cell copykat calls or the stemness score will be drawn here.",
                                  zh = "运行后，这里将绘制每个细胞的 copykat 判定或干性分数。"),
                     caption = list(en = "Cells on the embedding, coloured by copykat call (malignant / normal / no call) or by the stemness score.",
                                    zh = "嵌入图上的细胞，按 copykat 判定（恶性 / 正常 / 无判定）或干性分数着色。"))))
}

#' @rdname mod_malignancy
#' @keywords internal
mod_malignancy_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    res <- step_results(rv, "sc", cnv_done = FALSE, method = NULL, n_malignant = NA,
                        n_low = NA, n_normal = NA, n_nocall = NA, n_samples = NA,
                        skipped = character(0), ref_n = 0L, ref_flag = NA,
                        ref_groups = character(0), genome = NULL,
                        stem_done = FALSE, stem_cov = NULL)

    shiny::observe({
      obj <- rv$obj
      md <- obj_meta(obj)
      cols <- categorical_cols(md)
      def <- default_group_col(cols, obj_misc(obj, "omicone_cluster_col"))
      shiny::updateSelectInput(session, "ref_col", choices = cols,
                               selected = keep_selected(shiny::isolate(input$ref_col), cols, def))
      smp <- c("(one run, all cells)" = "", stats::setNames(cols, cols))
      guess <- guess_batch_col(md[, cols, drop = FALSE]) %||% ""
      shiny::updateSelectInput(session, "sample_col", choices = smp,
                               selected = keep_selected(shiny::isolate(input$sample_col),
                                                        unname(smp), guess))
    })

    shiny::observe({
      lv <- group_levels(obj_meta(rv$obj), input$ref_col)
      cur <- shiny::isolate(input$ref_groups)
      shiny::updateSelectInput(session, "ref_groups", choices = lv,
                               selected = intersect(cur, lv))
    })

    # ---- Run CNV / malignant-cell calling -----------------------------------
    shiny::observeEvent(input$run_cnv, {
      shiny::req(rv$obj)
      method <- input$method
      ref_col <- input$ref_col
      ref_groups <- input$ref_groups %||% character(0)
      sample_col <- input$sample_col %||% ""
      if (!identical(method, "copykat")) {
        shiny::showNotification(
          i18n(sprintf("%s needs gene positions or allele counts prepared outside the app; run it with scop::RunCNV() in R.", method),
               sprintf("%s 需要在应用外准备基因位置或等位基因计数；请在 R 中用 scop::RunCNV() 运行。", method)),
          type = "warning", duration = 12)
        return(NULL)
      }
      if (!require_pkgs(c("copykat", "Seurat"), "CNV / malignant cells")) return(NULL)
      md <- obj_meta(rv$obj)
      ref_cells <- if (length(ref_groups) && !is.null(ref_col) && ref_col %in% names(md)) {
        rownames(md)[as.character(md[[ref_col]]) %in% ref_groups]
      } else {
        character(0)
      }
      species <- guess_species(rv$obj)
      smp <- if (nzchar(sample_col)) sample_col else NULL
      out <- with_progress_notify({
        sc_copykat(rv$obj, ref_cells = ref_cells, sample_col = smp, species = species)
      }, message = "Running copykat per sample (this can take several minutes)...")
      if (is.null(out)) return(NULL)
      obj <- copykat_add_calls(rv$obj, out$pred)
      rv$obj <- obj
      calls <- stats::setNames(obj$malignant, colnames(obj))
      res$cnv_done    <- TRUE
      res$method      <- method
      res$genome      <- copykat_genome(species)
      res$n_malignant <- sum(calls == "malignant", na.rm = TRUE)
      res$n_low       <- sum(calls == "malignant" & obj$malignant_confidence == "low", na.rm = TRUE)
      res$n_normal    <- sum(calls == "normal", na.rm = TRUE)
      res$n_nocall    <- sum(is.na(calls))
      res$n_samples   <- length(unique(out$pred$sample))
      res$skipped     <- out$skipped
      res$ref_n       <- length(ref_cells)
      res$ref_groups  <- ref_groups
      res$ref_flag    <- copykat_ref_flagged(calls, ref_cells)
      mark_done(rv, "malignancy")
      log_step(log_rv, "Malignant cells / CNV",
               params = list(method = "copykat", genome = copykat_genome(species),
                             sample_col = smp, ref_col = if (length(ref_groups)) ref_col,
                             ref_groups = ref_groups),
               code = copykat_log_code(ref_col = ref_col, ref_groups = ref_groups,
                                       sample_col = smp, species = species))
      shiny::updateRadioButtons(session, "show", selected = "calls")
      if (!is.na(res$ref_flag) && res$ref_flag > 0.1) {
        shiny::showNotification(
          i18n(sprintf("%.0f%% of the reference cells were called malignant. Check the reference groups before using the calls.",
                       100 * res$ref_flag),
               sprintf("%.0f%% 的参考细胞被判为恶性。使用这些判定结果之前，请检查参考分组。",
                       100 * res$ref_flag)),
          type = "warning", duration = 15)
      }
      shiny::showNotification(i18n("copykat finished.", "copykat 运行完成。"), type = "message")
    })

    # ---- Stemness scoring ----------------------------------------------------
    shiny::observeEvent(input$run_stem, {
      shiny::req(rv$obj)
      genes <- parse_genes(input$stem_genes)
      cov <- geneset_coverage(list(stemness = genes), rownames(rv$obj))
      if (cov$n_found < 3) {
        shiny::showNotification(
          i18n(sprintf("The stemness set needs at least 3 genes present in the data; %d of %d found.",
                       cov$n_found, cov$n_input),
               sprintf("干性基因集在数据中至少需要 3 个基因；目前找到 %d / %d 个。",
                       cov$n_found, cov$n_input)),
          type = "warning", duration = 10)
        return(NULL)
      }
      if (!require_pkgs(c("Seurat", "UCell"), "Stemness scoring")) return(NULL)
      out <- with_progress_notify({
        sc_stemness(rv$obj, genes = genes)
      }, message = "Scoring the stemness gene set...")
      if (is.null(out)) return(NULL)
      rv$obj <- out$obj
      res$stem_done <- TRUE
      res$stem_cov  <- out$coverage
      mark_done(rv, "malignancy")
      log_step(log_rv, "Stemness score",
               params = list(method = "UCell", genes_found = out$coverage$found[[1]],
                             n_input = out$coverage$n_input),
               code = stemness_log_code(out$coverage$found[[1]]))
      shiny::updateRadioButtons(session, "show", selected = "stem")
      shiny::showNotification(
        i18n(sprintf("Stemness scored on %d gene(s).", out$coverage$n_found),
             sprintf("已基于 %d 个基因计算干性评分。", out$coverage$n_found)),
        type = "message")
    })

    output$summary <- shiny::renderUI({
      if (!isTRUE(res$cnv_done) && !isTRUE(res$stem_done)) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("Optional (cancer data). Run CNV or score stemness.",
                               "可选（癌症数据）。运行 CNV 或计算干性评分。")))
      }
      tiles <- list()
      if (isTRUE(res$cnv_done)) {
        tiles <- c(tiles, list(
          stat_tile(i18n("Malignant", "恶性"), format(res$n_malignant, big.mark = ",")),
          stat_tile(i18n("Normal", "正常"), format(res$n_normal, big.mark = ",")),
          stat_tile(i18n("No call", "无判定"), format(res$n_nocall, big.mark = ","))))
      }
      if (isTRUE(res$stem_done)) {
        tiles <- c(tiles, list(
          stat_tile(i18n("Stemness genes", "干性基因"),
                    sprintf("%d / %d", res$stem_cov$n_found, res$stem_cov$n_input))))
      }
      do.call(shiny::tagList, tiles)
    })

    output$insight <- shiny::renderUI({
      parts_en <- character(0)
      parts_zh <- character(0)
      if (isTRUE(res$cnv_done)) {
        skip_en <- if (length(res$skipped)) {
          sprintf(" Skipped (under 50 cells): %s.", paste(res$skipped, collapse = ", "))
        } else ""
        skip_zh <- if (length(res$skipped)) {
          sprintf("跳过（不足 50 个细胞）：%s。", paste(res$skipped, collapse = "、"))
        } else ""
        ref_en <- if (res$ref_n > 0) {
          sprintf(" Reference: %s cells of %s; %.0f%% of them were called malignant%s.",
                  format(res$ref_n, big.mark = ","), paste(res$ref_groups, collapse = ", "),
                  100 * res$ref_flag,
                  if (!is.na(res$ref_flag) && res$ref_flag > 0.1)
                    " (high: check the reference or the calls)" else "")
        } else " No reference given: copykat inferred the diploid baseline."
        ref_zh <- if (res$ref_n > 0) {
          sprintf("参考：%s 中的 %s 个细胞；其中 %.0f%% 被判为恶性%s。",
                  paste(res$ref_groups, collapse = "、"), format(res$ref_n, big.mark = ","),
                  100 * res$ref_flag,
                  if (!is.na(res$ref_flag) && res$ref_flag > 0.1)
                    "（偏高：请检查参考分组或判定结果）" else "")
        } else "未指定参考：由 copykat 自行推断二倍体基线。"
        parts_en <- c(parts_en, sprintf(
          "copykat (%s) in %d sample(s): <b>%s</b> malignant (%s low-confidence), %s normal, %s without a call.%s%s",
          res$genome, res$n_samples, format(res$n_malignant, big.mark = ","),
          format(res$n_low, big.mark = ","), format(res$n_normal, big.mark = ","),
          format(res$n_nocall, big.mark = ","), skip_en, ref_en))
        parts_zh <- c(parts_zh, sprintf(
          "copykat（%s），%d 个样本：<b>%s</b> 个恶性（其中低置信 %s 个），%s 个正常，%s 个无判定。%s%s",
          res$genome, res$n_samples, format(res$n_malignant, big.mark = ","),
          format(res$n_low, big.mark = ","), format(res$n_normal, big.mark = ","),
          format(res$n_nocall, big.mark = ","), skip_zh, ref_zh))
      }
      if (isTRUE(res$stem_done)) {
        parts_en <- c(parts_en, sprintf(
          "Stemness: UCell score of %d of %d supplied genes; a user-defined gene-set score, not mRNAsi or CytoTRACE.",
          res$stem_cov$n_found, res$stem_cov$n_input))
        parts_zh <- c(parts_zh, sprintf(
          "干性：%d / %d 个所给基因的 UCell 分数；这是自定义基因集分数，并非 mRNAsi 或 CytoTRACE。",
          res$stem_cov$n_found, res$stem_cov$n_input))
      }
      if (!length(parts_en)) return(NULL)
      insight_bar(paste(parts_en, collapse = " "), paste(parts_zh, collapse = ""))
    })

    render_step_plot(output, input, "preview", function() {
      obj <- rv$obj
      shiny::req(obj)
      cols <- obj_meta_cols(obj)
      if (identical(input$show, "stem")) {
        shiny::req("stemness_UCell" %in% cols)
        return(sc_featureplot(obj, features = "stemness_UCell"))
      }
      shiny::req("malignant" %in% cols, isTRUE(res$cnv_done))
      # show_na: cells copykat filtered or skipped stay visible as "no call"
      sc_dimplot(obj, group_by = "malignant", show_na = TRUE)
    }, name = "malignancy")
  })
}
