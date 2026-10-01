#' Module: Clinical data & survival
#'
#' Load a per-sample clinical table (follow-up time + outcome), then ask whether
#' anything in the analysis separates patients. Two ways to stratify: by a column
#' of the clinical table itself (stage, treatment, a score you already have), or
#' by the composition of the working object -- the fraction of a chosen cluster
#' or cell type in each sample, split high vs low.
#'
#' The normalised table is written to `rv$clinical`, the hub slot every omics
#' shares, so WES / bulk / integration reuse the same cohort without reloading
#' it. The maths lives in `fct_survival.R`.
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_clinical
NULL

#' @rdname mod_clinical
#' @keywords internal
mod_clinical_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = list(en = "Clinical data & survival", zh = "临床数据与生存分析"),
    what = list(
      en = "Attach patient follow-up (how long each patient was observed, and
            whether the event happened) and test whether a grouping separates
            their survival.",
      zh = "接入患者随访信息（每位患者被观察了多久、终点事件是否发生），并检验某种分组是否能区分生存差异。"),
    why  = list(
      en = "A cluster is only interesting if it means something. Linking the
            composition of your samples to outcome turns a descriptive atlas
            into a clinical claim -- and the same cohort is reused by the WES,
            bulk and integration pipelines.",
      zh = "一个细胞簇只有具备实际意义才值得关注。把样本组成与预后联系起来，能让描述性图谱变成有临床意义的结论——同一份队列还会被 WES、Bulk 和整合流程复用。"),
    how  = list(
      en = "Upload a CSV/TSV with <b>one row per patient</b>, then point the app
            at the id, time and outcome columns. Stratify either by a clinical
            column, or by <b>cell composition</b> (the fraction of one cell type
            per sample, split high vs low). Use <b>optimal</b> only as an
            exploratory cutpoint -- it is selected to maximise separation, so
            its p-value is optimistic.",
      zh = "上传<b>每位患者一行</b>的 CSV/TSV 文件，然后指定 ID、时间和终点列。可按临床列分组，也可按<b>细胞组成</b>分组（每个样本中某种细胞类型的占比，分高低两组）。<b>最优切点</b>仅供探索——它是为了让分离最大化而挑选的，因此 p 值偏乐观。"),
    read = list(
      en = "KM curves split patients by a high versus low score. Separation with
            heavily overlapping confidence bands is suggestive, not proof —
            check the group sizes before believing it.",
      zh = "KM 曲线按评分高低把患者分组。置信带大量重叠的分离只能算提示而非证明——采信前请核对各组样本量。"),
    example = list(
      en = "A table with <code>sample, os_months, os_status, stage</code>, then
               \"fraction of Exhausted CD8 T cells, high vs low\" as the grouping.",
      zh = "一张包含 <code>sample, os_months, os_status, stage</code> 的表格，再以「耗竭 CD8 T 细胞占比的高/低」作为分组。")
  )

  controls <- shiny::tagList(
    label_with_help("Clinical table",
                    "CSV or TSV, one row per patient. Sample ids must match the sample column in your object if you want to stratify by composition.",
                    label_zh = "临床表格",
                    tip_zh = "CSV 或 TSV，每位患者一行。若要按细胞组成分组，样本 ID 需与对象中的样本列一致。"),
    shiny::fileInput(ns("file"), i18n("Choose file", "选择文件"),
                     accept = c(".csv", ".tsv", ".txt")),
    shiny::selectInput(ns("sep"), i18n("Separator", "分隔符"),
                       choices = c("Comma" = ",", "Tab" = "\t"), selected = ","),
    shiny::uiOutput(ns("mapping")),
    shiny::hr(),
    label_with_help("Stratify by",
                    "A column of the clinical table, or the per-sample fraction of a cell type / cluster taken from the working object.",
                    label_zh = "分组依据",
                    tip_zh = "可选临床表中的某一列，或取自当前对象的每样本细胞类型/簇占比。"),
    shiny::radioButtons(ns("mode"), NULL,
                        c("Clinical column" = "clinical",
                          "Cell composition" = "composition"),
                        selected = "clinical"),
    shiny::uiOutput(ns("grouping")),
    shiny::hr(),
    shiny::uiOutput(ns("covars")),
    run_button(ns("run"), "Run survival analysis", "运行生存分析")
  )

  step_container(
    title     = list(en = "Clinical data & survival", zh = "临床数据与生存分析"),
    subtitle  = list(en = "Survival analysis on cell-level scores.",
                     zh = "基于细胞水平评分的生存分析。"),
    explainer = explainer,
    controls  = controls,
    summary   = shiny::uiOutput(ns("summary")),
    preview   = bslib::navset_card_tab(
      bslib::nav_panel(i18n("Kaplan-Meier", "生存曲线"), preview_plot_ui(ns("km"),
        guide = list(en = "Kaplan-Meier curves will be drawn here.",
                     zh = "运行后，这里将绘制 Kaplan-Meier 生存曲线。"),
        caption = list(en = "KM curves for the high- vs low-score groups.",
                       zh = "高分组与低分组的 KM 曲线。"))),
      bslib::nav_panel(i18n("Cox (univariable)", "Cox 单因素"), shiny::uiOutput(ns("cox_slot"))),
      bslib::nav_panel(i18n("Cohort", "队列表"), shiny::uiOutput(ns("cohort_slot")))
    )
  )
}

#' @rdname mod_clinical
#' @keywords internal
mod_clinical_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    ns  <- session$ns
    raw <- shiny::reactiveVal(NULL)   # clinical table as uploaded
    res <- shiny::reactiveValues(df = NULL, fit = NULL, lr = NULL, cox = NULL,
                                 label = NULL, matched = NA_integer_)

    # ---- 1. read the uploaded table ----------------------------------------
    shiny::observeEvent(input$file, {
      shiny::req(input$file)
      df <- tryCatch(read_clinical_table(input$file$datapath, sep = input$sep),
                     error = function(e) {
                       shiny::showNotification(paste("Could not read the table:",
                                                     conditionMessage(e)),
                                               type = "error", duration = 12)
                       NULL
                     })
      shiny::req(df)
      raw(df)
      shiny::showNotification(
        sprintf("Clinical table loaded: %d rows, %d columns.", nrow(df), ncol(df)),
        type = "message")
    })

    # Column mappers, populated from the uploaded table. Sensible guesses first
    # so the common case is a single click.
    output$mapping <- shiny::renderUI({
      df <- raw()
      if (is.null(df)) {
        return(shiny::div(class = "omicone-status-empty",
                          i18n("Upload a clinical table to continue.",
                               "请先上传临床表格。")))
      }
      cols <- names(df)
      pick <- function(cands, default = cols[1]) {
        hit <- cols[tolower(cols) %in% cands]
        if (length(hit)) hit[1] else default
      }
      shiny::tagList(
        label_with_help("Sample / patient id", "Column that identifies each patient.",
                        label_zh = "样本 / 患者 ID", tip_zh = "标识每位患者的列。"),
        shiny::selectInput(ns("id_col"), NULL, cols,
                           selected = pick(c("sample", "sample_id", "patient",
                                             "patient_id", "id", "case_id",
                                             "bcr_patient_barcode"))),
        label_with_help("Follow-up time", "Time from baseline to event or last contact.",
                        label_zh = "随访时间", tip_zh = "从基线到终点事件或最后一次随访的时间。"),
        shiny::selectInput(ns("time_col"), NULL, cols,
                           selected = pick(c("os_months", "os_time", "os", "time",
                                             "overall_survival", "days_to_death",
                                             "futime"))),
        shiny::selectInput(ns("time_unit"), i18n("Time unit", "时间单位"),
                           c("Months" = "months", "Days" = "days", "Years" = "years"),
                           selected = "months"),
        label_with_help("Outcome / event", "1 / Dead / TRUE = the event happened; 0 / Alive = censored.",
                        label_zh = "终点事件", tip_zh = "1 / Dead / TRUE 表示事件发生；0 / Alive 表示删失。"),
        shiny::selectInput(ns("event_col"), NULL, cols,
                           selected = pick(c("os_status", "status", "event",
                                             "vital_status", "death", "fustat"))),
        shiny::uiOutput(ns("event_level"))
      )
    })

    # Only ask which value means "event" when it cannot be inferred.
    output$event_level <- shiny::renderUI({
      df <- raw()
      shiny::req(df, input$event_col)
      x <- df[[input$event_col]]
      if (!all(is.na(encode_event(x)))) return(NULL)
      lv <- sort(unique(as.character(stats::na.omit(x))))
      shiny::tagList(
        label_with_help("Which value means the event happened?",
                        "The outcome coding was not recognised, so pick it here.",
                        label_zh = "哪个取值表示事件发生？",
                        tip_zh = "未能识别终点事件的编码方式，请在此选择。"),
        shiny::selectInput(ns("event_positive"), NULL, lv)
      )
    })

    # ---- 2. grouping controls ----------------------------------------------
    output$grouping <- shiny::renderUI({
      df <- raw()
      shiny::req(df)
      if (identical(input$mode, "clinical")) {
        cols <- setdiff(names(df), c(input$id_col, input$time_col, input$event_col))
        if (!length(cols)) {
          return(shiny::div(class = "omicone-status-empty",
                            i18n("No other column to group by.", "没有可用于分组的其他列。")))
        }
        return(shiny::tagList(
          shiny::selectInput(ns("clin_col"), i18n("Clinical column", "临床列"), cols),
          shiny::uiOutput(ns("cut_ui"))
        ))
      }
      # composition mode: needs the working object's metadata
      md <- obj_meta(rv$obj)
      if (!ncol(md)) {
        return(shiny::div(class = "omicone-status-empty",
                          i18n("Load and cluster your data first to group by composition.",
                               "请先加载并聚类数据，才能按细胞组成分组。")))
      }
      cats <- names(md)[vapply(md, function(c) is.factor(c) || is.character(c) ||
                                 length(unique(c)) <= 50, logical(1))]
      shiny::tagList(
        label_with_help("Sample column (in the object)",
                        "Which metadata column says what sample each cell came from.",
                        label_zh = "样本列（对象中）",
                        tip_zh = "对象元数据中标明每个细胞来自哪个样本的列。"),
        shiny::selectInput(ns("samp_col"), NULL, names(md),
                           selected = guess_batch_col(md) %||% names(md)[1]),
        label_with_help("Cell grouping", "Clusters or annotated cell types.",
                        label_zh = "细胞分组", tip_zh = "聚类结果或已注释的细胞类型。"),
        shiny::selectInput(ns("cell_col"), NULL, cats,
                           selected = if ("seurat_clusters" %in% cats) "seurat_clusters"
                                      else cats[1]),
        shiny::uiOutput(ns("cell_level")),
        shiny::uiOutput(ns("cut_ui"))
      )
    })

    output$cell_level <- shiny::renderUI({
      md <- obj_meta(rv$obj)
      shiny::req(input$cell_col, input$cell_col %in% names(md))
      lv <- sort(unique(as.character(stats::na.omit(md[[input$cell_col]]))))
      shiny::selectInput(ns("cell_level"), i18n("Cell type / cluster", "细胞类型 / 簇"), lv)
    })

    # The cutpoint control only makes sense for a numeric stratifier.
    output$cut_ui <- shiny::renderUI({
      numeric_mode <- identical(input$mode, "composition") || {
        df <- raw()
        !is.null(df) && !is.null(input$clin_col) && input$clin_col %in% names(df) &&
          is.numeric(df[[input$clin_col]])
      }
      if (!isTRUE(numeric_mode)) return(NULL)
      shiny::tagList(
        label_with_help("Cutpoint",
                        "Median = balanced halves. Tertiles = top vs bottom third (middle dropped). Optimal = the split with the strongest separation -- exploratory, its p-value is optimistic.",
                        label_zh = "切点",
                        tip_zh = "中位数 = 两组样本量均衡。三分位 = 上三分之一对下三分之一（中间组舍弃）。最优 = 分离度最强的切点——仅供探索，p 值偏乐观。"),
        shiny::selectInput(ns("cut"), NULL,
                           c("Median" = "median", "Tertiles" = "tertile",
                             "Optimal" = "optimal"),
                           selected = "median")
      )
    })

    output$covars <- shiny::renderUI({
      df <- raw()
      shiny::req(df)
      cols <- setdiff(names(df), c(input$id_col, input$time_col, input$event_col))
      shiny::tagList(
        label_with_help("Cox covariates (optional)",
                        "Each is fitted in its own univariable Cox model -- a screen, not an adjusted model.",
                        label_zh = "Cox 协变量（可选）",
                        tip_zh = "每个变量各自拟合一个单因素 Cox 模型——这是筛查，不是校正后的多因素模型。"),
        shiny::selectizeInput(ns("cox_vars"), NULL, choices = cols, multiple = TRUE,
                              options = list(placeholder = "stage, age, ..."))
      )
    })

    # ---- 3. run -------------------------------------------------------------
    shiny::observeEvent(input$run, {
      df <- raw()
      shiny::req(df, input$id_col, input$time_col, input$event_col)

      clin <- with_progress_notify(
        normalise_clinical(df, input$id_col, input$time_col, input$event_col,
                           event_positive = input$event_positive,
                           time_unit = input$time_unit %||% "months"),
        message = "Reading clinical table...")
      if (is.null(clin)) return(NULL)
      if (!nrow(clin)) {
        shiny::showNotification(
          "No usable rows: check the time and outcome columns.",
          type = "error", duration = 12)
        return(NULL)
      }
      if (isTRUE(attr(clin, "dropped") > 0)) {
        shiny::showNotification(
          sprintf("%d row(s) dropped: missing or invalid time/outcome.",
                  attr(clin, "dropped")),
          type = "warning", duration = 10)
      }

      # Attach the grouping.
      label <- NULL
      matched <- NA_integer_
      if (identical(input$mode, "clinical")) {
        shiny::req(input$clin_col)
        v <- clin[[input$clin_col]]
        if (is.numeric(v)) {
          clin$.group <- split_numeric(v, input$cut %||% "median",
                                       clin$.time, clin$.event)
          label <- sprintf("%s (%s split)", input$clin_col, input$cut %||% "median")
        } else {
          clin$.group <- factor(as.character(v))
          label <- input$clin_col
        }
      } else {
        shiny::req(rv$obj, input$samp_col, input$cell_col, input$cell_level)
        md <- obj_meta(rv$obj)
        comp <- tryCatch(composition_by_sample(md, input$samp_col, input$cell_col),
                         error = function(e) {
                           shiny::showNotification(conditionMessage(e), type = "error",
                                                   duration = 12); NULL
                         })
        if (is.null(comp)) return(NULL)
        lvl <- input$cell_level
        if (!lvl %in% names(comp)) {
          shiny::showNotification("That cell group is not present in the object.",
                                  type = "error", duration = 10)
          return(NULL)
        }
        clin$.frac <- comp[[lvl]][match(clin$.id, comp$.id)]
        matched <- sum(!is.na(clin$.frac))
        if (matched < 3) {
          shiny::showNotification(
            paste("Only", matched, "sample id(s) matched between the clinical",
                  "table and the object's sample column. Check that they use",
                  "the same identifiers."),
            type = "error", duration = 15)
          return(NULL)
        }
        clin <- clin[!is.na(clin$.frac), , drop = FALSE]
        clin$.group <- split_numeric(clin$.frac, input$cut %||% "median",
                                     clin$.time, clin$.event)
        label <- sprintf("%s fraction (%s split)", lvl, input$cut %||% "median")
      }

      fit <- tryCatch(km_fit(clin), error = function(e) {
        shiny::showNotification(paste("Survival fit failed:", conditionMessage(e)),
                                type = "error", duration = 12); NULL })
      if (is.null(fit)) return(NULL)

      res$df      <- clin
      res$fit     <- fit
      res$lr      <- logrank_test(clin)
      res$cox     <- if (length(input$cox_vars)) cox_univariable(clin, input$cox_vars)
                     else NULL
      res$label   <- label
      res$matched <- matched

      # Share the cohort with every other omics pipeline.
      rv$clinical <- clin
      mark_done(rv, "clinical")
      log_step(log_rv, "Clinical & survival",
               params = list(id = input$id_col, time = input$time_col,
                             event = input$event_col, unit = input$time_unit,
                             stratify = label,
                             covariates = paste(input$cox_vars, collapse = ", ")),
               code = c(
                 sprintf('clin <- read.csv("%s")', input$file$name %||% "clinical.csv"),
                 'fit  <- survival::survfit(survival::Surv(time, event) ~ group, data = clin)',
                 'survival::survdiff(survival::Surv(time, event) ~ group, data = clin)'))
      shiny::showNotification("Survival analysis complete.", type = "message")
    })

    # ---- 4. outputs ---------------------------------------------------------
    output$summary <- shiny::renderUI({
      d <- res$df
      if (is.null(d)) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("Upload a clinical table, choose a grouping, then click <b>Run survival analysis</b>.",
                               "上传临床表格，选择分组方式，然后点击<b>运行生存分析</b>。")))
      }
      med <- km_medians(res$fit)
      p <- res$lr$p
      bslib::layout_columns(
        col_widths = c(3, 3, 3, 3),
        stat_tile(i18n("Patients", "患者数"), format(nrow(d), big.mark = ",")),
        stat_tile(i18n("Events", "事件数"), format(sum(d$.event), big.mark = ",")),
        stat_tile(i18n("Median OS (mo)", "中位生存(月)"),
                  if (nrow(med)) paste(round(med$median, 1), collapse = " / ") else "-"),
        stat_tile(i18n("Log-rank p", "Log-rank p"),
                  if (is.null(p)) "-" else signif(p, 3))
      )
    })

    output$km <- render_scop_plot(function() {
      shiny::req(res$fit)
      km_plot(res$fit, res$lr,
              title = paste0("Overall survival",
                             if (!is.null(res$label)) paste0(" — ", res$label) else ""))
    })

    output$cox_slot <- shiny::renderUI({
      if (is.null(res$cox) || !nrow(res$cox)) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("Pick one or more Cox covariates, then run.",
                               "选择一个或多个 Cox 协变量后运行。")))
      }
      if (has_pkg("DT")) DT::dataTableOutput(ns("cox_tbl"))
      else shiny::verbatimTextOutput(ns("cox_txt"))
    })

    output$cohort_slot <- shiny::renderUI({
      shiny::req(res$df)
      if (has_pkg("DT")) DT::dataTableOutput(ns("cohort_tbl"))
      else shiny::verbatimTextOutput(ns("cohort_txt"))
    })

    cohort_view <- shiny::reactive({
      d <- res$df
      shiny::req(d)
      keep <- c(".id", ".time", ".event", ".frac", ".group")
      out <- d[, intersect(keep, names(d)), drop = FALSE]
      names(out) <- sub("^\\.", "", names(out))
      if ("time" %in% names(out)) out$time <- round(out$time, 1)
      if ("frac" %in% names(out)) out$frac <- round(out$frac, 4)
      out
    })

    cox_view <- shiny::reactive({
      d <- res$cox
      shiny::req(d)
      d$HR <- round(d$HR, 3); d$lower <- round(d$lower, 3); d$upper <- round(d$upper, 3)
      d$p  <- signif(d$p, 3)
      d
    })

    if (has_pkg("DT")) {
      output$cox_tbl <- DT::renderDataTable({
        DT::datatable(cox_view(), rownames = FALSE,
                      options = list(pageLength = 10, scrollX = TRUE))
      })
      output$cohort_tbl <- DT::renderDataTable({
        DT::datatable(cohort_view(), rownames = FALSE, filter = "top",
                      options = list(pageLength = 15, scrollX = TRUE))
      })
    } else {
      output$cox_txt    <- shiny::renderPrint({ utils::head(cox_view(), 20) })
      output$cohort_txt <- shiny::renderPrint({ utils::head(cohort_view(), 20) })
    }
  })
}
