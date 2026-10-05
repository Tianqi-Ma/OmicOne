#' Module: Clinical data & survival
#'
#' Load a per-patient clinical table (follow-up time + outcome), then ask
#' whether anything separates patients' survival. Two ways to stratify: by a
#' column of the clinical table itself (stage, treatment, a score you already
#' have), or by the composition of the working object -- the fraction of a
#' chosen cluster or cell type in each *patient*, computed after the cells are
#' filtered and the patients with too few cells are dropped.
#'
#' The order is fixed (AGENTS.md section 7): cells are aggregated per patient,
#' joined to the clinical table, missing values removed, and only then is any
#' cut-off derived. The normalised full cohort (no analysis columns) is written
#' to `rv$clinical`, the slot every omics shares; the analysis subset stays in
#' this module. The maths lives in `fct_survival.R`.
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
      en = "Linking what you see in the cells to outcome is how an atlas
            generates clinical hypotheses. Single-cell cohorts are small, so treat
            a result here as hypothesis-generating and validate it in a larger
            bulk cohort.",
      zh = "把细胞层面的发现与预后联系起来，是从图谱中提出临床假设的方式。单细胞队列通常很小，这里的结果应视为假设生成，并在更大的 bulk 队列中验证。"),
    how  = list(
      en = "Upload a table with <b>one row per patient</b> and map the id, time and
            outcome columns. To stratify by <b>cell composition</b>, pick the
            patient column of your object, the cell type, and optionally restrict
            the cells (e.g. tumour tissue at baseline) and the denominator (e.g.
            fraction <i>of T cells</i>). The continuous hazard ratio per SD is the
            primary result; the high/low split is an illustration.",
      zh = "上传<b>每位患者一行</b>的表格，并指定 ID、时间和终点列。若按<b>细胞组成</b>分组：选择对象中的患者列和细胞类型，可选地限定细胞（如基线时的肿瘤组织）和分母（如占 <i>T 细胞</i>的比例）。每个标准差的连续 HR 是主要结果，高低分组只是示意。"),
    read = list(
      en = "Curves show survival with 95% confidence bands; the table under the
            plot is the number still at risk. The log-rank p tests the curves; the
            HR (95% CI) is the effect size. With the <b>optimal</b> cut-point,
            report the selection-adjusted p, not the plain log-rank p.",
      zh = "曲线为生存率及其 95% 置信带，图下的表格是各时间点仍处于风险中的人数。log-rank p 检验曲线差异，HR（95% CI）是效应大小。使用<b>最优切点</b>时，应报告校正了切点选择的 p 值，而不是普通的 log-rank p。"),
    example = list(
      en = "Fraction of exhausted CD8 T cells among all T cells in baseline tumour
               samples, patients with at least 50 T cells: HR 1.8 per SD
               (95% CI 1.1–2.9).",
      zh = "基线肿瘤样本中，耗竭 CD8 T 细胞占全部 T 细胞的比例（每位患者至少 50 个 T 细胞）：每个标准差 HR 1.8（95% CI 1.1–2.9）。")
  )

  controls <- shiny::tagList(
    label_with_help("Clinical table",
                    "CSV or TSV, one row per patient. Ids must match the patient column of your object to stratify by composition.",
                    label_zh = "临床表格",
                    tip_zh = "CSV 或 TSV，每位患者一行。若要按细胞组成分组，ID 需与对象中的患者列一致。"),
    shiny::fileInput(ns("file"), i18n("Choose file", "选择文件"),
                     accept = c(".csv", ".tsv", ".txt")),
    shiny::selectInput(ns("sep"), i18n("Separator", "分隔符"),
                       choices = c("Auto" = "auto", "Comma" = ",", "Tab" = "\t",
                                   "Semicolon" = ";"),
                       selected = "auto"),
    shiny::uiOutput(ns("mapping")),
    shiny::hr(),
    label_with_help("Stratify by",
                    "A column of the clinical table, or the per-patient fraction of a cell type / cluster taken from the working object.",
                    label_zh = "分组依据",
                    tip_zh = "可选临床表中的某一列，或取自当前对象的每位患者细胞类型/簇占比。"),
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
    subtitle  = list(en = "Patient-level survival: one value per patient, then compare outcome.",
                     zh = "患者层面的生存分析：每位患者一个数值，再比较预后。"),
    explainer = explainer,
    controls  = controls,
    summary   = shiny::uiOutput(ns("summary")),
    preview   = shiny::tagList(
      shiny::uiOutput(ns("insight")),
      bslib::navset_card_tab(
        bslib::nav_panel(i18n("Kaplan-Meier", "生存曲线"), preview_plot_ui(ns("km"), download = TRUE,
          guide = list(en = "Kaplan-Meier curves with 95% bands and the number at risk will be drawn here.",
                       zh = "运行后，这里将绘制带 95% 置信带和风险人数表的 Kaplan-Meier 曲线。"),
          caption = list(en = "Survival by group: 95% confidence bands, censoring ticks (+), number at risk below.",
                         zh = "各组生存曲线：95% 置信带，删失标记（+），下方为风险人数。"))),
        bslib::nav_panel(i18n("Cox models", "Cox 模型"), shiny::uiOutput(ns("cox_slot"))),
        bslib::nav_panel(i18n("Cohort", "队列表"), shiny::uiOutput(ns("cohort_slot")))
      )
    )
  )
}

#' @rdname mod_clinical
#' @keywords internal
mod_clinical_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    ns  <- session$ns
    res <- step_results(rv, "sc", df = NULL, fit = NULL, lr = NULL, hr = NULL,
                        cont = NULL, ph = NA_real_, warn = character(0),
                        cox = NULL, mv = NULL, label = NULL, note = NULL,
                        endpoint = NULL, cut = NULL, matched = NA_integer_,
                        excluded = character(0))

    # ---- 1. read the uploaded table (re-read when the separator changes) ----
    raw <- shiny::reactive({
      shiny::req(input$file)
      tryCatch(read_clinical_table(input$file$datapath, sep = input$sep %||% "auto",
                                   name = input$file$name),
               error = function(e) {
                 shiny::showNotification(paste("Could not read the table:",
                                               conditionMessage(e)),
                                         type = "error", duration = 12)
                 NULL
               })
    })

    shiny::observeEvent(raw(), {
      df <- raw()
      shiny::req(df)
      shiny::showNotification(
        i18n(sprintf("Clinical table loaded: %d rows, %d columns.", nrow(df), ncol(df)),
             sprintf("临床表已加载：%d 行，%d 列。", nrow(df), ncol(df))),
        type = "message")
    })

    # Column mappers, populated from the uploaded table. Sensible guesses first
    # so the common case is a single click.
    output$mapping <- shiny::renderUI({
      df <- if (is.null(input$file)) NULL else raw()
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
        label_with_help("Patient id", "Column that identifies each patient.",
                        label_zh = "患者 ID", tip_zh = "标识每位患者的列。"),
        shiny::selectInput(ns("id_col"), NULL, cols,
                           selected = pick(c("patient", "patient_id", "sample", "sample_id",
                                             "id", "case_id", "bcr_patient_barcode"))),
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
        shiny::uiOutput(ns("event_level")),
        shiny::textInput(ns("endpoint"), i18n("Endpoint name", "终点名称"),
                         value = "Overall survival"),
        shiny::radioButtons(ns("dedup"), i18n("A patient listed twice", "同一患者出现多次时"),
                            c("Stop and tell me" = "error", "Keep the first row" = "first"),
                            selected = "error", inline = TRUE)
      )
    })

    needs_level <- shiny::reactive({
      df <- raw()
      !is.null(df) && !is.null(input$event_col) && input$event_col %in% names(df) &&
        all(is.na(encode_event(df[[input$event_col]])))
    })

    # Only ask which value means "event" when it cannot be inferred.
    output$event_level <- shiny::renderUI({
      shiny::req(needs_level())
      x <- raw()[[input$event_col]]
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
      df <- if (is.null(input$file)) NULL else raw()
      shiny::req(df)
      if (identical(input$mode, "clinical")) {
        cols <- setdiff(names(df), c(input$id_col, input$time_col, input$event_col))
        if (!length(cols)) {
          return(shiny::div(class = "omicone-status-empty",
                            i18n("No other column to group by.", "没有可用于分组的其他列。")))
        }
        return(shiny::tagList(
          shiny::selectInput(ns("clin_col"), i18n("Clinical column", "临床列"), cols,
                             selected = keep_selected(shiny::isolate(input$clin_col), cols)),
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
      cats <- categorical_cols(md, allow_na = TRUE)
      if (!length(cats)) {
        return(shiny::div(class = "omicone-status-empty",
                          i18n("The object has no categorical metadata to aggregate by.",
                               "对象中没有可用于汇总的分类元数据。")))
      }
      pat_default <- intersect(c("patient", "patient_id", "donor", "Patient"), cats)[1]
      if (is.na(pat_default)) pat_default <- guess_batch_col(md[cats]) %||% cats[1]
      cell_default <- intersect(c("celltype", "cell_type", "SingleR", "seurat_clusters"), cats)[1]
      if (is.na(cell_default)) cell_default <- cats[1]
      shiny::tagList(
        label_with_help("Patient column (in the object)",
                        "Which metadata column says which patient each cell came from. Use the patient, not the sample, when a patient has several samples.",
                        label_zh = "患者列（对象中）",
                        tip_zh = "对象元数据中标明每个细胞来自哪位患者的列。若一位患者有多个样本，请选患者列而不是样本列。"),
        shiny::selectInput(ns("pat_col"), NULL, cats,
                           selected = keep_selected(shiny::isolate(input$pat_col), cats, pat_default)),
        label_with_help("Cell grouping", "Clusters or annotated cell types.",
                        label_zh = "细胞分组", tip_zh = "聚类结果或已注释的细胞类型。"),
        shiny::selectInput(ns("cell_col"), NULL, cats,
                           selected = keep_selected(shiny::isolate(input$cell_col), cats, cell_default)),
        shiny::uiOutput(ns("cell_level_ui")),
        label_with_help("Use only cells where (optional)",
                        "Restrict the cells before counting, e.g. tissue = Tumor, timepoint = Baseline.",
                        label_zh = "仅使用满足条件的细胞（可选）",
                        tip_zh = "计数前先限定细胞，例如 tissue = Tumor、timepoint = Baseline。"),
        shiny::selectInput(ns("filter_col"), NULL, c("(all cells)" = "", cats),
                           selected = shiny::isolate(input$filter_col) %||% ""),
        shiny::uiOutput(ns("filter_val_ui")),
        shiny::numericInput(ns("min_cells"), i18n("Min. cells per patient", "每位患者最少细胞数"),
                            value = 50, min = 1, step = 10),
        shiny::uiOutput(ns("cut_ui"))
      )
    })

    output$cell_level_ui <- shiny::renderUI({
      md <- obj_meta(rv$obj)
      shiny::req(input$cell_col, input$cell_col %in% names(md))
      lv <- sort(unique(as.character(stats::na.omit(md[[input$cell_col]]))))
      shiny::tagList(
        shiny::selectInput(ns("cell_level"), i18n("Cell type / cluster", "细胞类型 / 簇"), lv,
                           selected = keep_selected(shiny::isolate(input$cell_level), lv)),
        shiny::selectizeInput(ns("denom"), i18n("Fraction of (denominator)", "占比的分母"),
                              choices = lv, multiple = TRUE,
                              options = list(placeholder = "all cells / 全部细胞"))
      )
    })

    output$filter_val_ui <- shiny::renderUI({
      md <- obj_meta(rv$obj)
      shiny::req(nzchar(input$filter_col %||% ""), input$filter_col %in% names(md))
      lv <- sort(unique(as.character(stats::na.omit(md[[input$filter_col]]))))
      shiny::selectizeInput(ns("filter_val"), NULL, lv, multiple = TRUE,
                            selected = lv[1])
    })

    # The cutpoint control only makes sense for a numeric stratifier.
    numeric_mode <- shiny::reactive({
      if (identical(input$mode, "composition")) return(TRUE)
      df <- raw()
      !is.null(df) && !is.null(input$clin_col) && input$clin_col %in% names(df) &&
        is.numeric(df[[input$clin_col]]) && !isTRUE(input$as_cat)
    })

    output$cut_ui <- shiny::renderUI({
      df <- raw()
      is_num_col <- identical(input$mode, "clinical") && !is.null(df) &&
        !is.null(input$clin_col) && input$clin_col %in% names(df) &&
        is.numeric(df[[input$clin_col]])
      shiny::tagList(
        if (is_num_col)
          shiny::checkboxInput(ns("as_cat"),
                               i18n("Treat as categories (e.g. stage coded 1–4)",
                                    "按分类变量处理（如编码为 1–4 的分期）"),
                               value = isTRUE(shiny::isolate(input$as_cat))),
        if (isTRUE(numeric_mode())) shiny::tagList(
          label_with_help("Cut-point for the curves",
                          "Median = balanced halves. Tertiles = top vs bottom third (middle dropped). Optimal = the split with the strongest separation: exploratory, reported with a selection-adjusted p-value.",
                          label_zh = "曲线的切点",
                          tip_zh = "中位数 = 两组样本量均衡。三分位 = 上三分之一对下三分之一（中间组舍弃）。最优 = 分离度最强的切点：仅供探索，同时报告校正了切点选择的 p 值。"),
          shiny::selectInput(ns("cut"), NULL,
                             c("Median" = "median", "Tertiles" = "tertile",
                               "Optimal (exploratory)" = "optimal"),
                             selected = shiny::isolate(input$cut) %||% "median"))
      )
    })

    output$covars <- shiny::renderUI({
      df <- if (is.null(input$file)) NULL else raw()
      shiny::req(df)
      cols <- setdiff(names(df), c(input$id_col, input$time_col, input$event_col))
      shiny::tagList(
        label_with_help("Cox covariates (optional)",
                        "Screened one by one (univariable, BH-adjusted) and, if ticked, together with the grouping in one adjusted model.",
                        label_zh = "Cox 协变量（可选）",
                        tip_zh = "逐个做单因素筛查（BH 校正）；若勾选，还会与分组一起放进一个校正后的多因素模型。"),
        shiny::selectizeInput(ns("cox_vars"), NULL, choices = cols, multiple = TRUE,
                              options = list(placeholder = "stage, age, ...")),
        shiny::checkboxInput(ns("adjusted"), i18n("Also fit the adjusted model",
                                                  "同时拟合多因素校正模型"), value = FALSE)
      )
    })

    # ---- 3. run -------------------------------------------------------------
    shiny::observeEvent(input$run, {
      df <- raw()
      shiny::req(df, input$id_col, input$time_col, input$event_col)
      mode <- input$mode
      cut_method <- input$cut %||% "median"
      endpoint <- trimws(input$endpoint %||% "")
      if (!nzchar(endpoint)) endpoint <- "Survival"
      positive <- if (isTRUE(needs_level())) input$event_positive else NULL

      clin <- tryCatch(
        normalise_clinical(df, input$id_col, input$time_col, input$event_col,
                           event_positive = positive,
                           time_unit = input$time_unit %||% "months",
                           dedup = input$dedup %||% "error"),
        error = function(e) {
          shiny::showNotification(conditionMessage(e), type = "error", duration = 15)
          NULL
        })
      if (is.null(clin)) return(NULL)
      if (!nrow(clin)) {
        shiny::showNotification(
          i18n("No usable rows: check the time and outcome columns.",
               "没有可用的行：请检查时间列和终点列。"),
          type = "error", duration = 12)
        return(NULL)
      }
      if (isTRUE(attr(clin, "dropped") > 0)) {
        shiny::showNotification(
          i18n(sprintf("%d row(s) dropped: missing or invalid time/outcome.", attr(clin, "dropped")),
               sprintf("已剔除 %d 行：时间或终点缺失/无效。", attr(clin, "dropped"))),
          type = "warning", duration = 10)
      }
      # the shared cohort is the full normalised table, never an analysis subset
      rv$clinical <- clin

      # Build the analysis set first; derive any cut-off from it afterwards.
      label <- NULL
      note <- NULL
      matched <- NA_integer_
      excluded <- character(0)
      value_col <- NULL
      comp_code <- character(0)
      if (identical(mode, "clinical")) {
        shiny::req(input$clin_col)
        v <- clin[[input$clin_col]]
        clin <- clin[!is.na(v) & nzchar(as.character(v)), , drop = FALSE]
        if (is.numeric(v) && !isTRUE(input$as_cat)) {
          value_col <- input$clin_col
          label <- input$clin_col
        } else {
          clin$.group <- factor(as.character(clin[[input$clin_col]]))
          label <- input$clin_col
        }
      } else {
        shiny::req(rv$obj, input$pat_col, input$cell_col, input$cell_level)
        md <- obj_meta(rv$obj)
        keep <- NULL
        if (nzchar(input$filter_col %||% "") && length(input$filter_val)) {
          keep <- as.character(md[[input$filter_col]]) %in% input$filter_val
        }
        min_cells <- max(1, as.integer(input$min_cells %||% 50), na.rm = TRUE)
        comp <- tryCatch(
          composition_by_patient(md, input$pat_col, input$cell_col, keep = keep,
                                 denom_levels = input$denom, min_cells = min_cells),
          error = function(e) {
            shiny::showNotification(conditionMessage(e), type = "error", duration = 12)
            NULL
          })
        if (is.null(comp)) return(NULL)
        excluded <- attr(comp, "excluded")
        lvl <- input$cell_level
        if (!lvl %in% names(comp)) {
          shiny::showNotification(
            i18n("That cell group has no cells after the filters.",
                 "按当前条件过滤后，该细胞分组没有细胞。"),
            type = "error", duration = 10)
          return(NULL)
        }
        clin$.frac <- comp[[lvl]][match(clin$.id, comp$.id)]
        matched <- sum(!is.na(clin$.frac))
        if (matched < 3) {
          shiny::showNotification(
            i18n(sprintf("Only %d patient id(s) matched between the clinical table and the object's patient column. Check that they use the same identifiers.", matched),
                 sprintf("临床表与对象的患者列只匹配上 %d 个 ID，请检查两者是否使用相同的标识。", matched)),
            type = "error", duration = 15)
          return(NULL)
        }
        clin <- clin[!is.na(clin$.frac), , drop = FALSE]
        value_col <- ".frac"
        denom_txt <- if (length(input$denom)) paste(input$denom, collapse = "+") else "all cells"
        label <- sprintf("%s fraction of %s", lvl, denom_txt)
        comp_code <- c(
          "md <- obj[[]]",
          if (!is.null(keep)) sprintf("md <- md[md[[%s]] %%in%% %s, , drop = FALSE]",
                                      deparse(input$filter_col), deparse(input$filter_val)),
          if (length(input$denom)) sprintf("md <- md[md[[%s]] %%in%% %s, , drop = FALSE]",
                                           deparse(input$cell_col), deparse(input$denom)),
          sprintf("tab <- table(md[[%s]], md[[%s]])", deparse(input$pat_col), deparse(input$cell_col)),
          sprintf("tab <- tab[rowSums(tab) >= %d, , drop = FALSE]", min_cells),
          sprintf("frac <- tab[, %s] / rowSums(tab)", deparse(lvl)),
          "clin$.frac <- unname(frac[match(clin$.id, rownames(tab))])",
          "clin <- clin[!is.na(clin$.frac), ]")
      }

      cut_value <- NULL
      used <- NULL
      if (!is.null(value_col)) {
        g <- tryCatch(split_numeric(clin[[value_col]], cut_method, clin$.time, clin$.event),
                      error = function(e) {
                        shiny::showNotification(conditionMessage(e), type = "error",
                                                duration = 12)
                        NULL
                      })
        if (is.null(g)) return(NULL)
        clin$.group <- g
        cut_value <- attr(g, "cutpoint")
        used <- attr(g, "method")
        label <- sprintf("%s (%s split)", label, used)
        if (identical(used, "optimal")) {
          note <- sprintf("Exploratory optimal cut-point %.3g; selection-adjusted p = %.3g",
                          cut_value, attr(g, "p_adjusted"))
        } else if (identical(used, "tertile")) {
          note <- "Middle third excluded from the curves"
        } else if (!identical(used, "median")) {
          note <- used
        }
      }

      analysed <- clin[!is.na(clin$.group), , drop = FALSE]
      if (length(unique(analysed$.group)) < 2) {
        shiny::showNotification(
          i18n("The grouping produced a single group, so there is nothing to compare.",
               "分组后只有一组，无法比较。"),
          type = "error", duration = 12)
        return(NULL)
      }

      fit <- tryCatch(km_fit(analysed), error = function(e) {
        shiny::showNotification(paste("Survival fit failed:", conditionMessage(e)),
                                type = "error", duration = 12)
        NULL
      })
      if (is.null(fit)) return(NULL)

      res$df       <- analysed
      res$fit      <- fit
      res$lr       <- logrank_test(analysed)
      res$hr       <- if (nlevels(droplevels(factor(analysed$.group))) == 2) cox_hr(analysed) else NULL
      res$cont     <- if (!is.null(value_col)) cox_continuous(clin, value_col) else NULL
      res$ph       <- if (!is.null(res$hr)) cox_ph_p(analysed) else NA_real_
      res$warn     <- survival_warnings(analysed)
      res$cox      <- if (length(input$cox_vars)) cox_univariable(clin, input$cox_vars) else NULL
      res$mv       <- if (isTRUE(input$adjusted) && length(input$cox_vars))
                        cox_multivariable(analysed, c(".group", input$cox_vars)) else NULL
      res$label    <- label
      res$note     <- note
      res$endpoint <- endpoint
      res$cut      <- cut_value
      res$matched  <- matched
      res$excluded <- excluded

      mark_done(rv, "clinical")
      log_step(log_rv, "Clinical & survival",
               params = list(id = input$id_col, time = input$time_col,
                             event = input$event_col, unit = input$time_unit,
                             endpoint = endpoint, stratify = label,
                             cutpoint = cut_value,
                             covariates = paste(input$cox_vars, collapse = ", ")),
               code = clinical_code(
                 file = input$file$name, sep = attr_sep(input), id = input$id_col,
                 time = input$time_col, event = input$event_col,
                 unit = input$time_unit %||% "months", positive = positive,
                 raw_event = df[[input$event_col]], comp_code = comp_code,
                 value_col = value_col, group_col = if (is.null(value_col)) input$clin_col,
                 cut_method = used %||% cut_method, cut_value = cut_value,
                 covars = input$cox_vars))
      shiny::showNotification(i18n("Survival analysis complete.", "生存分析完成。"),
                              type = "message")
    })

    # the separator actually used, for the log
    attr_sep <- function(input) {
      if (!identical(input$sep %||% "auto", "auto")) return(input$sep)
      guess_sep(input$file$datapath, input$file$name)
    }

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
      shiny::tagList(
        stat_tile(i18n("Patients (in curves)", "患者数（曲线中）"), format(nrow(d), big.mark = ",")),
        stat_tile(i18n("Events", "事件数"), format(sum(d$.event), big.mark = ",")),
        stat_tile(i18n("Median (mo)", "中位生存（月）"),
                  paste(sprintf("%s: %s", med$group, format_median(med$median)),
                        collapse = " · ")),
        stat_tile(i18n("Log-rank p", "Log-rank p"), if (is.null(p)) "-" else signif(p, 3)),
        if (!is.null(res$hr))
          stat_tile(i18n("HR", "HR"), sprintf("%.2f (%.2f–%.2f)", res$hr$hr,
                                               res$hr$lower, res$hr$upper))
      )
    })

    output$insight <- shiny::renderUI({
      d <- res$df
      shiny::req(d)
      parts_en <- character(0)
      parts_zh <- character(0)
      if (!is.null(res$cont)) {
        parts_en <- c(parts_en, sprintf("Continuous: HR %.2f per SD (95%% CI %.2f–%.2f), p = %.3g, n = %d.",
                                        res$cont$hr, res$cont$lower, res$cont$upper,
                                        res$cont$p, res$cont$n))
        parts_zh <- c(parts_zh, sprintf("连续变量：每个标准差 HR %.2f（95%% CI %.2f–%.2f），p = %.3g，n = %d。",
                                        res$cont$hr, res$cont$lower, res$cont$upper,
                                        res$cont$p, res$cont$n))
      }
      if (!is.null(res$hr)) {
        parts_en <- c(parts_en, sprintf("Curves: HR (%s) %.2f (95%% CI %.2f–%.2f).",
                                        res$hr$label, res$hr$hr, res$hr$lower, res$hr$upper))
        parts_zh <- c(parts_zh, sprintf("曲线：HR（%s）%.2f（95%% CI %.2f–%.2f）。",
                                        res$hr$label, res$hr$hr, res$hr$lower, res$hr$upper))
      }
      if (!is.null(res$note)) {
        parts_en <- c(parts_en, paste0(res$note, "."))
        parts_zh <- c(parts_zh, paste0(res$note, "。"))
      }
      if (is.finite(res$ph) && res$ph < 0.05) {
        parts_en <- c(parts_en, sprintf("Proportional hazards look violated (cox.zph p = %.3g): read the HR as an average.", res$ph))
        parts_zh <- c(parts_zh, sprintf("比例风险假设可能不成立（cox.zph p = %.3g）：HR 应理解为平均效应。", res$ph))
      }
      if (length(res$excluded)) {
        parts_en <- c(parts_en, sprintf("%d patient(s) excluded for having too few cells.", length(res$excluded)))
        parts_zh <- c(parts_zh, sprintf("%d 位患者因细胞数太少被排除。", length(res$excluded)))
      }
      if (length(res$warn)) {
        parts_en <- c(parts_en, paste("Small sample:", paste(res$warn, collapse = " "),
                                      "Treat as hypothesis-generating."))
        parts_zh <- c(parts_zh, paste0("样本量小：", paste(survival_warnings(d, "zh"), collapse = ""),
                                       "结果仅作为假设生成。"))
      }
      if (!length(parts_en)) return(NULL)
      insight_bar(paste(parts_en, collapse = " "), paste(parts_zh, collapse = ""))
    })

    render_step_plot(output, input, "km", function() {
      shiny::req(res$fit)
      km_plot(res$fit, res$lr,
              title = paste0(res$endpoint %||% "Survival",
                             if (!is.null(res$label)) paste0(" — ", res$label) else ""),
              hr = res$hr, note = res$note)
    }, name = "clinical_km", height = 8)

    output$cox_slot <- shiny::renderUI({
      if (is.null(res$cox) && is.null(res$mv)) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("Pick one or more Cox covariates, then run.",
                               "选择一个或多个 Cox 协变量后运行。")))
      }
      mv_note <- NULL
      if (!is.null(res$mv)) {
        epv <- attr(res$mv, "epv")
        mv_note <- shiny::div(
          class = "omicone-note",
          i18n(sprintf("Adjusted model: n = %d, %d events, %.1f events per coefficient%s.",
                       attr(res$mv, "n"), attr(res$mv, "events"), epv,
                       if (epv < 10) " (below 10: estimates unstable)" else ""),
               sprintf("多因素模型：n = %d，%d 个事件，每个系数 %.1f 个事件%s。",
                       attr(res$mv, "n"), attr(res$mv, "events"), epv,
                       if (epv < 10) "（低于 10：估计不稳定）" else "")))
      }
      shiny::tagList(
        shiny::h6(i18n("Univariable screen (BH across variables)", "单因素筛查（变量间 BH 校正）")),
        if (has_pkg("DT")) DT::dataTableOutput(ns("cox_tbl")) else shiny::verbatimTextOutput(ns("cox_txt")),
        if (!is.null(res$mv)) shiny::tagList(
          shiny::h6(class = "mt-3", i18n("Adjusted model (grouping + covariates)", "多因素模型（分组 + 协变量）")),
          mv_note,
          if (has_pkg("DT")) DT::dataTableOutput(ns("mv_tbl")) else shiny::verbatimTextOutput(ns("mv_txt")))
      )
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

    round_tbl <- function(d) {
      for (cl in intersect(c("HR", "lower", "upper"), names(d))) d[[cl]] <- round(d[[cl]], 3)
      for (cl in intersect(c("p", "p_overall", "q"), names(d))) d[[cl]] <- signif(d[[cl]], 3)
      d
    }
    cox_view <- shiny::reactive({
      shiny::req(res$cox)
      round_tbl(res$cox)
    })
    mv_view <- shiny::reactive({
      shiny::req(res$mv)
      round_tbl(as.data.frame(res$mv))
    })

    if (has_pkg("DT")) {
      output$cox_tbl <- DT::renderDataTable({
        DT::datatable(cox_view(), rownames = FALSE,
                      options = list(pageLength = 10, scrollX = TRUE))
      })
      output$mv_tbl <- DT::renderDataTable({
        DT::datatable(mv_view(), rownames = FALSE,
                      options = list(pageLength = 10, scrollX = TRUE, dom = "t"))
      })
      output$cohort_tbl <- DT::renderDataTable({
        DT::datatable(cohort_view(), rownames = FALSE, filter = "top",
                      options = list(pageLength = 15, scrollX = TRUE))
      })
    } else {
      output$cox_txt    <- shiny::renderPrint(utils::head(cox_view(), 20))
      output$mv_txt     <- shiny::renderPrint(mv_view())
      output$cohort_txt <- shiny::renderPrint(utils::head(cohort_view(), 20))
    }
  })
}

#' Runnable R code reproducing a clinical survival run
#'
#' Base R + survival only (no OmicOne internals), so the logged script runs
#' anywhere: read, encode time/event exactly as [normalise_clinical()] did,
#' build the grouping with the cut-off value actually used, then fit.
#' @keywords internal
clinical_code <- function(file, sep, id, time, event, unit, positive, raw_event,
                          comp_code, value_col, group_col, cut_method, cut_value,
                          covars) {
  q <- function(x) deparse(x)
  # full precision: a rounded cut-off moves boundary patients between groups
  exact <- function(x) format(x, digits = 17)
  scale <- switch(unit, days = " / 30.4375", years = " * 12", "")
  event_line <- if (!is.null(positive)) {
    sprintf("clin$.event <- as.integer(as.character(clin[[%s]]) == %s)", q(event), q(positive))
  } else if (is.logical(raw_event)) {
    sprintf("clin$.event <- as.integer(clin[[%s]])", q(event))
  } else if (is.numeric(raw_event)) {
    u <- sort(unique(stats::na.omit(raw_event)))
    if (length(u) == 2 && all(u == c(1, 2)))
      sprintf("clin$.event <- as.integer(clin[[%s]] == 2)", q(event))
    else sprintf("clin$.event <- as.integer(clin[[%s]] > 0)", q(event))
  } else {
    s <- tolower(trimws(as.character(raw_event)))
    dead <- intersect(unique(s), c("1", "true", "yes", "y", "dead", "deceased", "death",
                                   "event", "progressed", "progression", "recurrence", "relapse"))
    alive <- intersect(unique(s), c("0", "false", "no", "n", "alive", "living", "censored",
                                    "censor", "no event", "disease free", "disease-free"))
    c(sprintf("ev <- tolower(trimws(clin[[%s]]))", q(event)),
      sprintf("clin$.event <- ifelse(ev %%in%% %s, 1L, ifelse(ev %%in%% %s, 0L, NA))",
              q(dead), q(alive)))
  }
  group_lines <- if (!is.null(value_col)) {
    v <- if (value_col == ".frac") "clin$.frac" else sprintf("clin[[%s]]", q(value_col))
    c(sprintf("clin <- clin[!is.na(%s), ]", v),
      if (identical(cut_method, "tertile") && length(cut_value) == 2)
        sprintf('clin$.group <- factor(ifelse(%1$s <= %2$s, "Low", ifelse(%1$s >= %3$s, "High", NA)), levels = c("Low", "High"))',
                v, exact(cut_value[1]), exact(cut_value[2]))
      else
        sprintf('clin$.group <- factor(ifelse(%s > %s, "High", "Low"), levels = c("Low", "High"))  # %s cut-point',
                v, exact(cut_value), cut_method))
  } else {
    c(sprintf("clin <- clin[!is.na(clin[[%s]]), ]", q(group_col)),
      sprintf("clin$.group <- factor(clin[[%s]])", q(group_col)))
  }
  c(sprintf("clin <- read.delim(%s, sep = %s, check.names = FALSE)", q(file %||% "clinical.csv"), q(sep)),
    sprintf("clin$.id <- trimws(as.character(clin[[%s]]))", q(id)),
    sprintf("clin$.time <- as.numeric(clin[[%s]])%s   # months", q(time), scale),
    event_line,
    "clin <- clin[!is.na(clin$.time) & clin$.time >= 0 & !is.na(clin$.event), ]",
    comp_code,
    group_lines,
    "fit <- survival::survfit(survival::Surv(.time, .event) ~ .group, data = clin)",
    "survival::survdiff(survival::Surv(.time, .event) ~ .group, data = clin)",
    "summary(survival::coxph(survival::Surv(.time, .event) ~ .group, data = clin))",
    if (length(covars))
      sprintf("summary(survival::coxph(survival::Surv(.time, .event) ~ .group + %s, data = clin))",
              paste(sprintf("`%s`", covars), collapse = " + ")))
}
