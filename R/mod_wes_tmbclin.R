#' Module: TMB versus outcome
#'
#' Is tumour mutational burden associated with survival or with response? TMB
#' is used as a continuous variable on a log2 scale (no cut-off to choose): a
#' Cox hazard ratio or a logistic odds ratio per doubling of TMB, plus the ROC
#' AUC for a binary response. The median split is drawn as an illustration
#' only. Sequenced patients with no variant record have TMB 0, not missing.
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_wes_tmbclin
NULL

#' @rdname mod_wes_tmbclin
#' @keywords internal
mod_wes_tmbclin_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = list(en = "TMB vs outcome", zh = "TMB 与临床结局"),
    what = list(
      en = "Test whether patients with a higher mutational burden live longer (or shorter), or respond
            more often to treatment.",
      zh = "检验突变负荷更高的患者生存更长（或更短），或者对治疗的应答率更高。"),
    why = list(
      en = "TMB is a candidate predictor of response to immune checkpoint inhibitors. Whether it carries
            information in your cohort is a question for the data, and the honest way to ask it is with
            TMB as a continuous variable: a cut-off picked after looking at the data inflates the
            evidence, and published cut-offs (10 mut/Mb) come from panel assays, not from WES.",
      zh = "TMB 是免疫检查点抑制剂疗效的候选预测指标。它在你的队列里是否有信息量，要由数据回答；诚实的问法是把 TMB 当作连续变量：看过数据后再挑的阈值会夸大证据，而文献中的阈值（10 mut/Mb）来自 panel 检测，不是 WES。"),
    how = list(
      en = "Enter the <b>capture size</b>, then pick the outcome from the clinical table imported with the
            MAF: <b>survival</b> (follow-up time and event) or a <b>binary response</b> (a column and the
            level that means responder). TMB enters the models as log2(TMB + 1 mutation / capture), so a
            hazard or odds ratio reads as the change per doubling of TMB.",
      zh = "填写<b>捕获区域大小</b>，再从随 MAF 导入的临床表中选择结局：<b>生存</b>（随访时间和事件）或<b>二分类疗效</b>（一列以及代表应答的取值）。TMB 以 log2(TMB + 1 个突变 / 捕获大小) 进入模型，因此风险比或比值比表示 TMB 每翻一倍的变化。"),
    read = list(
      en = "The headline is the HR or OR per doubling with its 95% CI; an interval that spans 1 means no
            clear association. For survival, the Kaplan-Meier plot splits at the median TMB only to
            illustrate. For response, the box plot shows TMB by group (rank-sum test) and the ROC
            curve how well TMB alone separates the groups (AUC 0.5 = no better than chance).",
      zh = "主要结果是每翻一倍的 HR 或 OR 及其 95% CI；区间跨过 1 表示没有明确的关联。生存分析中，Kaplan-Meier 图按 TMB 中位数分组，只用于直观展示。疗效分析中，箱线图按组显示 TMB（秩和检验），ROC 曲线显示仅凭 TMB 能把两组分开到什么程度（AUC 0.5 = 与随机无异）。"),
    example = list(
      en = "In TCGA-LAML (188 patients with follow-up, 35.8 Mb), the HR per doubling of TMB is about 0.99
            (95% CI 0.83–1.18): no association, as expected in a low-burden leukaemia.",
      zh = "在 TCGA-LAML 中（188 名有随访的患者，35.8 Mb），TMB 每翻一倍的 HR 约为 0.99（95% CI 0.83–1.18）：没有关联，这符合低突变负荷白血病的预期。")
  )
  controls <- shiny::tagList(
    label_with_help("Capture size (Mb) — required",
                    "The target territory of the capture kit, e.g. 35.8 = TCGA / Agilent SureSelect.",
                    "捕获区域大小（Mb）——必填", "捕获试剂盒的目标区域大小，例如 35.8 = TCGA / Agilent SureSelect。"),
    shiny::numericInput(ns("capture"), NULL, value = NA, min = 0.1, max = 3200, step = 0.1),
    label_with_help("Outcome", "Survival (time + event) or a binary response column.",
                    "结局", "生存（时间 + 事件）或二分类疗效列。"),
    shiny::radioButtons(ns("outcome"), NULL,
                        choiceNames = list(i18n("Survival", "生存"), i18n("Binary response", "二分类疗效")),
                        choiceValues = c("survival", "response"), selected = "survival", inline = TRUE),
    shiny::uiOutput(ns("mapping_ui")),
    shiny::checkboxInput(ns("unmatched_zero"),
                         i18n("Clinical samples with no MAF record were sequenced: TMB = 0",
                              "临床表中无 MAF 记录的样本已测序：TMB = 0"), value = TRUE),
    run_button(ns("run"), "Test TMB association", "检验 TMB 关联")
  )
  step_container(id = id, 
    title     = list(en = "TMB vs outcome", zh = "TMB 与临床结局"),
    subtitle  = list(en = "Survival or response per doubling of TMB, with no cut-off to choose.",
                     zh = "以 TMB 每翻一倍衡量与生存或疗效的关联，无需选择阈值。"),
    explainer = explainer,
    controls  = controls,
    summary   = shiny::uiOutput(ns("summary")),
    preview   = shiny::tagList(
      shiny::uiOutput(ns("insight")),
      bslib::navset_card_tab(
        bslib::nav_panel(i18n("Main figure", "主图"),
          preview_plot_ui(ns("main"), download = TRUE, scene = "wes_tmbclin",
            guide = list(en = "Survival by TMB, or TMB by response, will be drawn here.",
                         zh = "运行后，这里将绘制按 TMB 分组的生存曲线，或按疗效分组的 TMB。"),
            caption = list(en = "Survival: Kaplan-Meier split at the median TMB (illustration; the primary result is the Cox HR per doubling). Response: TMB by group, Wilcoxon rank-sum test.",
                           zh = "生存：按 TMB 中位数分组的 Kaplan-Meier 曲线（仅作展示；主要结果是每翻一倍的 Cox HR）。疗效：按组显示 TMB，Wilcoxon 秩和检验。"))),
        bslib::nav_panel(i18n("ROC", "ROC"),
          preview_plot_ui(ns("roc"), download = TRUE,
            caption = list(en = "Binary response only: sensitivity vs 1 - specificity of TMB as the only predictor; dashed = chance.",
                           zh = "仅用于二分类疗效：以 TMB 作为唯一预测指标时的灵敏度对 1 - 特异度；虚线 = 随机。"))),
        bslib::nav_panel(i18n("Patients", "患者"), shiny::uiOutput(ns("tbl_slot")))
      )
    )
  )
}

#' @rdname mod_wes_tmbclin
#' @keywords internal
mod_wes_tmbclin_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    res <- step_results(rv, "wes", d = NULL, fit = NULL, outcome = NULL, positive = NULL,
                        km = NULL, warn = character(0))

    output$mapping_ui <- shiny::renderUI({
      if (is.null(rv$maf)) return(wes_no_maf())
      df <- rv$wes_clin_raw
      if (is.null(df) || !ncol(df)) {
        return(shiny::div(class = "omicone-status-empty",
                          i18n("No clinical table was imported with this MAF: attach one on the Import step.",
                               "导入该 MAF 时未附带临床表：请在导入步骤附上临床表。")))
      }
      cols <- names(df)
      pick <- function(cands, default = cols[1]) {
        hit <- cols[tolower(cols) %in% cands]
        if (length(hit)) hit[1] else default
      }
      id_ui <- wes_col_select(ns, "id_col", label = list(en = "Sample id", zh = "样本编号"),
                              tip = list(en = "Matched to the MAF's Tumor_Sample_Barcode after upper-casing and trimming.",
                                         zh = "转为大写并去除空格后，与 MAF 的 Tumor_Sample_Barcode 匹配。"),
                              choices = cols,
                              selected = keep_selected(shiny::isolate(input$id_col), cols, pick("tumor_sample_barcode")))
      if (identical(input$outcome, "response")) {
        rc <- keep_selected(shiny::isolate(input$resp_col), cols,
                            pick(c("response", "best_response", "orr", "recist", "responder", "benefit", "dcb")))
        return(shiny::tagList(id_ui,
          wes_col_select(ns, "resp_col", label = list(en = "Response column", zh = "疗效列"),
                         tip = list(en = "Any column with the response category per patient; pick the level that means responder below.",
                                    zh = "每位患者疗效类别所在的列；在下方选择代表应答的取值。"),
                         choices = cols, selected = rc),
          shiny::uiOutput(ns("positive_ui"))))
      }
      time_sel <- pick(c("days_to_last_followup", "os_months", "os_days", "os.time", "pfs_months", "time", "futime", "os"))
      event_sel <- pick(c("overall_survival_status", "os_status", "os.event", "vital_status", "status", "fustat", "pfs_status"))
      unit_guess <- if (grepl("month", time_sel, ignore.case = TRUE)) "months"
                    else if (grepl("year", time_sel, ignore.case = TRUE)) "years" else "days"
      shiny::tagList(id_ui,
        wes_col_select(ns, "time_col", label = list(en = "Follow-up time", zh = "随访时间"),
                       tip = list(en = "Time to event or last contact.", zh = "到终点事件或最后随访的时间。"),
                       choices = cols, selected = keep_selected(shiny::isolate(input$time_col), cols, time_sel)),
        shiny::selectInput(ns("time_unit"), i18n("Time unit", "时间单位"),
                           c("Days" = "days", "Months" = "months", "Years" = "years"),
                           selected = keep_selected(shiny::isolate(input$time_unit), c("days", "months", "years"), unit_guess)),
        wes_col_select(ns, "event_col", label = list(en = "Event", zh = "终点事件"),
                       tip = list(en = "1 / Dead / TRUE = the event happened.", zh = "1 / Dead / TRUE 表示事件发生。"),
                       choices = cols, selected = keep_selected(shiny::isolate(input$event_col), cols, event_sel)))
    })

    output$positive_ui <- shiny::renderUI({
      df <- rv$wes_clin_raw
      shiny::req(df, input$resp_col %in% names(df))
      lv <- sort(unique(as.character(df[[input$resp_col]])))
      lv <- lv[!is.na(lv) & nzchar(trimws(lv))]
      guess <- lv[grepl("^(cr|pr|cr/pr|r|resp|responder|yes|1|true|dcb|benefit)$", lv, ignore.case = TRUE)]
      shiny::selectInput(ns("positive"), i18n("Level meaning responder", "代表应答的取值"), choices = lv,
                         selected = keep_selected(shiny::isolate(input$positive), lv, if (length(guess)) guess[1] else lv[1]))
    })

    shiny::observeEvent(input$run, {
      shiny::req(rv$maf)
      if (!require_pkgs(c("maftools", "survival"), "TMB vs outcome")) return(NULL)
      cap <- num_input(input$capture, 0.001)
      if (is.na(cap)) {
        wes_notify("Enter the capture size in Mb (e.g. 35.8 = TCGA/SureSelect).",
                   "请填写捕获区域大小（Mb，例如 35.8 = TCGA/SureSelect）。", duration = 12)
        return(NULL)
      }
      clin <- rv$wes_clin_raw
      if (is.null(clin) || !nrow(clin)) {
        wes_notify("No clinical table: attach one on the Import step.", "没有临床表：请在导入步骤附上临床表。")
        return(NULL)
      }
      outcome <- if (identical(input$outcome, "response")) "response" else "survival"
      unmatched_zero <- isTRUE(input$unmatched_zero)
      id_col <- input$id_col
      p <- list(capture = cap, id_col = id_col, unmatched_zero = unmatched_zero, outcome = outcome)
      if (outcome == "survival") {
        shiny::req(input$time_col, input$event_col)
        p <- c(p, list(time_col = input$time_col, event_col = input$event_col, time_unit = input$time_unit %||% "days"))
      } else {
        shiny::req(input$resp_col, input$positive)
        p <- c(p, list(resp_col = input$resp_col, positive = input$positive))
      }
      out <- with_progress_notify({
        tmb <- wes_tmb_table(rv$maf, cap, universe = if (unmatched_zero) clin[[id_col]], sequenced = rv$wes_sequenced)
        d <- wes_tmb_assoc_data(clin, tmb, id_col, outcome, time_col = p$time_col, event_col = p$event_col,
                                time_unit = p$time_unit %||% "days", resp_col = p$resp_col, positive = p$positive,
                                unmatched_zero = unmatched_zero)
        if (nrow(d) < 10) stop(sprintf("only %d patients with a usable outcome; check the id and outcome columns", nrow(d)))
        fit <- if (outcome == "survival") wes_tmb_cox(d) else wes_tmb_logit(d)
        if (is.null(fit)) stop("the model could not be fitted (TMB is constant, or too few patients)")
        km <- NULL
        if (outcome == "survival") {
          dk <- d
          dk$.group <- split_numeric(dk$tmb, "median")
          levels(dk$.group) <- c("TMB low (<= median)", "TMB high (> median)")
          if (length(unique(dk$.group)) == 2) km <- list(d = dk, fit = km_fit(dk), lr = logrank_test(dk))
        }
        list(d = d, fit = fit, km = km)
      }, message = "Testing TMB against the outcome...")
      if (is.null(out)) return(NULL)
      res$d <- out$d
      res$fit <- out$fit
      res$km <- out$km
      res$outcome <- outcome
      res$positive <- p$positive
      res$warn <- wes_tmb_assoc_warnings(out$d, outcome)
      mark_done(rv, "wes_tmbclin")
      fl <- attr(out$d, "flow")
      if (outcome == "survival") p$event_code <- wes_event_code(clin[[p$event_col]], p$event_col)
      log_step(log_rv, "WES TMB vs outcome",
               params = c(p[setdiff(names(p), "event_code")], list(n = fl$n, events = fl$events)),
               code = wes_tmb_assoc_code(p))
    })

    output$summary <- shiny::renderUI({
      if (is.null(rv$maf)) return(wes_no_maf())
      if (is.null(res$fit)) {
        return(wes_prompt("Enter the capture size, choose the outcome and click <b>Test TMB association</b>.",
                          "填写捕获区域大小、选择结局后点击<b>检验 TMB 关联</b>。"))
      }
      f <- res$fit
      fl <- attr(res$d, "flow")
      if (res$outcome == "survival") {
        return(shiny::tagList(
          stat_tile(i18n("Patients (events)", "患者数（事件数）"), sprintf("%s (%s)", wes_fmt(fl$n), wes_fmt(fl$events))),
          stat_tile(i18n("HR per doubling (95% CI)", "每翻一倍 HR（95% CI）"), sprintf("%.2f (%.2f–%.2f)", f$hr, f$lower, f$upper)),
          stat_tile(i18n("Wald p", "Wald p"), format(signif(f$p, 3))),
          stat_tile(i18n("PH check (cox.zph p)", "比例风险检验（cox.zph p）"),
                    if (is.na(f$zph_p)) "-" else format(signif(f$zph_p, 2)))))
      }
      shiny::tagList(
        stat_tile(i18n("Responders / others", "应答 / 其他"), sprintf("%d / %d", f$n_resp, f$n_non)),
        stat_tile(i18n("OR per doubling (95% CI)", "每翻一倍 OR（95% CI）"), sprintf("%.2f (%.2f–%.2f)", f$or, f$lower, f$upper)),
        stat_tile(i18n("Logistic p", "Logistic p"), format(signif(f$p, 3))),
        stat_tile(i18n("AUC (95% CI)", "AUC（95% CI）"), sprintf("%.2f (%.2f–%.2f)", f$auc, f$auc_lower, f$auc_upper)))
    })

    output$insight <- shiny::renderUI({
      f <- res$fit
      if (is.null(f)) return(NULL)
      fl <- attr(res$d, "flow")
      warn <- res$warn
      flow_en <- sprintf(" %d clinical sample(s) had no MAF record%s.", fl$unmatched,
                         if (fl$unmatched_zero) " and enter with TMB 0" else " and are left out")
      flow_zh <- sprintf("%d 个临床样本没有 MAF 记录%s。", fl$unmatched, if (fl$unmatched_zero) "，按 TMB 0 纳入" else "，未纳入")
      w_en <- if (length(warn)) sprintf(" <b>Caution:</b> %s.", paste(warn, collapse = "; ")) else ""
      w_zh <- if (length(warn)) sprintf("<b>注意：</b>%s。", paste(warn, collapse = "；")) else ""
      clear <- function(lo, hi) lo > 1 || hi < 1
      if (res$outcome == "survival") {
        insight_bar(
          sprintf("HR per doubling of TMB %.2f (95%% CI %.2f–%.2f, p = %s; %d patients, %d events): %s.%s%s",
                  f$hr, f$lower, f$upper, format(signif(f$p, 3)), f$n, f$events,
                  if (clear(f$lower, f$upper)) (if (f$hr < 1) "higher TMB goes with longer survival" else "higher TMB goes with shorter survival")
                  else "no clear association (the interval spans 1)", flow_en, w_en),
          sprintf("TMB 每翻一倍的 HR 为 %.2f（95%% CI %.2f–%.2f，p = %s；%d 名患者，%d 个事件）：%s。%s%s",
                  f$hr, f$lower, f$upper, format(signif(f$p, 3)), f$n, f$events,
                  if (clear(f$lower, f$upper)) (if (f$hr < 1) "TMB 越高生存越长" else "TMB 越高生存越短")
                  else "没有明确关联（区间跨过 1）", flow_zh, w_zh))
      } else {
        insight_bar(
          sprintf("OR per doubling of TMB %.2f (95%% CI %.2f–%.2f, p = %s); AUC %.2f (%.2f–%.2f); median TMB %.2f vs %.2f mut/Mb in %s vs others: %s.%s%s",
                  f$or, f$lower, f$upper, format(signif(f$p, 3)), f$auc, f$auc_lower, f$auc_upper,
                  f$median_resp, f$median_non, res$positive,
                  if (clear(f$lower, f$upper)) "TMB is associated with response" else "no clear association",
                  flow_en, w_en),
          sprintf("TMB 每翻一倍的 OR 为 %.2f（95%% CI %.2f–%.2f，p = %s）；AUC %.2f（%.2f–%.2f）；%s 组与其他组的 TMB 中位数为 %.2f 对 %.2f mut/Mb：%s。%s%s",
                  f$or, f$lower, f$upper, format(signif(f$p, 3)), f$auc, f$auc_lower, f$auc_upper,
                  res$positive, f$median_resp, f$median_non,
                  if (clear(f$lower, f$upper)) "TMB 与疗效相关" else "没有明确关联", flow_zh, w_zh))
      }
    })

    render_step_plot(output, input, "main", function() {
      shiny::req(res$fit)
      if (res$outcome == "survival") {
        shiny::validate(shiny::need(!is.null(res$km), "All patients have the same TMB group; no split to draw."))
        f <- res$fit
        km_plot(res$km$fit, res$km$lr, title = "Survival by TMB (median split, illustration)",
                note = sprintf("Primary result: Cox HR per doubling of TMB %.2f (%.2f-%.2f), p = %s",
                               f$hr, f$lower, f$upper, format(signif(f$p, 3))))
      } else {
        wes_tmb_resp_plot(res$d, res$fit, res$positive)
      }
    }, name = "wes_tmb_outcome", width = 8, height = 7)

    render_step_plot(output, input, "roc", function() {
      shiny::req(res$fit)
      shiny::validate(shiny::need(identical(res$outcome, "response"), "The ROC curve applies to a binary response."))
      wes_tmb_roc_plot(res$fit)
    }, name = "wes_tmb_roc", width = 6.5, height = 6.5)

    view <- shiny::reactive({
      d <- res$d
      shiny::req(d)
      out <- data.frame(sample = d$.id, TMB_per_Mb = round(d$tmb, 4), log2_TMB = round(d$x, 4),
                        stringsAsFactors = FALSE)
      if (identical(res$outcome, "survival")) {
        out$time_months <- round(d$.time, 2)
        out$event <- d$.event
      } else {
        out$responder <- d$.resp
      }
      out
    })
    tb <- wes_table(ns, "tbl", view, "wes_tmb_outcome",
                    has_maf = function() !is.null(rv$maf),
                    ready = function() !is.null(res$d))
    output$tbl_slot <- tb$slot
    output$tbl <- tb$table
    output$tbl_dl <- tb$download
  })
}
