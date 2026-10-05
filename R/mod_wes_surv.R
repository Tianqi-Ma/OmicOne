#' WES module 11: Mutation vs survival
#'
#' Does carrying a mutation in a gene (or any gene of a set) change outcome? This
#' step reuses the shared survival layer (`fct_survival.R`) for the curves, the
#' log-rank test and the Cox hazard ratio, so they match the single-cell side.
#' The analysis set is built in [wes_surv_data()]: ids are normalised, samples
#' are collapsed to one row per patient, and — by default — clinical rows with
#' no MAF record are kept as sequenced wild-type, because `read.maf()` drops
#' exactly those patients and the WT arm would otherwise lose them. The shared
#' `rv$clinical` is only read, never written.
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_wes_surv
NULL

#' @rdname mod_wes_surv
#' @keywords internal
mod_wes_surv_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = list(en = "Mutation vs survival", zh = "突变与预后"),
    what = list(
      en = "Split patients into mutant and wild-type for a gene (or a set of
            genes) and compare their survival with Kaplan-Meier curves, the
            log-rank test and a Cox hazard ratio (mutant vs WT).",
      zh = "按某个基因（或一组基因）把患者分为突变型与野生型，用 Kaplan-Meier 曲线、log-rank 检验和 Cox 风险比（突变型 vs 野生型）比较两者的生存差异。"),
    why  = list(
      en = "Frequency tells you a gene is mutated; survival tells you whether that
            matters to the patient. The comparison is only fair if every
            sequenced patient is in it: a MAF lists samples with variants, so
            patients with none are missing from it — they belong to the WT arm.",
      zh = "频率只能说明某个基因发生了突变；生存分析才能说明这对患者是否有意义。只有把每一位已测序的患者都纳入，比较才公平：MAF 只列出有变异的样本，没有任何变异的患者不在其中——他们属于野生型组。"),
    how  = list(
      en = "Pick one gene to start; several genes form a set (<i>mutant</i> =
            a non-synonymous mutation in <b>any</b> of them). Survival data come
            from the clinical table imported with the MAF (default) or from the
            cohort of the <b>Clinical &amp; survival</b> step. Keep <b>no MAF record
            = WT</b> ticked when the clinical table lists sequenced patients
            only; untick it if it also lists patients who were never sequenced.
            Ids are matched without regard to case or spaces; use the TCGA
            option or a patient-id column when one patient has several samples
            (mutant if any sample is).",
      zh = "先从单个基因开始；选择多个基因时作为基因集处理（<i>突变型</i>＝其中<b>任意一个</b>发生非同义突变）。生存数据来自随 MAF 导入的临床表（默认），或<b>临床与生存</b>步骤中的队列。若临床表只包含已测序患者，请保持勾选<b>无 MAF 记录＝野生型</b>；若其中也有从未测序的患者，请取消勾选。编号匹配不区分大小写与空格；一名患者有多个样本时，请使用 TCGA 选项或患者编号列（任一样本突变即算突变型）。"),
    read = list(
      en = "Each curve is the fraction of patients still event-free over time;
            every step down is an event and ticks are censored patients. The
            shaded band is the 95% CI, and the risk table under the plot shows
            how many patients each arm still has at each time. The log-rank p
            tests the whole curves; the HR (95% CI) comes from a Cox model with
            WT as reference. A separation resting on a handful of late patients
            is fragile whatever the p-value; small-sample warnings appear above
            the plot.",
      zh = "每条曲线是随时间推移仍未发生事件的患者比例；每一次下降是一个事件，刻度线为删失患者。阴影带为 95% 置信区间，图下的风险人数表显示各组在各时间点的剩余人数。log-rank p 检验整条曲线；HR（95% CI）来自以野生型为参照的 Cox 模型。只靠少数晚期患者支撑的曲线分离并不牢靠，无论 p 值多小；小样本警告显示在图上方。"),
    example = list(
      en = "In TCGA LAML (188 patients with follow-up, 7 without any MAF record
               counted as WT), <code>TP53</code> mutants do clearly worse
               (log-rank p ≈ 1e-5); <code>DNMT3A</code> mutants also do worse
               (p ≈ 0.001), while <code>FLT3</code> does not separate (p ≈ 0.14).",
      zh = "在 TCGA LAML 中（188 名有随访的患者，其中 7 名无任何 MAF 记录者计为野生型），<code>TP53</code> 突变型明显更差（log-rank p ≈ 1e-5）；<code>DNMT3A</code> 突变型同样更差（p ≈ 0.001），而 <code>FLT3</code> 没有分离（p ≈ 0.14）。")
  )
  controls <- shiny::tagList(
    shiny::uiOutput(ns("gene_ui")),
    label_with_help("Endpoint name", "Used in the plot title and the insight, e.g. Overall survival, Progression-free survival.",
                    label_zh = "终点名称", tip_zh = "用于图标题和结论栏，例如 Overall survival、Progression-free survival。"),
    shiny::textInput(ns("endpoint"), NULL, value = "Overall survival"),
    shiny::hr(),
    shiny::uiOutput(ns("source_ui")),
    shiny::uiOutput(ns("mapping_ui")),
    shiny::uiOutput(ns("options_ui")),
    run_button(ns("run"), "Run survival analysis", "运行生存分析")
  )
  step_container(
    title     = list(en = "Mutation vs survival", zh = "突变与预后"),
    subtitle  = list(en = "Kaplan-Meier curves, log-rank test and Cox HR for mutant versus wild-type.",
                     zh = "突变型与野生型的 Kaplan-Meier 曲线、log-rank 检验与 Cox 风险比。"),
    explainer = explainer,
    controls  = controls,
    summary   = shiny::uiOutput(ns("summary")),
    preview   = shiny::tagList(
      shiny::uiOutput(ns("insight")),
      bslib::navset_card_tab(
        bslib::nav_panel(i18n("Kaplan-Meier", "生存曲线"),
                         preview_plot_ui(ns("km"), download = TRUE,
                                         guide = list(en = "The Kaplan-Meier curves will be drawn here.",
                                                      zh = "运行后，这里将绘制 Kaplan-Meier 生存曲线。"),
                                         caption = list(en = "Curves: fraction event-free over time (months); band = 95% CI; ticks = censored; table below = patients still at risk. p = log-rank; HR = Cox, mutant vs WT.",
                                                        zh = "曲线：随时间（月）的无事件比例；阴影带＝95% CI；刻度＝删失；下表＝各时点风险人数。p＝log-rank；HR＝Cox，突变型 vs 野生型。"))),
        bslib::nav_panel(i18n("Cohort", "队列表"), shiny::uiOutput(ns("tbl_slot")))
      ))
  )
}

#' @rdname mod_wes_surv
#' @keywords internal
mod_wes_surv_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    ns  <- session$ns
    res <- step_results(rv, "wes", df = NULL, fit = NULL, lr = NULL, hr = NULL,
                        warn = character(0), label = NULL, endpoint = NULL, flow = NULL,
                        source = NULL)

    output$gene_ui <- shiny::renderUI({
      shiny::req(rv$maf)
      g <- wes_genes(rv$maf, n = 300)
      wes_col_select(ns, "genes",
                     label = list(en = "Gene(s)", zh = "基因"),
                     tip = list(en = "One gene, or several treated as a set (mutant = a non-synonymous mutation in any of them).",
                                zh = "可选单个基因，也可选多个作为基因集（突变型＝其中任一发生非同义突变）。"),
                     choices = g, selected = if (length(g)) g[1] else NULL, multiple = TRUE)
    })

    has_shared <- shiny::reactive(!is.null(rv$clinical) && nrow(rv$clinical) > 0)

    output$source_ui <- shiny::renderUI({
      shiny::req(rv$maf)
      choices <- c("Clinical table imported with the MAF" = "maf")
      if (has_shared()) choices <- c(choices, "Cohort from Clinical & survival" = "shared")
      shiny::tagList(
        label_with_help("Survival data from",
                        "The shared cohort is the table loaded in the Clinical & survival step; this step only reads it.",
                        label_zh = "生存数据来源",
                        tip_zh = "共享队列即「临床与生存」步骤中加载的表；本步骤只读取，不修改。"),
        shiny::radioButtons(ns("src"), NULL, choices,
                            selected = keep_selected(shiny::isolate(input$src), choices, "maf"))
      )
    })

    src_kind <- shiny::reactive(if (identical(input$src, "shared") && has_shared()) "shared"
                                else "maf")

    # the table the analysis set is built from, and its column mapping
    source_spec <- shiny::reactive({
      shiny::req(rv$maf)
      if (identical(src_kind(), "shared")) {
        return(list(df = rv$clinical, id = ".id", time = ".time", event = ".event",
                    unit = "months", patient = NULL))
      }
      df <- rv$wes_clin_raw
      shiny::req(df, input$id_col, input$time_col, input$event_col)
      shiny::req(all(c(input$id_col, input$time_col, input$event_col) %in% names(df)))
      pcol <- input$patient_col
      list(df = df, id = input$id_col, time = input$time_col, event = input$event_col,
           unit = input$time_unit %||% "days",
           patient = if (!is.null(pcol) && nzchar(pcol) && pcol %in% names(df)) pcol else NULL)
    })

    overlap <- shiny::reactive({
      sp <- source_spec()
      tc <- isTRUE(input$tcga12)
      ids <- wes_norm_id(sp$df[[sp$id]], tc)
      ids <- unique(ids[!is.na(ids) & nzchar(ids)])
      maf_ids <- unique(wes_norm_id(wes_samples(rv$maf), tc))
      list(rows = nrow(sp$df), ids = length(ids), in_maf = sum(ids %in% maf_ids),
           no_maf = sum(!ids %in% maf_ids), maf = length(maf_ids))
    })

    output$mapping_ui <- shiny::renderUI({
      shiny::req(rv$maf)
      if (identical(src_kind(), "shared")) {
        ov <- tryCatch(overlap(), error = function(e) NULL)
        txt_en <- if (is.null(ov)) "Using the cohort from the Clinical &amp; survival step."
                  else sprintf("Shared cohort: %s rows; %s of its ids are MAF samples.",
                               wes_fmt(ov$rows), wes_fmt(ov$in_maf))
        txt_zh <- if (is.null(ov)) "正在使用「临床与生存」步骤中的队列。"
                  else sprintf("共享队列：%s 行；其中 %s 个编号是 MAF 样本。",
                               wes_fmt(ov$rows), wes_fmt(ov$in_maf))
        return(shiny::div(class = "omicone-status-empty", i18n(txt_en, txt_zh)))
      }
      df <- rv$wes_clin_raw
      if (is.null(df) || !ncol(df)) {
        return(shiny::div(class = "omicone-status-empty",
                          i18n("No clinical table was imported with this MAF. Attach one on the Import step, or load a cohort in the Clinical &amp; survival step.",
                               "导入该 MAF 时未附带临床表。请在导入步骤附上临床表，或在「临床与生存」步骤加载队列。")))
      }
      cols <- names(df)
      pick <- function(cands, default = cols[1]) {
        hit <- cols[tolower(cols) %in% cands]
        if (length(hit)) hit[1] else default
      }
      time_sel <- pick(c("days_to_last_followup", "os_months", "os_days", "os.time",
                         "overall_survival_time", "time", "futime", "os"))
      unit_guess <- if (grepl("month", time_sel, ignore.case = TRUE)) "months"
                    else if (grepl("year", time_sel, ignore.case = TRUE)) "years" else "days"
      event_sel <- pick(c("overall_survival_status", "os_status", "os.event", "vital_status",
                          "status", "fustat"))
      shiny::tagList(
        wes_col_select(ns, "id_col",
                       label = list(en = "Sample id", zh = "样本编号"),
                       tip = list(en = "Matched to the MAF's Tumor_Sample_Barcode after upper-casing and trimming.",
                                  zh = "转为大写并去除空格后，与 MAF 的 Tumor_Sample_Barcode 匹配。"),
                       choices = cols, selected = keep_selected(shiny::isolate(input$id_col), cols,
                                                                pick("tumor_sample_barcode"))),
        wes_col_select(ns, "time_col",
                       label = list(en = "Follow-up time", zh = "随访时间"),
                       tip = list(en = "Time to event or last contact. Missing or non-finite values (Inf) are dropped and counted.",
                                  zh = "到终点事件或最后随访的时间。缺失或非有限值（Inf）会被剔除并计数。"),
                       choices = cols, selected = keep_selected(shiny::isolate(input$time_col),
                                                                cols, time_sel)),
        shiny::selectInput(ns("time_unit"), i18n("Time unit", "时间单位"),
                           c("Days" = "days", "Months" = "months", "Years" = "years"),
                           selected = keep_selected(shiny::isolate(input$time_unit),
                                                    c("days", "months", "years"), unit_guess)),
        wes_col_select(ns, "event_col",
                       label = list(en = "Outcome", zh = "终点事件"),
                       tip = list(en = "1 / Dead / TRUE = the event happened.",
                                  zh = "1 / Dead / TRUE 表示事件发生。"),
                       choices = cols,
                       selected = keep_selected(shiny::isolate(input$event_col), cols, event_sel)),
        wes_col_select(ns, "patient_col",
                       label = list(en = "Patient id (optional)", zh = "患者编号（可选）"),
                       tip = list(en = "When several samples belong to one patient: one row per patient, mutant if any sample is. Follow-up must agree across a patient's rows.",
                                  zh = "一名患者有多个样本时使用：每名患者一行，任一样本突变即为突变型。同一患者各行的随访信息必须一致。"),
                       choices = cols, selected = keep_selected(shiny::isolate(input$patient_col),
                                                                c("", cols), ""),
                       none = "(none — one row per sample)")
      )
    })

    output$options_ui <- shiny::renderUI({
      shiny::req(rv$maf)
      ov <- tryCatch(overlap(), error = function(e) NULL)
      k <- if (is.null(ov)) "?" else wes_fmt(ov$no_maf)
      tcga_like <- mean(grepl("^TCGA-", wes_samples(rv$maf))) > 0.5
      shiny::tagList(
        shiny::checkboxInput(
          ns("unmatched_wt"),
          i18n(sprintf("The %s clinical sample(s) with no MAF record were sequenced: count them as WT", k),
               sprintf("临床表中无 MAF 记录的 %s 个样本已测序、视为 WT", k)),
          value = isTRUE(shiny::isolate(input$unmatched_wt) %||% TRUE)),
        shiny::div(class = "omicone-muted",
                   i18n("read.maf keeps only samples with at least one variant, so a sequenced patient with none is absent from the MAF; dropping them would bias the WT arm. Untick if the table also lists unsequenced patients.",
                        "read.maf 只保留至少有一个变异的样本，没有任何变异的已测序患者因此不在 MAF 中；剔除他们会使野生型组产生偏倚。若临床表中还有未测序的患者，请取消勾选。")),
        if (tcga_like) {
          shiny::checkboxInput(ns("tcga12"),
                               i18n("Match TCGA barcodes on the first 12 characters (patient)",
                                    "按 TCGA 条码前 12 位（患者）匹配"),
                               value = isTRUE(shiny::isolate(input$tcga12)))
        }
      )
    })

    shiny::observeEvent(input$run, {
      shiny::req(rv$maf, input$genes)
      if (!require_pkgs("maftools", "Mutation vs survival")) return(NULL)
      genes <- input$genes
      endpoint <- trimws(input$endpoint %||% "")
      if (!nzchar(endpoint)) endpoint <- "Overall survival"
      kind <- src_kind()
      sp <- tryCatch(source_spec(), error = function(e) NULL)
      if (is.null(sp) || is.null(sp$df) || !nrow(sp$df)) {
        wes_notify("No usable clinical table: attach one at import, or load a cohort in Clinical & survival.",
                   "没有可用的临床表：请在导入时附上临床表，或在「临床与生存」中加载队列。")
        return(NULL)
      }
      tcga12 <- isTRUE(input$tcga12)
      unmatched_wt <- isTRUE(input$unmatched_wt)

      st <- tryCatch(wes_mutation_status(rv$maf, genes,
                                         universe = if (unmatched_wt) sp$df[[sp$id]],
                                         tcga12 = tcga12),
                     error = function(e) {
                       wes_notify(conditionMessage(e), conditionMessage(e), duration = 12)
                       NULL
                     })
      if (is.null(st)) return(NULL)
      d <- tryCatch(wes_surv_data(sp$df, st, sp$id, sp$time, sp$event, time_unit = sp$unit,
                                  patient_col = sp$patient, tcga12 = tcga12,
                                  unmatched_wt = unmatched_wt),
                    error = function(e) {
                      wes_notify(conditionMessage(e), conditionMessage(e), duration = 15)
                      NULL
                    })
      if (is.null(d)) return(NULL)
      fl <- attr(d, "flow")
      if (nrow(d) < 3 || fl$matched < 3) {
        wes_notify(sprintf("Only %d patient(s) with follow-up (%d matched to the MAF). Check that both tables use the same sample ids.",
                           nrow(d), fl$matched),
                   sprintf("只有 %d 名有随访的患者（%d 名与 MAF 匹配）。请检查两张表是否使用相同的样本编号。",
                           nrow(d), fl$matched), duration = 15)
        return(NULL)
      }
      if (length(unique(d$.group)) < 2) {
        wes_notify("Every patient falls in one group; nothing to compare.",
                   "所有患者都落在同一组，无法比较。", duration = 12)
        return(NULL)
      }
      fit <- tryCatch(km_fit(d), error = function(e) {
        wes_notify(paste("Survival fit failed:", conditionMessage(e)),
                   paste("生存拟合失败：", conditionMessage(e)), duration = 12)
        NULL
      })
      if (is.null(fit)) return(NULL)
      lr <- logrank_test(d)
      hr <- tryCatch(cox_hr(d), error = function(e) NULL)
      warn <- tryCatch(survival_warnings(d), error = function(e) character(0))

      res$df       <- d
      res$fit      <- fit
      res$lr       <- lr
      res$hr       <- hr
      res$warn     <- warn
      res$label    <- paste(genes, collapse = " / ")
      res$endpoint <- endpoint
      res$flow     <- fl
      res$source   <- kind
      mark_done(rv, "wes_surv")
      log_step(log_rv, "WES mutation vs survival",
               params = list(genes = paste(genes, collapse = ", "), endpoint = endpoint,
                             source = kind, id = sp$id, time = sp$time, unit = sp$unit,
                             event = sp$event, patient = sp$patient %||% "(none)",
                             tcga12 = tcga12, unmatched_as_WT = unmatched_wt,
                             n = fl$n, events = fl$events),
               code = wes_surv_code(list(genes = genes, source = kind, id_col = sp$id,
                                         time_col = sp$time, event_col = sp$event,
                                         time_unit = sp$unit, patient_col = sp$patient,
                                         tcga12 = tcga12, unmatched_wt = unmatched_wt,
                                         raw_event = sp$df[[sp$event]])))
    })

    stats <- shiny::reactive({
      shiny::req(res$df, res$fit)
      med <- tryCatch(km_medians(res$fit), error = function(e) NULL)
      med_txt <- "-"
      if (!is.null(med) && nrow(med)) {
        med <- med[match(c("WT", "Mutant"), med$group), , drop = FALSE]
        med_txt <- paste(trimws(format_median(med$median)), collapse = " / ")
      }
      hr <- res$hr
      hr_txt <- if (!is.null(hr) && is.finite(hr$hr))
        sprintf("%.2f (%.2f–%.2f)", hr$hr, hr$lower, hr$upper) else "-"
      p <- res$lr$p
      list(flow = res$flow, med = med_txt, hr = hr_txt,
           hr_p = if (!is.null(hr)) hr$p else NA_real_,
           p = if (is.null(p) || !is.finite(p)) NA_real_ else p)
    })

    output$summary <- shiny::renderUI({
      if (is.null(rv$maf)) return(wes_no_maf())
      if (is.null(res$df)) {
        return(wes_prompt("Pick a gene and click <b>Run survival analysis</b>.",
                          "选择基因后点击<b>运行生存分析</b>。"))
      }
      s <- stats()
      shiny::tagList(
        stat_tile(i18n("Patients (events)", "患者数（事件数）"),
                  sprintf("%s (%s)", wes_fmt(s$flow$n), wes_fmt(s$flow$events))),
        stat_tile(i18n("Mutant / WT", "突变型 / 野生型"),
                  sprintf("%s / %s", wes_fmt(s$flow$n_mut), wes_fmt(s$flow$n_wt))),
        stat_tile(i18n("Median, mo (WT / Mut)", "中位生存，月（野生/突变）"), s$med),
        stat_tile(i18n("Log-rank p", "Log-rank p"),
                  if (is.na(s$p)) "-" else format(signif(s$p, 3))),
        stat_tile(i18n("HR mut vs WT (Cox)", "HR 突变 vs 野生（Cox）"), s$hr)
      )
    })

    output$insight <- shiny::renderUI({
      if (is.null(res$df)) return(NULL)
      s <- stats()
      fl <- s$flow
      flow_en <- sprintf("Clinical table %s rows (%s ids%s); MAF %s samples; %s matched; %s without a MAF record %s; %s dropped for missing time or event.",
                         wes_fmt(fl$clin_rows), wes_fmt(fl$clin_ids),
                         if (fl$dup_rows > 0) sprintf(", %s duplicate rows collapsed", wes_fmt(fl$dup_rows)) else "",
                         wes_fmt(fl$maf_samples), wes_fmt(fl$matched), wes_fmt(fl$unmatched),
                         if (fl$unmatched_wt) "counted as WT" else "left out",
                         wes_fmt(fl$dropped_na))
      flow_zh <- sprintf("临床表 %s 行（%s 个编号%s）；MAF %s 个样本；匹配 %s 个；无 MAF 记录 %s 个%s；因时间或事件缺失剔除 %s 个。",
                         wes_fmt(fl$clin_rows), wes_fmt(fl$clin_ids),
                         if (fl$dup_rows > 0) sprintf("，合并重复行 %s 行", wes_fmt(fl$dup_rows)) else "",
                         wes_fmt(fl$maf_samples), wes_fmt(fl$matched), wes_fmt(fl$unmatched),
                         if (fl$unmatched_wt) "，计为野生型" else "，未纳入",
                         wes_fmt(fl$dropped_na))
      sig <- !is.na(s$p) && s$p < 0.05
      verdict_en <- if (sig) "the curves separate (log-rank p < 0.05); check the risk table for how many patients carry the late part"
                    else "no significant separation; with small groups that often means underpowered rather than equal"
      verdict_zh <- if (sig) "曲线分离（log-rank p < 0.05）；请看风险人数表中曲线后段还剩多少患者"
                    else "未见显著分离；组小时这往往意味着效能不足，而非两组真的相同"
      warn <- res$warn
      warn_txt <- if (length(warn)) paste0(" <b>Caution:</b> ", paste(warn, collapse = "; "), ".") else ""
      insight_bar(
        sprintf("<b>%s</b>, %s: %s mutant vs %s WT patients (%s events); median %s months (WT / mutant, NR = not reached); log-rank p = %s; Cox HR %s — %s. %s%s",
                res$label, res$endpoint, wes_fmt(fl$n_mut), wes_fmt(fl$n_wt), wes_fmt(fl$events),
                s$med, if (is.na(s$p)) "-" else format(signif(s$p, 3)), s$hr, verdict_en,
                flow_en, warn_txt),
        sprintf("<b>%s</b>，%s：突变型 %s 人 vs 野生型 %s 人（%s 个事件）；中位生存 %s 个月（野生型/突变型，NR＝未达到）；log-rank p = %s；Cox HR %s——%s。%s%s",
                res$label, res$endpoint, wes_fmt(fl$n_mut), wes_fmt(fl$n_wt), wes_fmt(fl$events),
                s$med, if (is.na(s$p)) "-" else format(signif(s$p, 3)), s$hr, verdict_zh,
                flow_zh, if (length(warn)) paste0("<b>注意：</b>", paste(warn, collapse = "；"), "。")
                         else ""))
    })

    km_gg <- function() {
      shiny::req(res$fit)
      note <- if (length(res$warn)) paste("Caution:", paste(res$warn, collapse = "; ")) else NULL
      km_plot(res$fit, res$lr,
              title = paste0(res$endpoint %||% "Overall survival", " — ",
                             res$label %||% "mutation status"),
              hr = res$hr, note = note)
    }
    render_step_plot(output, input, "km", km_gg, name = "wes_kaplan_meier",
                     width = 8, height = 7)

    view <- shiny::reactive({
      d <- res$df
      shiny::req(d)
      out <- data.frame(patient = d$.id, time_months = round(d$.time, 2), event = d$.event,
                        group = as.character(d$.group), in_MAF = d$.in_maf,
                        stringsAsFactors = FALSE)
      out
    })
    tb <- wes_table(ns, "tbl", view, "wes_survival_cohort",
                    has_maf = function() !is.null(rv$maf),
                    ready = function() !is.null(res$df))
    output$tbl_slot <- tb$slot
    output$tbl <- tb$table
    output$tbl_dl <- tb$download
  })
}
