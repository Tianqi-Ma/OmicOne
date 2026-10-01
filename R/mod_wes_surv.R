#' WES module 11: Mutation vs survival
#'
#' Does carrying a mutation in a gene (or any gene of a set) change outcome? This
#' step reuses the shared survival layer (`fct_survival.R`) rather than a
#' WES-specific one, so the curves, the log-rank test and the Cox screen are
#' identical to the single-cell side — and a cohort loaded once in the
#' *Clinical & survival* step is picked up here automatically.
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
            genes) and compare their survival with a Kaplan-Meier curve and the
            log-rank test.",
      zh = "按某个基因（或一组基因）把患者分为突变型与野生型，用 Kaplan-Meier 曲线和 log-rank 检验比较两者的生存差异。"),
    why  = list(
      en = "Frequency tells you a gene is mutated; survival tells you whether that
            matters to the patient. This is the step that turns a mutation
            landscape into a prognostic claim.",
      zh = "频率只能说明某个基因发生了突变；生存分析才能说明这对患者是否有意义。正是这一步把突变全景变成预后结论。"),
    how  = list(
      en = "Pick one gene to start. Selecting several treats them as a set —
            <i>mutant</i> means a mutation in <b>any</b> of them, which is the
            right framing for a pathway. Survival data comes from the clinical
            table attached to the MAF, or from a cohort you already loaded in the
            single-cell <b>Clinical &amp; survival</b> step.",
      zh = "先从单个基因开始。选择多个基因时会作为一个基因集处理——<i>突变型</i>指其中<b>任意一个</b>发生突变，这正适合用于通路层面的分析。生存数据来自随 MAF 附带的临床表，或你已在单细胞<b>临床与生存</b>步骤中加载的队列。"),
    read = list(
      en = "Each curve is the fraction of patients still event-free over time;
            every step down is an event. The shaded band is the 95% CI — where
            two bands overlap heavily, the visual gap is noise. The p-value is
            the log-rank test over the whole curve, and the risk table below
            shows how many patients each arm still has late on: a separation
            resting on five patients is fragile, whatever the p-value says.",
      zh = "每条曲线是随时间推移仍未发生事件的患者比例；每一次下降就是一个事件。阴影带为 95% 置信区间——若两条带大面积重叠，视觉上的差距只是噪声。p 值来自对整条曲线的 log-rank 检验；图下方的风险人数表显示各组在晚期时间点的剩余人数：如果曲线分离只靠五名患者支撑，无论 p 值多小都不牢靠。"),
    example = list(
      en = "In TCGA LAML, <code>TP53</code> mutants have clearly worse overall
               survival; <code>DNMT3A</code> is the classic borderline case.",
      zh = "在 TCGA LAML 中，<code>TP53</code> 突变型的总生存明显更差；<code>DNMT3A</code> 则是经典的临界案例。")
  )
  controls <- shiny::tagList(
    shiny::uiOutput(ns("gene_ui")),
    shiny::hr(),
    shiny::uiOutput(ns("source_ui")),
    shiny::uiOutput(ns("mapping_ui")),
    run_button(ns("run"), "Run survival analysis", "运行生存分析")
  )
  step_container(
    title     = list(en = "Mutation vs survival", zh = "突变与预后"),
    subtitle  = list(en = "Kaplan-Meier curves for mutant versus wild-type.",
                     zh = "突变型与野生型的 Kaplan-Meier 生存曲线。"),
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
                                       caption = list(en = "Curves: fraction event-free over time; band = 95% CI; table below = patients still at risk.",
                                                      zh = "曲线：随时间的无事件比例；阴影带＝95% CI；下表＝各时点风险人数。"))),
      bslib::nav_panel(i18n("Cohort", "队列表"), shiny::uiOutput(ns("tbl_slot")))
    ))
  )
}

#' @rdname mod_wes_surv
#' @keywords internal
mod_wes_surv_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    ns  <- session$ns
    res <- shiny::reactiveValues(df = NULL, fit = NULL, lr = NULL, label = NULL,
                                 matched = NA_integer_)

    output$gene_ui <- shiny::renderUI({
      shiny::req(rv$maf)
      g <- wes_genes(rv$maf, n = 300)
      shiny::tagList(
        label_with_help("Gene(s)",
                        "One gene, or several treated as a set (mutant = a mutation in any of them).",
                        label_zh = "基因",
                        tip_zh = "可选单个基因，也可选多个作为基因集（突变型 = 其中任一发生突变）。"),
        shiny::selectizeInput(ns("genes"), NULL, choices = g, multiple = TRUE,
                              selected = if (length(g)) g[1] else NULL)
      )
    })

    # Survival data can come from the MAF's own clinical table, or from the
    # cohort the shared Clinical & survival step already normalised.
    has_shared <- shiny::reactive(!is.null(rv$clinical) && nrow(rv$clinical) > 0)

    output$source_ui <- shiny::renderUI({
      shiny::req(rv$maf)
      choices <- c("Clinical table in the MAF" = "maf")
      if (has_shared()) {
        choices <- c("Cohort loaded in Clinical & survival" = "shared", choices)
      }
      shiny::tagList(
        label_with_help("Survival data from",
                        "The shared cohort is whatever you loaded in the Clinical & survival step; it is already normalised.",
                        label_zh = "生存数据来源",
                        tip_zh = "共享队列即你在「临床与生存」步骤中加载的数据，已完成标准化。"),
        shiny::radioButtons(ns("src"), NULL, choices,
                            selected = if (has_shared()) "shared" else "maf")
      )
    })

    output$mapping_ui <- shiny::renderUI({
      shiny::req(rv$maf)
      if (identical(input$src, "shared")) {
        return(shiny::div(class = "omicone-status-empty",
                          i18n("Using the cohort from the Clinical &amp; survival step.",
                               "正在使用「临床与生存」步骤中的队列。")))
      }
      cols <- wes_clinical_cols(rv$maf)
      if (!length(cols)) {
        return(shiny::div(class = "omicone-status-empty",
                          i18n("This MAF carries no clinical columns. Attach a clinical table on the Import step, or load a cohort in the Clinical &amp; survival step.",
                               "该 MAF 未附带临床列。请在导入步骤附上临床表，或在「临床与生存」步骤加载队列。")))
      }
      pick <- function(cands, default = cols[1]) {
        hit <- cols[tolower(cols) %in% cands]
        if (length(hit)) hit[1] else default
      }
      shiny::tagList(
        label_with_help("Follow-up time", "Column holding time to event or last contact.",
                        label_zh = "随访时间", tip_zh = "存放到终点事件或最后随访时间的列。"),
        shiny::selectInput(ns("time_col"), NULL, cols,
                           selected = pick(c("days_to_last_followup", "os_months",
                                             "overall_survival_time", "time", "futime"))),
        shiny::selectInput(ns("time_unit"), i18n("Time unit", "时间单位"),
                           c("Days" = "days", "Months" = "months", "Years" = "years"),
                           selected = "days"),
        label_with_help("Outcome", "1 / Dead / TRUE = the event happened.",
                        label_zh = "终点事件", tip_zh = "1 / Dead / TRUE 表示事件发生。"),
        shiny::selectInput(ns("event_col"), NULL, cols,
                           selected = pick(c("overall_survival_status", "os_status",
                                             "vital_status", "status", "fustat")))
      )
    })

    shiny::observeEvent(input$run, {
      shiny::req(rv$maf, input$genes)
      if (!require_pkgs("maftools", "Mutation vs survival")) return(NULL)

      # 1. per-sample mutation status
      st <- tryCatch(wes_mutation_status(rv$maf, input$genes),
                     error = function(e) {
                       shiny::showNotification(conditionMessage(e), type = "error",
                                               duration = 12); NULL })
      if (is.null(st)) return(NULL)

      # 2. the cohort
      if (identical(input$src, "shared")) {
        clin <- rv$clinical
        if (is.null(clin) || !nrow(clin)) {
          shiny::showNotification("No shared cohort loaded yet.", type = "error",
                                  duration = 10)
          return(NULL)
        }
      } else {
        shiny::req(input$time_col, input$event_col)
        cd <- as.data.frame(maftools::getClinicalData(rv$maf))
        clin <- tryCatch(
          normalise_clinical(cd, "Tumor_Sample_Barcode", input$time_col,
                             input$event_col, time_unit = input$time_unit),
          error = function(e) {
            shiny::showNotification(conditionMessage(e), type = "error", duration = 12)
            NULL
          })
        if (is.null(clin) || !nrow(clin)) {
          shiny::showNotification(
            "No usable follow-up rows: check the time and outcome columns.",
            type = "error", duration = 12)
          return(NULL)
        }
      }

      # 3. join
      clin$.group <- st$status[match(clin$.id, st$.id)]
      matched <- sum(!is.na(clin$.group))
      if (matched < 3) {
        shiny::showNotification(
          paste("Only", matched, "sample id(s) matched between the survival",
                "table and the MAF. Check that both use the same barcodes."),
          type = "error", duration = 15)
        return(NULL)
      }
      clin <- clin[!is.na(clin$.group), , drop = FALSE]
      clin$.group <- factor(clin$.group, levels = c("WT", "Mutant"))
      if (length(unique(clin$.group[!is.na(clin$.group)])) < 2) {
        shiny::showNotification(
          "Every matched patient falls in one group; nothing to compare.",
          type = "error", duration = 12)
        return(NULL)
      }

      fit <- tryCatch(km_fit(clin), error = function(e) {
        shiny::showNotification(paste("Survival fit failed:", conditionMessage(e)),
                                type = "error", duration = 12); NULL })
      if (is.null(fit)) return(NULL)

      res$df      <- clin
      res$fit     <- fit
      res$lr      <- logrank_test(clin)
      res$label   <- paste(input$genes, collapse = " / ")
      res$matched <- matched
      rv$clinical <- clin        # keep the joined cohort for the other pipelines
      mark_done(rv, "wes_surv")
      log_step(log_rv, "WES mutation vs survival",
               params = list(genes = paste(input$genes, collapse = ", "),
                             source = input$src, matched = matched),
               code = c(
                 sprintf('mut <- maftools::subsetMaf(maf, genes = c("%s"))',
                         paste(input$genes, collapse = '", "')),
                 'clin$group <- ifelse(clin$id %in% mut_samples, "Mutant", "WT")',
                 'survival::survdiff(survival::Surv(time, event) ~ group, data = clin)'))
      shiny::showNotification("Survival analysis complete.", type = "message")
    })

    output$summary <- shiny::renderUI({
      if (is.null(rv$maf)) return(wes_no_maf())
      d <- res$df
      if (is.null(d)) {
        return(wes_prompt("Pick a gene and click <b>Run survival analysis</b>.",
                          "选择基因后点击<b>运行生存分析</b>。"))
      }
      med <- km_medians(res$fit)
      p <- res$lr$p
      n_mut <- sum(d$.group == "Mutant", na.rm = TRUE)
      bslib::layout_columns(
        col_widths = c(3, 3, 3, 3),
        stat_tile(i18n("Patients", "患者数"), format(nrow(d), big.mark = ",")),
        stat_tile(i18n("Mutant", "突变型"), format(n_mut, big.mark = ",")),
        stat_tile(i18n("Median OS (mo)", "中位生存(月)"),
                  if (nrow(med)) paste(round(med$median, 1), collapse = " / ") else "-"),
        stat_tile(i18n("Log-rank p", "Log-rank p"),
                  if (is.null(p)) "-" else signif(p, 3))
      )
    })

    output$insight <- shiny::renderUI({
      d <- res$df
      if (is.null(d)) return(NULL)
      p <- res$lr$p
      med <- tryCatch(km_medians(res$fit), error = function(e) NULL)
      med_txt <- if (!is.null(med) && nrow(med)) paste(round(med$median, 1), collapse = " / ") else "?"
      n_mut <- sum(d$.group == "Mutant", na.rm = TRUE)
      n_wt  <- sum(d$.group == "WT", na.rm = TRUE)
      sig_p <- !is.null(p) && is.finite(p) && p < 0.05
      verdict_en <- if (sig_p)
        "the separation is statistically significant — still, check the risk table for how many patients support the late part of the curves"
      else "the separation is not statistically significant; with small groups that often means underpowered rather than equal"
      verdict_zh <- if (sig_p)
        "分离具有统计学显著性——但仍请看风险人数表：曲线后段还剩多少患者支撑"
      else "分离不具统计学显著性；组小时这往往意味着效能不足，而非两组真的相同"
      insight_bar(
        sprintf("<b>%s</b>: %s mutant vs %s wild-type patients; median OS %s months (WT / mutant); log-rank p = %s — %s.",
                res$label, format(n_mut, big.mark = ","), format(n_wt, big.mark = ","),
                med_txt, if (is.null(p)) "?" else signif(p, 3), verdict_en),
        sprintf("<b>%s</b>：突变型 %s 人 vs 野生型 %s 人；中位生存 %s 个月（野生型/突变型）；log-rank p = %s——%s。",
                res$label, format(n_mut, big.mark = ","), format(n_wt, big.mark = ","),
                med_txt, if (is.null(p)) "?" else signif(p, 3), verdict_zh))
    })

    km_gg <- function() {
      shiny::req(res$fit)
      km_plot(res$fit, res$lr,
              title = paste0("Overall survival — ", res$label %||% "mutation status"))
    }
    output$km <- render_scop_plot(km_gg)
    register_figure_download(output, input, "km", function() print(km_gg()),
                             "wes_kaplan_meier", width = 8, height = 6)

    output$tbl_slot <- shiny::renderUI({
      if (is.null(res$df)) return(wes_no_maf())
      if (has_pkg("DT")) DT::dataTableOutput(ns("tbl"))
      else shiny::verbatimTextOutput(ns("tbl_txt"))
    })
    view <- shiny::reactive({
      d <- res$df; shiny::req(d)
      out <- d[, intersect(c(".id", ".time", ".event", ".group"), names(d)), drop = FALSE]
      names(out) <- sub("^\\.", "", names(out))
      out$time <- round(out$time, 1)
      out
    })
    if (has_pkg("DT")) {
      output$tbl <- DT::renderDataTable(
        DT::datatable(view(), rownames = FALSE, filter = "top",
                      options = list(pageLength = 15, scrollX = TRUE)))
    } else {
      output$tbl_txt <- shiny::renderPrint(utils::head(view(), 20))
    }
  })
}
