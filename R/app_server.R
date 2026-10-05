#' Top-level server: omics routing, per-omics stepper, shared state
#'
#' `rv$omics` selects the active pipeline (NULL = landing page). The left
#' navigator and the main workspace are rendered from that omics' step registry
#' (see steps.R). A single `rv` hub carries the working object, per-step status,
#' the clinical/survival table (shared across omics), and the reproducibility log.
#'
#' @param input,output,session Standard Shiny server arguments.
#' @keywords internal
app_server <- function(input, output, session) {

  # One slot per omics for the working object (`obj` = single-cell Seurat,
  # `maf` = maftools MAF, ...) so switching pipelines cannot leave the sidebar
  # reporting the previous pipeline's data. `clinical` is deliberately shared:
  # a cohort loaded once is reused by every pipeline that needs outcome.
  # `epoch_*`, `ckpt` and the stale entries of `status` are the pipeline-state
  # bookkeeping described in fct_state.R.
  rv <- shiny::reactiveValues(
    omics = NULL, status = list(),
    obj = NULL, source = NULL,          # single-cell
    maf = NULL, maf_source = NULL, wes_sequenced = NULL,   # WES
    wes_clin_raw = NULL,                # clinical table as read with the MAF
    clinical = NULL,                    # shared across omics
    epoch_sc = 0L, epoch_wes = 0L,
    ckpt = new.env(parent = emptyenv())
  )
  log_rv <- shiny::reactiveVal(list())

  # --- Shared "export current data" (available on every step) ----------------
  register_exports(input, output, session, rv,
                   omics_items = list(wes = wes_export_items))
  register_wes_exports(output, rv)

  # --- Omics routing ---------------------------------------------------------
  shiny::observeEvent(input$omics, { rv$omics <- input$omics })
  shiny::observeEvent(input$switch_omics, { rv$omics <- NULL })

  # Main body: landing page, or the selected omics' hidden tabset.
  output$main_body <- shiny::renderUI({
    if (is.null(rv$omics)) return(app_landing())
    steps <- steps_for(rv$omics)
    panels <- lapply(steps, function(s) {
      bslib::nav_panel(title = s$en, value = s$v, do.call(s$ui, list(s$v)))
    })
    do.call(bslib::navset_hidden, c(list(id = "steps"), panels))
  })

  # Left navigator: collapsible phase groups, status-coloured; empty on the
  # landing page. Only the current step's group (and the first) stays expanded.
  output$step_nav <- shiny::renderUI({
    if (is.null(rv$omics)) {
      return(shiny::div(class = "omicone-status-empty",
                        i18n("Pick an analysis to begin.", "选择一个分析开始。")))
    }
    current <- input$steps %||% steps_for(rv$omics)[[1]]$v
    status  <- rv$status
    phases  <- app_phases()
    steps   <- steps_for(rv$omics)
    children <- list()
    for (ph in phase_order_for(rv$omics)) {
      lab <- phases[[ph]]
      ph_steps <- Filter(function(x) identical(x$phase, ph), steps)
      items <- lapply(ph_steps, function(s) {
        st <- step_state(rv, s$v)
        state <- if (identical(s$v, current)) "current" else st
        shiny::tags$a(
          class = paste("omicone-navitem", state, if (st == "stale") "is-stale"),
          title = if (st == "stale") "Needs a re-run: an earlier step changed",
          onclick = sprintf("Shiny.setInputValue('goto','%s',{priority:'event'})", s$v),
          shiny::span(class = "omicone-navdot"),
          shiny::span(class = "omicone-navnum", s$n),
          shiny::span(class = "omicone-navlabel", i18n(s$en, s$zh))
        )
      })
      n_done <- sum(vapply(ph_steps, function(s) isTRUE(status[[s$v]]), logical(1)))
      has_current <- any(vapply(ph_steps, function(s) identical(s$v, current),
                                logical(1)))
      children[[length(children) + 1]] <- shiny::tags$details(
        class = "omicone-phasegroup",
        open = if (has_current || !length(children)) NA else NULL,
        shiny::tags$summary(
          class = "omicone-phase",
          i18n(lab$en, lab$zh),
          shiny::span(class = "omicone-phase-count",
                      sprintf("%d/%d", n_done, length(ph_steps)))
        ),
        items
      )
    }
    shiny::div(class = "omicone-nav", children)
  })

  shiny::observeEvent(input$goto, { bslib::nav_select("steps", input$goto) })

  # --- Progress chip (topbar): done/total + the next unfinished step ---------
  # The "next step" link lives here, in the always-visible top bar, rather than
  # floating over the workspace where it covered the figure-export row.
  output$progress_chip <- shiny::renderUI({
    if (is.null(rv$omics)) return(NULL)
    steps <- steps_for(rv$omics)
    n <- length(steps)
    if (!n) return(NULL)
    keys <- vapply(steps, function(s) s$v, character(1))
    is_done <- vapply(keys, function(k) step_state(rv, k) == "done", logical(1))
    done <- sum(is_done)
    nxt <- NULL
    cur <- input$steps
    if (!is.null(cur) && isTRUE(is_done[cur]) && !all(is_done)) {
      i <- match(cur, keys)
      undone <- which(!is_done)
      j <- undone[undone > i][1]
      if (is.na(j)) j <- undone[1]
      s <- steps[[j]]
      nxt <- shiny::tags$a(
        class = "omicone-nextlink",
        onclick = sprintf("Shiny.setInputValue('goto','%s',{priority:'event'})", s$v),
        shiny::span(class = "omicone-next-label", i18n("Next:", "下一步：")),
        shiny::strong(i18n(s$en, s$zh)), " \u2192")
    }
    shiny::div(class = "omicone-progress",
               title = i18n("Steps completed", "已完成步骤"),
               shiny::span(class = "omicone-progress-num",
                           sprintf("%d/%d", done, n)),
               shiny::div(class = "omicone-progress-bar",
                          shiny::div(class = "omicone-progress-fill",
                                     style = sprintf("width:%.0f%%", 100 * done / n))),
               nxt)
  })

  # --- Banner above the workspace when the current step is stale ------------
  output$stale_banner <- shiny::renderUI({
    if (is.null(rv$omics)) return(NULL)
    stale_banner(rv, input$steps)
  })

  # --- Dataset status (bottom of sidebar), for the active omics --------------
  output$global_status <- shiny::renderUI({
    empty <- function() shiny::div(class = "omicone-status-empty",
                                   i18n("No data loaded.", "尚未加载数据"))
    cohort <- if (!is.null(rv$clinical) && nrow(rv$clinical))
      stat_line(i18n("Cohort", "队列"),
                format(nrow(rv$clinical), big.mark = ",")) else NULL

    if (identical(rv$omics, "wes")) {
      if (is.null(rv$maf)) return(empty())
      ov <- wes_overview(rv$maf)
      return(shiny::tagList(
        stat_line(i18n("Samples", "样本"), format(ov$samples, big.mark = ",")),
        stat_line(i18n("Genes", "基因"), format(ov$genes, big.mark = ",")),
        stat_line(i18n("Variants", "变异"), format(ov$variants, big.mark = ",")),
        cohort
      ))
    }

    obj <- rv$obj
    if (is.null(obj)) return(if (is.null(cohort)) empty() else shiny::tagList(cohort))
    dims <- obj_dims(obj)
    advice <- memory_advice(dims$cells)
    shiny::tagList(
      stat_line(i18n("Cells", "细胞"), format(dims$cells, big.mark = ",")),
      stat_line(i18n("Genes", "基因"), format(dims$genes, big.mark = ",")),
      cohort,
      if (nzchar(advice))
        shiny::div(class = "omicone-warn", shiny::icon("triangle-exclamation"), " ", advice)
    )
  })

  # --- Single-cell module servers (fully implemented) ------------------------
  mod_import_server("import", rv, log_rv, parent = session)
  mod_qc_server("qc", rv, log_rv)
  mod_doublet_server("doublet", rv, log_rv)
  mod_normalize_server("normalize", rv, log_rv)
  mod_reduce_server("reduce", rv, log_rv)
  mod_integrate_server("integrate", rv, log_rv)
  mod_cluster_server("cluster", rv, log_rv)
  mod_embed_server("embed", rv, log_rv)
  mod_markers_server("markers", rv, log_rv)
  mod_annotate_server("annotate", rv, log_rv)
  mod_enrichment_server("enrichment", rv, log_rv)
  mod_pseudobulk_server("pseudobulk", rv, log_rv)
  mod_abundance_server("abundance", rv, log_rv)
  mod_trajectory_server("trajectory", rv, log_rv)
  mod_velocity_server("velocity", rv, log_rv)
  mod_dynamic_server("dynamic", rv, log_rv)
  mod_cellcycle_signatures_server("cellcycle", rv, log_rv)
  mod_cellcomm_server("cellcomm", rv, log_rv)
  mod_malignancy_server("malignancy", rv, log_rv)
  mod_clinical_server("clinical", rv, log_rv)
  mod_viz_server("viz", rv, log_rv)
  mod_report_server("report", rv, log_rv)
  mod_export_server("export", rv, log_rv)

  # --- WES module servers (fully implemented) --------------------------------
  mod_wes_import_server("wes_import", rv, log_rv)
  mod_wes_filter_server("wes_filter", rv, log_rv)
  mod_wes_summary_server("wes_summary", rv, log_rv)
  mod_wes_onco_server("wes_onco", rv, log_rv)
  mod_wes_titv_server("wes_titv", rv, log_rv)
  mod_wes_tmb_server("wes_tmb", rv, log_rv)
  mod_wes_lolli_server("wes_lolli", rv, log_rv)
  mod_wes_driver_server("wes_driver", rv, log_rv)
  mod_wes_sig_server("wes_sig", rv, log_rv)
  mod_wes_clin_server("wes_clin", rv, log_rv)
  mod_wes_compare_server("wes_compare", rv, log_rv)
  mod_wes_surv_server("wes_surv", rv, log_rv)
  mod_wes_tmbclin_server("wes_tmbclin", rv, log_rv)
  mod_wes_hetero_server("wes_hetero", rv, log_rv)
  mod_wes_report_server("wes_report", rv, log_rv)

  # --- Placeholder servers for the still-planned omics -----------------------
  ph_ids <- unlist(lapply(c("bulk", "spatial", "integration"),
                          function(o) vapply(steps_for(o), function(s) s$v, character(1))))
  for (pid in ph_ids) mod_placeholder_server(pid, rv, log_rv)
}
