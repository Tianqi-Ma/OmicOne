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
  rv <- shiny::reactiveValues(
    omics = NULL, status = list(),
    obj = NULL, source = NULL,          # single-cell
    maf = NULL, maf_source = NULL,      # WES
    clinical = NULL                     # shared across omics
  )
  log_rv <- shiny::reactiveVal(list())

  # --- Shared "export current data" (available on every step) ----------------
  register_exports(input, output, session, rv)

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
        state <- if (identical(s$v, current)) "current"
                 else if (isTRUE(status[[s$v]])) "done" else "todo"
        shiny::tags$a(
          class = paste("omicone-navitem", state),
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
  mod_wes_hetero_server("wes_hetero", rv, log_rv)

  # --- Placeholder servers for the still-planned omics -----------------------
  ph_ids <- unlist(lapply(c("bulk", "spatial", "integration"),
                          function(o) vapply(steps_for(o), function(s) s$v, character(1))))
  for (pid in ph_ids) mod_placeholder_server(pid, rv, log_rv)
}

#' NULL-coalescing helper
#' @keywords internal
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Mark a pipeline step completed (colours the navigator)
#' @param rv Shared hub. @param step Step key.
#' @keywords internal
mark_done <- function(rv, step) {
  st <- rv$status; st[[step]] <- TRUE; rv$status <- st; invisible(TRUE)
}

#' Small labelled status line for the sidebar
#' @keywords internal
stat_line <- function(label, value) {
  shiny::div(class = "omicone-statline",
             shiny::span(class = "omicone-statlabel", label),
             shiny::span(class = "omicone-statvalue", value))
}
