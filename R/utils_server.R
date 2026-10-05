#' Reusable server-side helpers for OmicOne
#'
#' @name utils_server
#' @keywords internal
NULL

#' NULL-coalescing helper (base R has it only from 4.4)
#' @keywords internal
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Check that suggested packages are installed, notify if not
#'
#' Heavy compute packages (Seurat, scDblFinder, harmony, SingleR, ...) live in
#' Suggests so the app installs and the UI loads even without them. Each compute
#' step calls this first and refuses to run (with a helpful install hint) if a
#' required package is missing.
#'
#' @param pkgs Character vector of package names required by the step.
#' @param what Character. Short description of the step, used in the message.
#' @return `TRUE` if all present; otherwise `FALSE` (and shows a notification).
#' @keywords internal
require_pkgs <- function(pkgs, what = "this step") {
  missing <- pkgs[!vapply(pkgs, function(p) {
    requireNamespace(p, quietly = TRUE)
  }, logical(1))]
  if (length(missing) == 0) return(TRUE)
  msg <- sprintf(
    "%s needs package(s) not installed: %s. Install them, then retry.",
    what, paste(missing, collapse = ", ")
  )
  if (shiny::isRunning()) {
    shiny::showNotification(
      i18n(msg, sprintf("%s 需要尚未安装的包：%s。请先安装再重试。", what,
                        paste(missing, collapse = ", "))),
      type = "error", duration = 10)
  } else {
    warning(msg, call. = FALSE)
  }
  FALSE
}

#' Is a package available?
#' @param pkg Package name.
#' @keywords internal
has_pkg <- function(pkg) requireNamespace(pkg, quietly = TRUE)

#' Render a data.frame as an interactive table (DT) or a plain fallback
#'
#' Pairs with `DT::dataTableOutput(id)` / `shiny::verbatimTextOutput(id)` chosen
#' at UI-build time by whether DT is installed.
#'
#' Large tables are capped at `max_rows` for the browser (AGENTS.md section 3);
#' offer the full table as a CSV download next to it.
#'
#' @param data_fn A function returning a data.frame (may call req()).
#' @param max_rows Rows sent to the browser.
#' @keywords internal
render_tbl_wrap <- function(data_fn, max_rows = 5000) {
  if (has_pkg("DT")) {
    DT::renderDataTable({
      df <- data_fn()
      shiny::req(df)
      DT::datatable(utils::head(df, max_rows), options = list(pageLength = 15, scrollX = TRUE),
                    rownames = TRUE, class = "compact stripe")
    })
  } else {
    shiny::renderPrint({
      df <- data_fn()
      shiny::req(df)
      utils::head(df, 20)
    })
  }
}

#' Run an expression with a Shiny progress bar and graceful error notify
#'
#' Wraps a (possibly slow) compute call so the UI shows progress and any error
#' surfaces as a notification instead of crashing the session.
#'
#' @param expr Expression to evaluate.
#' @param message Progress message shown to the user.
#' @param session Shiny session (for progress).
#' @return The value of `expr`, or `NULL` on error.
#' @keywords internal
with_progress_notify <- function(expr, message = "Working...", session = shiny::getDefaultReactiveDomain()) {
  # The progress panel is drawn by Shiny's own client code, outside the i18n
  # swap, so its message stays a plain string; errors below are bilingual.
  prog <- shiny::Progress$new(session)
  prog$set(message = message, value = 0.1)
  on.exit(prog$close(), add = TRUE)
  out <- tryCatch({
    prog$set(value = 0.5)
    force(expr)
  }, error = function(e) {
    shiny::showNotification(i18n(paste("Error:", conditionMessage(e)),
                                 paste("出错：", conditionMessage(e))),
                            type = "error", duration = 12)
    NULL
  })
  prog$set(value = 1)
  out
}

#' Record a reproducibility log entry (step + parameters + equivalent R code)
#'
#' Called from inside a module server, the entry is tagged with the module's
#' step key and omics (read off the session namespace), so the report and the
#' exported script can keep each pipeline's entries apart. Re-running a step
#' *replaces* its earlier entry in place rather than appending a second one:
#' the log describes the analysis as it now stands, which is what a script
#' replaying it must reproduce. Entries are cleared by [start_epoch()].
#'
#' @param log_rv A `reactiveVal` holding a list of log entries.
#' @param step Character step name (a module may log several distinct names).
#' @param params Named list of chosen parameters.
#' @param code Character vector of equivalent R code lines.
#' @param key Step key; defaults to the calling module's id.
#' @keywords internal
log_step <- function(log_rv, step, params = list(), code = character(0),
                     key = NULL) {
  if (is.null(key)) {
    dom <- shiny::getDefaultReactiveDomain()
    ns <- if (!is.null(dom)) tryCatch(dom$ns(""), error = function(e) "") else ""
    key <- sub("-$", "", ns)
  }
  entry <- list(
    step = step,
    key = key,
    omics = (if (nzchar(key)) step_omics(key) else NULL) %||% "sc",
    time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    params = params,
    code = code
  )
  cur <- shiny::isolate(log_rv())
  same <- vapply(cur, function(e) identical(e$key, entry$key) &&
                   identical(e$step, entry$step), logical(1))
  if (any(same)) cur[[which(same)[1]]] <- entry else cur[[length(cur) + 1]] <- entry
  log_rv(cur)
  invisible(entry)
}

#' Log entries belonging to one omics pipeline
#' @param entries List of log entries. @param omics Omics key.
#' @keywords internal
log_entries_for <- function(entries, omics = "sc") {
  Filter(function(e) identical(e$omics %||% "sc", omics), entries %||% list())
}

#' Guess the sample/batch column in an object's metadata
#'
#' @param meta A data.frame of cell metadata.
#' @return A best-guess column name, or NULL.
#' @keywords internal
guess_batch_col <- function(meta) {
  if (is.null(meta) || ncol(meta) == 0) return(NULL)
  candidates <- c("sample", "orig.ident", "batch", "donor", "patient",
                  "Sample", "Batch", "sampleID", "sample_id")
  hit <- candidates[candidates %in% colnames(meta)]
  if (length(hit)) return(hit[1])
  NULL
}

#' Human-readable size of the current dataset for memory guardrails
#'
#' @param n_cells Integer number of cells.
#' @return A short advisory string, or "" if within comfortable range.
#' @keywords internal
memory_advice <- function(n_cells) {
  if (is.null(n_cells) || is.na(n_cells)) return("")
  if (n_cells > 5e5) {
    return(paste0("This dataset has ", format(n_cells, big.mark = ","),
                  " cells. Expect high memory use (>32 GB) and slow steps; ",
                  "consider downsampling for exploration."))
  }
  if (n_cells > 1e5) {
    return(paste0(format(n_cells, big.mark = ","),
                  " cells: plan for 16-32 GB RAM; heavy steps may take minutes."))
  }
  ""
}
