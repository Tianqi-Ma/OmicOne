#' Pipeline state: step dependencies, staleness, data epochs, checkpoints
#'
#' Modules share one working object per omics (`rv$obj`, `rv$maf`) and write
#' their result back into it. That keeps each module small, but on its own it
#' has three failure modes, which this file closes:
#'
#' 1. **Stale downstream results.** Re-running Normalize after Cluster leaves
#'    the clusters computed on the old normalisation. Every registry entry
#'    declares the steps it consumes (`deps`); when a step completes, each
#'    finished step downstream of it becomes `"stale"`, and the navigator and
#'    a banner above the workspace say so.
#' 2. **A new dataset inheriting the old one's results.** A successful import
#'    starts a new *epoch* for that omics: its statuses, log entries and
#'    checkpoints are cleared, and every module store created with
#'    [step_results()] / [step_result()] resets to its initial value.
#' 3. **Non-idempotent re-runs.** Filtering steps (QC, doublet removal) read
#'    their input through [step_input()], which pins the object they first
#'    received. Running QC twice therefore filters the imported cells with the
#'    new thresholds, instead of filtering the already-filtered cells again.
#'
#' `rv$status[[step]]` is `NULL` (not run), `TRUE` (done) or `"stale"`.
#'
#' @name fct_state
#' @keywords internal
NULL

#' Which omics registry a step key belongs to
#' @param step Step key, e.g. "qc" or "wes_tmb".
#' @return The omics key ("sc", "wes", ...) or NULL.
#' @keywords internal
step_omics <- function(step) {
  regs <- all_step_registries()
  for (om in names(regs)) {
    if (step %in% vapply(regs[[om]], function(s) s$v, character(1))) return(om)
  }
  NULL
}

#' Every step that (transitively) consumes the output of `step`
#'
#' Walks the `deps` declared in the step's own registry. A step never lists
#' itself, and the walk is cycle-safe.
#' @param step Step key.
#' @return Character vector of dependent step keys (possibly empty).
#' @keywords internal
step_dependents <- function(step) {
  om <- step_omics(step)
  if (is.null(om)) return(character(0))
  steps <- all_step_registries()[[om]]
  deps <- stats::setNames(lapply(steps, function(s) s$deps %||% character(0)),
                          vapply(steps, function(s) s$v, character(1)))
  out <- character(0)
  frontier <- step
  while (length(frontier)) {
    nxt <- names(deps)[vapply(deps, function(d) any(d %in% frontier), logical(1))]
    nxt <- setdiff(nxt, c(out, step))
    out <- c(out, nxt)
    frontier <- nxt
  }
  out
}

#' Status of one step: "done", "stale" or "todo"
#' @param rv Shared hub. @param step Step key.
#' @keywords internal
step_state <- function(rv, step) {
  s <- rv$status[[step]]
  if (isTRUE(s)) "done" else if (identical(s, "stale")) "stale" else "todo"
}

#' Mark a pipeline step completed, and its finished dependents stale
#'
#' Also drops the checkpoints of every dependent: their input has changed, so
#' the next run must take the fresh object rather than a pinned old one.
#' @param rv Shared hub. @param step Step key.
#' @keywords internal
mark_done <- function(rv, step) {
  st <- rv$status
  if (!is.list(st)) st <- as.list(st)
  st[[step]] <- TRUE
  for (d in step_dependents(step)) {
    if (isTRUE(st[[d]])) st[[d]] <- "stale"
    ckpt_drop(rv, d)
  }
  rv$status <- st
  invisible(TRUE)
}

#' Start a new data epoch for one omics (call right after a successful import)
#'
#' Clears that omics' step statuses, checkpoints and log entries, then bumps
#' `rv$epoch_<omics>`, which resets every module store built with
#' [step_results()] / [step_result()]. The other omics are left untouched.
#' @param rv Shared hub. @param omics Omics key ("sc", "wes").
#' @param log_rv The reproducibility log (`reactiveVal`), or NULL.
#' @keywords internal
start_epoch <- function(rv, omics, log_rv = NULL) {
  keys <- vapply(all_step_registries()[[omics]] %||% list(),
                 function(s) s$v, character(1))
  st <- rv$status
  if (!is.list(st)) st <- as.list(st)
  st[intersect(names(st), keys)] <- NULL
  rv$status <- st
  for (k in keys) ckpt_drop(rv, k)
  if (!is.null(log_rv)) {
    log_rv(Filter(function(e) !identical(e$omics %||% "sc", omics),
                  shiny::isolate(log_rv())))
  }
  nm <- paste0("epoch_", omics)
  rv[[nm]] <- (shiny::isolate(rv[[nm]]) %||% 0L) + 1L
  invisible(TRUE)
}

#' Module result store that resets when its omics starts a new epoch
#'
#' Drop-in replacement for `shiny::reactiveValues(...)` inside a module
#' server: same fields, same initial values, restored whenever a new dataset
#' is imported for `omics`. Do not use it in the import module itself (its
#' own epoch bump would wipe the result it just stored).
#' @param rv Shared hub. @param omics Omics key the module belongs to.
#' @param ... Named initial values.
#' @keywords internal
step_results <- function(rv, omics, ...) {
  init <- list(...)
  res <- do.call(shiny::reactiveValues, init)
  shiny::observeEvent(rv[[paste0("epoch_", omics)]], {
    for (nm in names(init)) res[[nm]] <- init[[nm]]
  }, ignoreInit = TRUE)
  res
}

#' Single-value version of [step_results()] (a resetting `reactiveVal`)
#' @param rv Shared hub. @param omics Omics key. @param value Initial value.
#' @keywords internal
step_result <- function(rv, omics, value = NULL) {
  val <- shiny::reactiveVal(value)
  shiny::observeEvent(rv[[paste0("epoch_", omics)]], val(value),
                      ignoreInit = TRUE)
  val
}

#' The object a filtering step should start from
#'
#' The first run pins `current` as the step's input; later runs get that same
#' object back, so re-running with new thresholds never compounds. The pin is
#' dropped by [mark_done()] on any upstream step and by [start_epoch()].
#' @param rv Shared hub (`rv$ckpt` is a plain environment).
#' @param step Step key. @param current The object the step would otherwise use.
#' @keywords internal
step_input <- function(rv, step, current) {
  ck <- shiny::isolate(rv$ckpt)
  if (!is.environment(ck)) return(current)
  if (!exists(step, envir = ck, inherits = FALSE)) assign(step, current, envir = ck)
  get(step, envir = ck, inherits = FALSE)
}

#' Drop a step's pinned input (see [step_input()])
#' @keywords internal
ckpt_drop <- function(rv, step) {
  ck <- shiny::isolate(rv$ckpt)
  if (is.environment(ck) && exists(step, envir = ck, inherits = FALSE)) {
    rm(list = step, envir = ck)
  }
  invisible(NULL)
}

#' Banner shown above a step whose upstream was re-run after it
#' @param rv Shared hub. @param current Current step key.
#' @keywords internal
stale_banner <- function(rv, current) {
  if (is.null(current) || !identical(step_state(rv, current), "stale")) return(NULL)
  shiny::div(
    class = "omicone-stalebar", role = "status",
    shiny::icon("rotate"),
    i18n("An earlier step was re-run after this one, so the results shown here may no longer match the data. Run this step again to refresh them.",
         "在这一步之后，前面的某个步骤被重新运行过，这里显示的结果可能已与当前数据不符。请重新运行这一步以更新结果。")
  )
}
