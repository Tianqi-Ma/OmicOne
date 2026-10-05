# Regression guard for the gap that let a stale function name reach main.
#
# Since the omics landing page moved the pipeline into `output$main_body`'s
# renderUI, `app_ui()` no longer contains any module -- so the shell smoke test
# stopped exercising module UI *or* module outputs. A call to a function that
# does not exist (e.g. a half-applied scStudio -> OmicOne rename) therefore
# passed every test and only failed in the browser.
#
# These tests close that hole two ways:
#   1. build every step's UI for every omics;
#   2. force every module's outputs to evaluate, with and without a working
#      object, asserting nothing raises an unexpected error.
# Outputs are allowed to be *empty* (req()/validate() with no data, or a missing
# Suggests package) -- they are not allowed to raise anything else.

# Conditions that are a *correct* degraded path rather than a defect:
#   - req()/validate(): nothing to show until the user runs the step;
#   - "hasn't been defined yet": an output that is only wired when an optional
#     package is installed (DT, plotly) -- its UI slot is conditional too;
#   - a missing Suggests package, reported by the step's require_pkgs() gate;
#   - the counts accessor refusing the non-Seurat stand-in from helper-fake.R.
expected_blank <- function(e) {
  inherits(e, "shiny.silent.error") ||
    inherits(e, "validation") ||
    grepl(paste("hasn't been defined yet",
                "Seurat|SeuratObject|scop|there is no package",
                "Could not access the counts matrix",
                sep = "|"),
          conditionMessage(e))
}

# Output ids a server function assigns to, read from its source.
#
# `names(output)` inside testServer does not list them (a shinyoutput proxy only
# exposes its internals), and the MockShinySession registry is private, so read
# the `output$<id> <-` assignments off the function instead. That stays correct
# across shiny versions and is enough here: no module builds output ids
# dynamically.
module_output_ids <- function(server_fn) {
  src <- paste(deparse(body(server_fn)), collapse = "\n")
  ids <- regmatches(src, gregexpr("output\\$[A-Za-z0-9_.]+", src))[[1]]
  ids <- sub("^output\\$", "", ids)
  # plots registered through render_step_plot(output, input, "<id>", ...)
  # assign output[[id]] (and its "<id>_dl" download) inside the helper
  reg <- regmatches(src, gregexpr('render_step_plot\\(output, input, "[^"]+"', src))[[1]]
  reg <- sub('.*"([^"]+)"$', "\\1", reg)
  unique(c(ids, reg))
}

# Evaluate every output of a module server, collecting only unexpected errors.
force_outputs <- function(server_fn, id, rv, log_rv) {
  ids <- module_output_ids(server_fn)
  bad <- character(0)
  # require_pkgs() warns instead of notifying outside a running app; that is the
  # expected path here, so it should not clutter the test report.
  suppressWarnings(shiny::testServer(server_fn, args = list(rv = rv, log_rv = log_rv), {
    session$flushReact()
    for (nm in ids) {
      err <- tryCatch({ output[[nm]]; NULL }, error = function(e) e)
      if (!is.null(err) && !expected_blank(err)) {
        bad <<- c(bad, sprintf("%s$%s: %s", id, nm, conditionMessage(err)))
      }
    }
  }))
  bad
}

# The single-cell modules, paired with their step key (= module namespace).
sc_module_servers <- function() {
  list(
    import     = mod_import_server,
    qc         = mod_qc_server,
    doublet    = mod_doublet_server,
    normalize  = mod_normalize_server,
    reduce     = mod_reduce_server,
    integrate  = mod_integrate_server,
    cluster    = mod_cluster_server,
    embed      = mod_embed_server,
    markers    = mod_markers_server,
    annotate   = mod_annotate_server,
    enrichment = mod_enrichment_server,
    trajectory = mod_trajectory_server,
    velocity   = mod_velocity_server,
    dynamic    = mod_dynamic_server,
    cellcycle  = mod_cellcycle_signatures_server,
    cellcomm   = mod_cellcomm_server,
    malignancy = mod_malignancy_server,
    clinical   = mod_clinical_server,
    viz        = mod_viz_server,
    report     = mod_report_server,
    export     = mod_export_server
  )
}

# The WES modules, paired with their step key.
wes_module_servers <- function() {
  list(
    wes_import  = mod_wes_import_server,
    wes_summary = mod_wes_summary_server,
    wes_onco    = mod_wes_onco_server,
    wes_titv    = mod_wes_titv_server,
    wes_tmb     = mod_wes_tmb_server,
    wes_lolli   = mod_wes_lolli_server,
    wes_driver  = mod_wes_driver_server,
    wes_sig     = mod_wes_sig_server,
    wes_clin    = mod_wes_clin_server,
    wes_compare = mod_wes_compare_server,
    wes_surv    = mod_wes_surv_server,
    wes_hetero  = mod_wes_hetero_server
  )
}

test_that("every step's UI builds, for every omics", {
  for (omics in c("sc", "wes", "bulk", "spatial", "integration")) {
    for (s in steps_for(omics)) {
      html <- as.character(do.call(s$ui, list(s$v)))
      expect_true(nchar(html) > 200,
                  info = paste(omics, s$v, "produced a suspiciously small UI"))
    }
  }
})

test_that("the single-cell registry and the wired servers agree", {
  keys <- vapply(steps_sc(), function(s) s$v, character(1))
  expect_setequal(keys, names(sc_module_servers()))
})

test_that("the WES registry and the wired servers agree", {
  keys <- vapply(steps_wes(), function(s) s$v, character(1))
  expect_setequal(keys, names(wes_module_servers()))
  # and none of them is still the placeholder
  for (s in steps_wes()) {
    expect_false(identical(s$ui, mod_placeholder_ui), info = s$v)
  }
})

test_that("only the implemented omics are advertised as available", {
  ready <- vapply(omics_catalogue(), function(o) isTRUE(o$ready), logical(1))
  keys  <- vapply(omics_catalogue(), function(o) o$v, character(1))
  expect_setequal(keys[ready], c("sc", "wes"))
  # every "planned" omics really is still all placeholders
  for (om in keys[!ready]) {
    for (s in steps_for(om)) {
      expect_true(identical(s$ui, mod_placeholder_ui), info = paste(om, s$v))
    }
  }
})

test_that("every WES module's outputs evaluate with no MAF loaded", {
  bad <- character(0)
  for (id in names(wes_module_servers())) {
    rv <- shiny::reactiveValues(omics = "wes", maf = NULL, obj = NULL,
                                status = list(), clinical = NULL)
    bad <- c(bad, force_outputs(wes_module_servers()[[id]], id, rv,
                                shiny::reactiveVal(list())))
  }
  expect_equal(bad, character(0))
})

test_that("every module's outputs evaluate with no data loaded", {
  bad <- character(0)
  for (id in names(sc_module_servers())) {
    rv <- shiny::reactiveValues(obj = NULL, source = NULL, status = list(),
                                clinical = NULL)
    bad <- c(bad, force_outputs(sc_module_servers()[[id]], id, rv,
                                shiny::reactiveVal(list())))
  }
  expect_equal(bad, character(0))
})

test_that("every module's outputs evaluate against a working object", {
  # A minimal stand-in for a Seurat object: modules must not blow up on it, they
  # should degrade to "can't do that yet" the same way they do on real data
  # without the Suggests packages installed.
  obj <- fake_obj()
  bad <- character(0)
  for (id in names(sc_module_servers())) {
    rv <- shiny::reactiveValues(obj = obj, source = "test", status = list(),
                                clinical = NULL)
    bad <- c(bad, force_outputs(sc_module_servers()[[id]], id, rv,
                                shiny::reactiveVal(list())))
  }
  expect_equal(bad, character(0))
})

test_that("the shared export handlers evaluate without data", {
  ids <- module_output_ids(register_exports)
  expect_true(length(ids) >= 4)   # object / metadata / counts / embeddings
  bad <- character(0)
  shiny::testServer(function(input, output, session) {
    rv <- shiny::reactiveValues(obj = NULL)
    register_exports(input, output, session, rv)
  }, {
    session$flushReact()
    for (nm in ids) {
      err <- tryCatch({ output[[nm]]; NULL }, error = function(e) e)
      if (!is.null(err) && !expected_blank(err)) {
        bad <<- c(bad, sprintf("%s: %s", nm, conditionMessage(err)))
      }
    }
  })
  expect_equal(bad, character(0))
})

test_that("no package function calls a name that does not exist", {
  # The direct guard against a half-applied rename: walk every function in the
  # package and check each unqualified call resolves. This is what would have
  # caught `scstudio_theme()` surviving in fct_overview.R.
  skip_if_not_installed("codetools")
  env <- environment(overview_plots)          # the package/test namespace
  missing <- character(0)
  for (nm in ls(env, all.names = TRUE)) {
    f <- get(nm, envir = env)
    if (!is.function(f)) next
    globs <- tryCatch(codetools::findGlobals(f, merge = FALSE)$functions,
                      error = function(e) character(0))
    for (g in globs) {
      if (grepl("^[.]|::", g)) next
      if (!exists(g, envir = env)) {
        missing <- c(missing, sprintf("%s() calls undefined %s()", nm, g))
      }
    }
  }
  expect_equal(unique(missing), character(0))
})
