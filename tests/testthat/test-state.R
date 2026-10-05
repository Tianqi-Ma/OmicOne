# Pipeline-state rules (fct_state.R): step dependencies, staleness, data
# epochs and the pinned inputs of filtering steps. Pure reactive logic -- no
# Seurat or maftools needed.

test_that("every declared dependency is an earlier step of the same registry", {
  for (om in names(all_step_registries())) {
    steps <- all_step_registries()[[om]]
    keys <- vapply(steps, function(s) s$v, character(1))
    for (i in seq_along(steps)) {
      d <- steps[[i]]$deps
      expect_true(is.character(d), info = paste(om, keys[i], "has no deps field"))
      expect_equal(setdiff(d, keys[seq_len(i - 1)]), character(0),
                   info = paste(om, keys[i]))
    }
  }
})

test_that("dependents are transitive and output steps are never stale", {
  expect_true(all(c("doublet", "normalize", "cluster", "markers", "annotate") %in%
                    step_dependents("qc")))
  expect_true("dynamic" %in% step_dependents("normalize"))   # via trajectory
  expect_false("qc" %in% step_dependents("normalize"))       # never upstream
  expect_false(any(c("report", "export") %in% step_dependents("import")))
  expect_equal(step_dependents("viz"), character(0))
})

test_that("re-running a step marks finished dependents stale, not unfinished ones", {
  rv <- shiny::reactiveValues(status = list(), epoch_sc = 0L,
                              ckpt = new.env(parent = emptyenv()))
  shiny::isolate({
    for (k in c("import", "qc", "normalize", "reduce", "cluster")) mark_done(rv, k)
    mark_done(rv, "normalize")
    expect_equal(step_state(rv, "normalize"), "done")
    expect_equal(step_state(rv, "reduce"), "stale")
    expect_equal(step_state(rv, "cluster"), "stale")
    expect_equal(step_state(rv, "qc"), "done")
    expect_equal(step_state(rv, "markers"), "todo")      # never ran: stays todo
    # finishing the stale step makes it done again
    mark_done(rv, "reduce")
    expect_equal(step_state(rv, "reduce"), "done")
  })
})

test_that("filtering steps re-run from their pinned input, dropped upstream", {
  rv <- shiny::reactiveValues(status = list(), ckpt = new.env(parent = emptyenv()))
  shiny::isolate({
    expect_equal(step_input(rv, "qc", "imported"), "imported")
    # second run: the current object is QC's own output, but QC gets its input
    expect_equal(step_input(rv, "qc", "already-filtered"), "imported")
    expect_equal(step_input(rv, "doublet", "after-qc-1"), "after-qc-1")
    mark_done(rv, "qc")                              # QC re-ran: doublet's input changed
    expect_equal(step_input(rv, "doublet", "after-qc-2"), "after-qc-2")
  })
})

test_that("a new epoch clears only its own omics and resets module stores", {
  rv <- shiny::reactiveValues(status = list(), epoch_sc = 0L, epoch_wes = 0L,
                              ckpt = new.env(parent = emptyenv()))
  log_rv <- shiny::reactiveVal(list(
    list(step = "QC", key = "qc", omics = "sc"),
    list(step = "WES import", key = "wes_import", omics = "wes")))
  shiny::isolate({
    mark_done(rv, "qc")
    step_input(rv, "qc", "obj")
    if ("wes" %in% names(all_step_registries())) mark_done(rv, "wes_import")
    start_epoch(rv, "sc", log_rv)
    expect_equal(step_state(rv, "qc"), "todo")
    expect_false(exists("qc", envir = rv$ckpt))
    expect_equal(rv$epoch_sc, 1L)
    expect_equal(vapply(log_rv(), function(e) e$omics, ""), "wes")
    if ("wes" %in% names(all_step_registries())) {
      expect_equal(step_state(rv, "wes_import"), "done")
    }
  })
})

test_that("step_results() stores reset when their omics starts a new epoch", {
  shiny::testServer(function(input, output, session) {
    rv <- shiny::reactiveValues(status = list(), epoch_sc = 0L)
    res <- step_results(rv, "sc", done = FALSE, n = NA)
    one <- step_result(rv, "sc", NULL)
    session$userData$rv <- rv
    session$userData$res <- res
    session$userData$one <- one
  }, {
    rv <- session$userData$rv
    res <- session$userData$res
    one <- session$userData$one
    res$done <- TRUE
    res$n <- 42
    one("table")
    session$flushReact()
    rv$epoch_sc <- 1L
    session$flushReact()
    expect_false(res$done)
    expect_true(is.na(res$n))
    expect_null(one())
  })
})

test_that("log_step replaces a re-run step's entry instead of appending", {
  log_rv <- shiny::reactiveVal(list())
  shiny::isolate({
    log_step(log_rv, "QC", params = list(nmad = 5), key = "qc")
    log_step(log_rv, "Normalization", key = "normalize")
    log_step(log_rv, "QC", params = list(nmad = 3), key = "qc")
    entries <- log_rv()
    expect_length(entries, 2)
    expect_equal(entries[[1]]$step, "QC")              # keeps its position
    expect_equal(entries[[1]]$params$nmad, 3)          # with the new parameters
    expect_equal(entries[[1]]$omics, "sc")
    expect_length(log_entries_for(entries, "wes"), 0)
  })
})
