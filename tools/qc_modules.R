#!/usr/bin/env Rscript
# Module QC: use every module like a user would -- set every input to its UI
# default (also the inputs created by renderUI), press every Run / Apply button,
# then render every output at 600x400 down to 8x8. Reports unexpected errors,
# plot errors (the renderer draws them as a message and a notification) and
# graphics devices left open. Run from the repository root:
#
#   Rscript tools/qc_modules.R [module-regex] [--pbmc3k]
#
# --pbmc3k adds a pass on the bundled pbmc3k (slow: minutes) for the steps the
# tiny demo cannot run (cell cycle, markers on real clusters). Steps whose
# Suggests packages are missing are reported as "not run", which is expected.
# Exit 0 = no issues, 1 = issues printed.
args <- commandArgs(trailingOnly = TRUE)
use_pbmc <- "--pbmc3k" %in% args
args <- setdiff(args, "--pbmc3k")
only <- if (length(args)) args[1] else "."
repo <- "."
if (nzchar(Sys.getenv("OMICONE_QC_LIB"))) .libPaths(c(Sys.getenv("OMICONE_QC_LIB"), .libPaths()))
suppressMessages(pkgload::load_all(repo, quiet = TRUE))
pkg <- read.dcf("DESCRIPTION", "Package")[1, 1]
pkgenv <- asNamespace(pkg)

# helpers copied from tests/testthat/test-outputs.R
expected_blank <- function(e) {
  inherits(e, "shiny.silent.error") || inherits(e, "validation") ||
    grepl("hasn't been defined yet|Seurat|SeuratObject|scop|there is no package|Could not access the counts matrix",
          conditionMessage(e))
}
module_output_ids <- function(server_fn) {
  src <- paste(deparse(body(server_fn)), collapse = "\n")
  ids <- sub("^output\\$", "", regmatches(src, gregexpr("output\\$[A-Za-z0-9_.]+", src))[[1]])
  reg <- regmatches(src, gregexpr('render_step_plot\\(output, input, "[^"]+"', src))[[1]]
  unique(c(ids, sub('.*"([^"]+)"$', "\\1", reg)))
}

# ---- instrumentation ------------------------------------------------------------
SIZE <- c(w = 600, h = 400)
cd_get <- function(x, name) {
  if (grepl("_width$", name)) return(SIZE[["w"]])
  if (grepl("_height$", name)) return(SIZE[["h"]])
  if (identical(name, "pixelratio")) return(1)
  if (identical(name, "allowDataUriScheme")) return(TRUE)
  NULL
}
registerS3method("$", "mockclientdata", cd_get, envir = asNamespace("shiny"))
registerS3method("[[", "mockclientdata", cd_get, envir = asNamespace("shiny"))
notes <- character(0)
assignInNamespace("showNotification", function(ui, ...) {
  txt <- if (is.character(ui)) ui else paste(as.character(ui), collapse = " ")
  notes <<- c(notes, gsub("<[^>]+>", "", txt))
  invisible("id")
}, "shiny")

# ---- reading defaults off the UI ---------------------------------------------------
ui_inputs <- function(html, prefix) {
  if (!nzchar(html)) return(list(values = list(), buttons = character(0)))
  doc <- xml2::read_html(paste0("<div>", html, "</div>"))
  strip <- function(id) sub(paste0("^(", prefix, "|proxy[0-9]+)-"), "", id)
  vals <- list()
  for (n in xml2::xml_find_all(doc, "//input[@id]")) {
    id <- xml2::xml_attr(n, "id")
    type <- xml2::xml_attr(n, "type")
    if (is.na(type)) type <- "text"
    if (type %in% c("file", "hidden", "radio")) next
    cls <- xml2::xml_attr(n, "class")
    if (!is.na(cls) && grepl("js-range-slider", cls)) {
      fr <- suppressWarnings(as.numeric(xml2::xml_attr(n, "data-from")))
      to <- suppressWarnings(as.numeric(xml2::xml_attr(n, "data-to")))
      vals[[strip(id)]] <- if (!is.na(to) && identical(xml2::xml_attr(n, "data-type"), "double")) c(fr, to) else fr
      next
    }
    if (type == "checkbox") vals[[strip(id)]] <- !is.na(xml2::xml_attr(n, "checked"))
    else if (type == "number") vals[[strip(id)]] <- suppressWarnings(as.numeric(xml2::xml_attr(n, "value")))
    else {
      v <- xml2::xml_attr(n, "value")
      vals[[strip(id)]] <- if (is.na(v)) "" else v
    }
  }
  for (n in xml2::xml_find_all(doc, "//select[@id]")) {
    id <- xml2::xml_attr(n, "id")
    opts <- xml2::xml_find_all(n, ".//option")
    v <- xml2::xml_attr(opts, "value")
    sel <- v[!is.na(xml2::xml_attr(opts, "selected"))]
    multiple <- !is.na(xml2::xml_attr(n, "multiple"))
    vals[[strip(id)]] <- if (length(sel)) sel else if (!multiple && length(v)) v[1] else NULL
  }
  for (n in xml2::xml_find_all(doc, "//textarea[@id]")) vals[[strip(xml2::xml_attr(n, "id"))]] <- xml2::xml_text(n)
  for (g in xml2::xml_find_all(doc, "//div[contains(@class,'shiny-input-radiogroup') or contains(@class,'shiny-input-checkboxgroup')][@id]")) {
    id <- xml2::xml_attr(g, "id")
    chk <- xml2::xml_find_all(g, ".//input[@checked]")
    v <- xml2::xml_attr(chk, "value")
    if (grepl("radiogroup", xml2::xml_attr(g, "class")) && !length(v)) v <- xml2::xml_attr(xml2::xml_find_first(g, ".//input"), "value")
    vals[[strip(id)]] <- v
  }
  btn <- xml2::xml_find_all(doc, "//button[@id][contains(@class,'action-button') or contains(@class,'bslib-task-button')]")
  list(values = vals, buttons = strip(xml2::xml_attr(btn, "id")))
}

# ---- one module -------------------------------------------------------------------
run_module <- function(key, server, ui_fn, rv, label) {
  ids <- module_output_ids(server)
  ui_html <- paste(as.character(ui_fn(key)), collapse = "\n")
  found <- character(0)
  rec <- function(stage, nm, msg) found <<- c(found, sprintf("%-28s %-9s %-14s %s", label, stage, nm, msg))
  check <- function(nm, stage) {
    for (sz in list(c(600, 400), c(320, 60), c(60, 30), c(8, 8))) {
      SIZE[["w"]] <<- sz[1]
      SIZE[["h"]] <<- sz[2]
      notes <<- character(0)
      d0 <- length(grDevices::dev.list())
      err <- tryCatch({ output[[nm]]; NULL }, error = function(e) e)
      d1 <- length(grDevices::dev.list())
      if (!is.null(err) && !expected_blank(err)) rec(stage, nm, sprintf("@%dx%d ERROR: %s", sz[1], sz[2], conditionMessage(err)))
      pe <- unique(grep("Plot error|作图出错", notes, value = TRUE))
      if (length(pe)) rec(stage, nm, sprintf("@%dx%d %s", sz[1], sz[2], substr(pe[1], 1, 160)))
      if (d1 != d0) rec(stage, nm, sprintf("@%dx%d DEVICES %d -> %d", sz[1], sz[2], d0, d1))
    }
  }
  output <- NULL
  suppressWarnings(suppressMessages(shiny::testServer(server, args = list(rv = rv, log_rv = shiny::reactiveVal(list())), {
    output <<- output
    session$flushReact()
    set <- function(html) {
      u <- ui_inputs(html, key)
      if (length(u$values)) do.call(session$setInputs, u$values)
      u$buttons
    }
    btns <- set(ui_html)
    if (length(OVERRIDE[[key]])) do.call(session$setInputs, OVERRIDE[[key]])
    for (round in 1:2) {                          # inputs created by renderUI
      for (nm in ids) {
        h <- tryCatch(output[[nm]], error = function(e) NULL)
        if (is.list(h) && !is.null(h$html)) btns <- unique(c(btns, set(as.character(h$html))))
      }
      session$flushReact()
    }
    if (length(OVERRIDE[[key]])) do.call(session$setInputs, OVERRIDE[[key]])
    session$flushReact()
    for (nm in ids) check(nm, "opened")
    why <- character(0)
    for (b in btns) {
      notes <<- character(0)
      d0 <- length(grDevices::dev.list())
      e <- tryCatch(withCallingHandlers({ do.call(session$setInputs, stats::setNames(list(1), b)); session$flushReact(); NULL },
                    warning = function(w) { notes <<- c(notes, conditionMessage(w)); invokeRestart("muffleWarning") }),
                    error = function(e) e)
      if (!is.null(e) && !expected_blank(e)) rec("click", b, paste("ERROR:", conditionMessage(e)))
      d1 <- length(grDevices::dev.list())
      if (d1 != d0) rec("click", b, sprintf("DEVICES %d -> %d", d0, d1))
      why <- c(why, notes)
    }
    for (nm in ids) check(nm, "after-run")
    ok <- isTRUE(shiny::isolate(rv$status[[key]]))
    done <<- c(done, sprintf("%-28s %s", label, if (ok) "DONE" else paste("not run:", substr(paste(unique(why), collapse = " / "), 1, 150))))
  })))
  found
}

# ---- data -------------------------------------------------------------------------
set.seed(1)
sc_raw <- suppressWarnings(as_seurat(readRDS(demo_bundled_path())))
sc_raw$sample <- rep(c("A", "B", "C", "D"), length.out = ncol(sc_raw))
sc_raw$condition <- ifelse(sc_raw$sample %in% c("A", "B"), "ctrl", "trt")
sc <- qc_add_metrics(sc_raw, "human")
sc <- normalize_obj(sc, "LogNormalize")
sc <- suppressWarnings(reduce_obj(sc, n_hvg = 300, npcs = 15))
sc <- suppressWarnings(cluster_obj(sc, dims = 10, resolutions = 0.8, algorithm = 1))
sc <- suppressWarnings(embed_obj(sc, "umap", dims = 10))
sc$celltype <- paste0("type", as.character(sc$seurat_clusters))
sc <- color_celltypes(sc)

all <- character(0)
done <- character(0)
# inputs a user would have to type (no usable default in the UI)
OVERRIDE <- list(wes_tmb = list(capture = 35.8), wes_tmbclin = list(capture = 35.8),
                 wes_lolli = list(gene = "DNMT3A"), wes_surv = list(genes = "DNMT3A"),
                 wes_hetero = list(sample = "TCGA-AB-2972"),
                 pseudobulk = list(sample = "sample", condition = "condition", group = "celltype",
                                   ref = "ctrl", alt = "trt", min_cells = 3),
                 abundance = list(sample = "sample", condition = "condition", group = "celltype"),
                 cluster = list(method = "louvain"), clinical = list(), viz = list(genes = "CD3E, LYZ"))
# (no package name here: the sync script renames it for scStudio)
steps <- if (exists("steps_for", pkgenv)) get("steps_for", pkgenv)("sc") else get("app_steps", pkgenv)()
# the module lists live in test-outputs.R: take the two list functions from it
src <- parse(file.path(repo, "tests/testthat/test-outputs.R"))
src_outputs <- new.env(parent = pkgenv)
for (e in src) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) && as.character(e[[2]]) %in% c("sc_module_servers", "wes_module_servers")) {
    eval(e, src_outputs)
  }
}
servers <- src_outputs$sc_module_servers()
for (s in steps) {
  if (!grepl(only, s$v) || is.null(servers[[s$v]])) next
  for (st in c("raw", "processed")) {
    obj <- if (st == "raw") sc_raw else sc
    rv <- shiny::reactiveValues(obj = obj, source = "qc", status = list(), clinical = NULL, ckpt = new.env())
    all <- c(all, run_module(s$v, servers[[s$v]], s$ui, rv, paste0("sc:", s$v, ":", st)))
  }
}
# a pass on real data (pbmc3k, processed) for the steps the synthetic demo cannot run
if (use_pbmc) {
  pb <- suppressWarnings(as_seurat(readRDS(app_sys("extdata", "pbmc3k.rds"))))
  pb <- qc_add_metrics(pb, "human")
  pb <- normalize_obj(pb, "LogNormalize")
  pb <- suppressWarnings(reduce_obj(pb, n_hvg = 2000, npcs = 20))
  pb <- suppressWarnings(cluster_obj(pb, dims = 10, resolutions = 0.5, algorithm = 1))
  pb <- suppressWarnings(embed_obj(pb, "umap", dims = 10))
  for (s in steps) {
    if (!grepl("cellcycle|malignancy|viz|markers|embed|cluster|annotate", s$v) || !grepl(only, s$v) || is.null(servers[[s$v]])) next
    rv <- shiny::reactiveValues(obj = pb, source = "qc", status = list(), clinical = NULL, ckpt = new.env())
    all <- c(all, run_module(s$v, servers[[s$v]], s$ui, rv, paste0("sc:", s$v, ":pbmc3k")))
  }
}
if (exists("wes_demo_paths", pkgenv) && requireNamespace("maftools", quietly = TRUE)) {
  p <- wes_demo_paths()
  clin <- wes_read_clinical(p$clinical, p$clinical)
  maf <- suppressMessages(wes_read_maf(p$maf, clinical = clin))
  ws <- src_outputs$wes_module_servers()
  for (s in steps_for("wes")) {
    if (!grepl(only, s$v) || is.null(ws[[s$v]])) next
    rv <- shiny::reactiveValues(omics = "wes", maf = maf, maf_source = "Demo: TCGA LAML", wes_clin_raw = clin,
                                wes_sequenced = wes_all_samples(maf), obj = NULL, status = list(),
                                clinical = NULL, epoch_wes = 1L, ckpt = new.env())
    all <- c(all, run_module(s$v, ws[[s$v]], s$ui, rv, paste0("wes:", s$v)))
  }
}
cat("\n==== coverage (step marked done after pressing its buttons) ====\n")
cat(done, sep = "\n")
cat("\n==== QC report ====\n")
cat("devices open at the end:", length(grDevices::dev.list()), "\n")
cat(unique(all), sep = "\n")
cat("\nissues:", length(unique(all)), "\n")
quit(status = if (length(all)) 1 else 0)
