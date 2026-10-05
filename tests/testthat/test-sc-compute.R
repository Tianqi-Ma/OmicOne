skip_if_not_installed("Seurat")
skip_if_not_installed("SeuratObject")

# Compute wrappers on the bundled 500-gene x 300-cell demo matrix: object
# standardisation, idempotent QC, the preview helpers, reduction bookkeeping,
# SCT handling, and that every reproducibility log is runnable R (the
# Seurat-4-compatible parts are replayed and must give the app's result).

demo_counts <- function() readRDS(demo_bundled_path())

demo_obj <- function() {
  obj <- suppressWarnings(as_seurat(demo_counts()))
  obj$sample <- rep(c("A", "B"), length.out = ncol(obj))
  obj
}

prepped <- function() {
  obj <- qc_add_metrics(demo_obj(), "human")
  obj <- normalize_obj(obj, "LogNormalize")
  suppressWarnings(reduce_obj(obj, n_hvg = 300, npcs = 20))
}

test_that("standardize_sc_object rebuilds an RNA assay with nCount_RNA (SCE-style object)", {
  counts <- demo_counts()
  sce_like <- suppressWarnings(Seurat::CreateSeuratObject(counts, assay = "originalexp"))
  expect_null(obj_meta(sce_like)$nCount_RNA)
  obj <- standardize_sc_object(sce_like)
  expect_true("RNA" %in% obj_assays(obj))
  expect_equal(obj_default_assay(obj), "RNA")
  expect_equal(unname(obj_meta(obj)$nCount_RNA), unname(Matrix::colSums(counts)))
  expect_equal(unname(obj_meta(obj)$nFeature_RNA), unname(Matrix::colSums(counts > 0)))
  acts <- obj_misc(obj, "omicone_import_actions")
  expect_true("rebuild_rna:originalexp" %in% acts)
  expect_silent(parse(text = import_log_code("seurat", acts)))
  # downstream QC now works on it
  expect_silent(qc_add_metrics(obj, "human"))
})

test_that("standardize_sc_object resets a non-RNA default assay (v4 integrated objects)", {
  obj <- suppressWarnings(Seurat::CreateSeuratObject(demo_counts()))
  obj[["integrated"]] <- SeuratObject::CreateAssayObject(counts = demo_counts())
  SeuratObject::DefaultAssay(obj) <- "integrated"
  out <- standardize_sc_object(obj)
  expect_equal(obj_default_assay(out), "RNA")
  expect_true("default_assay:integrated" %in% obj_misc(out, "omicone_import_actions"))
})

test_that("as_seurat converts a SingleCellExperiment from counts + colData", {
  skip_if_not(requireNamespace("SingleCellExperiment", quietly = TRUE),
              "SingleCellExperiment cannot be loaded")
  counts <- demo_counts()
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
  SummarizedExperiment::colData(sce)$donor <- rep(c("d1", "d2"), length.out = ncol(sce))
  obj <- as_seurat(sce)
  expect_equal(obj_default_assay(obj), "RNA")
  expect_false(is.null(obj_meta(obj)$nCount_RNA))
  expect_true("donor" %in% obj_meta_cols(obj))
})

test_that("QC is idempotent and its log replays to the same cells", {
  base <- demo_obj()
  run_qc <- function(o) {
    o <- qc_add_metrics(o, "human")
    keep <- qc_mad_keep(o, 3, 3, batch = "sample")
    list(obj = o[, as.vector(keep)], keep = keep)
  }
  a <- run_qc(base)
  b <- run_qc(base)
  expect_identical(as.vector(a$keep), as.vector(b$keep))
  expect_identical(colnames(a$obj), colnames(b$obj))
  # recomputing the metrics on their own output changes nothing
  md1 <- obj_meta(qc_add_metrics(base, "human"))
  md2 <- obj_meta(qc_add_metrics(qc_add_metrics(base, "human"), "human"))
  expect_equal(md1, md2)
  code <- qc_log_code("human", attr(a$keep, "thresholds"), "sample", "rule")
  e <- new.env()
  e$obj <- base
  eval(parse(text = code), envir = e)
  expect_identical(colnames(e$obj), colnames(a$obj))
})

test_that("the doublet preview keeps two fill levels when the class is a factor", {
  md <- obj_meta(demo_obj())
  set.seed(4)
  md$doublet_score <- runif(nrow(md))
  # scDblFinder returns a factor; ifelse() on it used to yield integer codes
  md$doublet_class <- factor(ifelse(md$doublet_score > 0.9, "doublet", "singlet"),
                             levels = c("singlet", "doublet"))
  pd <- doublet_plot_data(md, "sample")
  expect_equal(levels(pd$class), c("singlet", "doublet"))
  b <- ggplot2::ggplot_build(doublet_plot(pd, threshold = 0.9))
  expect_length(unique(b$data[[1]]$fill), 2)
  expect_equal(length(unique(b$layout$layout$PANEL)), 2)       # one panel per sample
  expect_true(any(vapply(b$plot$layers, function(l) inherits(l$geom, "GeomVline"), logical(1))))
  # a class vector without doublets still maps to the fixed levels
  md$doublet_class <- factor("singlet")
  expect_equal(levels(doublet_plot_data(md)$class), c("singlet", "doublet"))
})

test_that("graph_reduction_choices offers integrated.dr and follows the recorded default", {
  obj <- prepped()
  expect_equal(graph_reduction_choices(obj), "pca")
  emb <- SeuratObject::Embeddings(obj, "pca")
  obj[["integrated.dr"]] <- SeuratObject::CreateDimReducObject(
    embeddings = emb, key = "integrateddr_", assay = "RNA")
  obj[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = emb[, 1:2], key = "UMAP_", assay = "RNA")
  expect_equal(graph_reduction_choices(obj), c("pca", "integrated.dr"))
  obj@misc$omicone_reduction <- "integrated.dr"
  expect_equal(default_graph_reduction(obj), "integrated.dr")
  # integrate "none" removes the old integrated space and points back to PCA
  none <- integrate_obj(obj, "sample", "none")
  expect_false(has_reduction(none, "integrated.dr"))
  expect_equal(default_graph_reduction(none), "pca")
  # re-running PCA drops every reduction built on the old one
  again <- suppressWarnings(reduce_obj(obj, n_hvg = 300, npcs = 20))
  expect_equal(obj_reductions(again), "pca")
  expect_equal(default_graph_reduction(again), "pca")
})

test_that("Harmony becomes the default reduction and its log code is valid", {
  skip_if_not_installed("harmony")
  obj <- suppressWarnings(integrate_obj(prepped(), "sample", "harmony", dims = 15))
  expect_true(has_reduction(obj, "harmony"))
  expect_equal(default_graph_reduction(obj), "harmony")
  expect_silent(parse(text = integrate_log_code("harmony", "sample", 15)))
  skip_if_not_installed("RANN")
  pd <- integrate_plot_data(obj, "sample", "harmony", 15)
  expect_equal(names(pd$entropy), c("pca", "harmony"))
  expect_true(all(pd$entropy >= 0 & pd$entropy <= 1.05))
  expect_s3_class(ggplot2::ggplot_build(integrate_plot(pd, "sample")), "ggplot_built")
  expect_error(integrate_obj(prepped(), "orig.ident", "harmony"), "single level")
})

test_that("reduce keeps SCT's residuals and variable genes", {
  obj <- qc_add_metrics(demo_obj(), "human")
  sct <- tryCatch(suppressWarnings(normalize_obj(obj, "SCT")), error = function(e) NULL)
  skip_if(is.null(sct), "SCTransform cannot run here")
  expect_equal(obj_default_assay(sct), "SCT")
  before <- obj_layer(sct, "scale.data", "SCT")
  vf <- SeuratObject::VariableFeatures(sct)
  red <- suppressWarnings(reduce_obj(sct, n_hvg = 100, npcs = 15, hvg_method = "mvp"))
  expect_identical(obj_layer(red, "scale.data", "SCT"), before)
  expect_identical(SeuratObject::VariableFeatures(red), vf)
  expect_true(has_reduction(red, "pca"))
  expect_false(is.null(hvg_plot_data(red, "SCT")))
  # switching back to LogNormalize normalises RNA and drops the SCT assay
  back <- normalize_obj(red, "LogNormalize")
  expect_equal(obj_default_assay(back), "RNA")
  expect_false("SCT" %in% obj_assays(back))
})

test_that("clustering passes k, records every resolution and switches the active one", {
  obj <- prepped()
  out <- cluster_obj(obj, "pca", dims = 50, resolutions = c(0.3, 0.8), algorithm = 1,
                     k = 12, seed = 0)
  cols <- obj_misc(out, "omicone_cluster_cols")
  expect_equal(unname(cols), c("RNA_snn_res.0.3", "RNA_snn_res.0.8"))
  expect_equal(obj_misc(out, "omicone_cluster_col"), "RNA_snn_res.0.3")
  expect_equal(as.character(SeuratObject::Idents(out)),
               as.character(obj_meta(out)[["RNA_snn_res.0.3"]]))
  cmd <- out@commands[[grep("^FindNeighbors", names(out@commands))[1]]]
  expect_equal(cmd$k.param, 12)
  expect_equal(max(cmd$dims), 20)                      # capped at the 20 PCs available
  sw <- cluster_set_active(out, "RNA_snn_res.0.8")
  expect_equal(as.character(obj_meta(sw)$seurat_clusters),
               as.character(obj_meta(sw)[["RNA_snn_res.0.8"]]))
  expect_equal(obj_misc(sw, "omicone_cluster_col"), "RNA_snn_res.0.8")
  expect_error(cluster_set_active(out, "nope"), "not found")
})

test_that("each step's log code parses and the Seurat-4-compatible chain replays the app", {
  base <- demo_obj()
  # app path
  o <- qc_add_metrics(base, "human")
  keep <- qc_mad_keep(o, 5, 3, batch = "sample")
  o <- o[, as.vector(keep)]
  o <- normalize_obj(o, "LogNormalize", 1e4)
  dropped <- obj_reductions(o)
  o <- suppressWarnings(reduce_obj(o, n_hvg = 300, npcs = 20, hvg_method = "vst"))
  o <- cluster_obj(o, "pca", dims = 15, resolutions = c(0.5, 1), algorithm = 1, k = 20)
  pr <- embed_params(ncol(o), 20, 15, 30, 30)
  o <- suppressWarnings(embed_obj(o, "umap", "pca", dims = pr$dims, n_neighbors = pr$n_neighbors))
  codes <- list(
    qc = qc_log_code("human", attr(keep, "thresholds"), "sample", "rule"),
    normalize = normalize_log_code("LogNormalize", 1e4),
    reduce = reduce_log_code(FALSE, "vst", 300, 20, 42, dropped),
    cluster = cluster_log_code("pca", 15, c(0.5, 1), 1, 20, 0, "RNA", "RNA_snn_res.0.5"),
    embed = embed_log_code("umap", "pca", pr, 0.3))
  for (cd in codes) expect_silent(parse(text = cd))
  # replay the exported script body on the imported object
  e <- new.env()
  e$obj <- base
  # the script runs Seurat with its default (verbose) output; keep the test log quiet
  utils::capture.output(
    suppressMessages(suppressWarnings(eval(parse(text = unlist(codes)), envir = e))))
  expect_identical(colnames(e$obj), colnames(o))
  expect_equal(SeuratObject::Embeddings(e$obj, "pca"), SeuratObject::Embeddings(o, "pca"))
  expect_identical(as.character(SeuratObject::Idents(e$obj)),
                   as.character(SeuratObject::Idents(o)))
  expect_equal(SeuratObject::Embeddings(e$obj, "umap"), SeuratObject::Embeddings(o, "umap"))
  # the remaining steps' code (v5-only calls or optional packages) must still parse
  rest <- list(doublet_log_code("sample", 42, NULL, "remove"),
               integrate_log_code("RPCA", "sample", 15, 42, 99),
               markers_log_code("seurat_clusters"),
               annotate_singler_log_code("MonacoImmuneData"),
               import_log_code("h5", file = "f.h5"))
  for (cd in rest) expect_silent(parse(text = cd))
})

test_that("preview builders draw from real objects", {
  o <- qc_add_metrics(demo_obj(), "human")
  keep <- qc_mad_keep(o, 3, 3, batch = "sample")
  qp <- qc_plot(obj_meta(o), keep, attr(keep, "thresholds"), "sample")
  qb <- ggplot2::ggplot_build(qp)
  expect_equal(length(unique(qb$layout$layout$PANEL)), 3)        # genes / UMIs / mito
  o <- normalize_obj(o[, as.vector(keep)], "LogNormalize")
  dd <- normalize_diag_data(o)
  expect_equal(names(attr(dd, "rho")), c("before", "after"))
  expect_s3_class(ggplot2::ggplot_build(normalize_diag_plot(dd, "LogNormalize")), "ggplot_built")
  o <- suppressWarnings(reduce_obj(o, n_hvg = 300, npcs = 20))
  pv <- pca_variance(o)
  expect_equal(pv$basis[1], "total")
  expect_true(all(diff(pv$cum) >= 0) && max(pv$cum) <= 100.0001)
  expect_false(is.null(reduce_plot(hvg_plot_data(o, "vst"), pv)))
  o <- cluster_obj(o, "pca", dims = 15, resolutions = 0.8, algorithm = 1)
  expect_s3_class(cluster_plot(o, obj_misc(o, "omicone_cluster_col"), "0.8"), "ggplot")
  mk <- suppressWarnings(markers_obj(o, group_by = "seurat_clusters"))
  top <- marker_top_n(mk, 3)
  if (nrow(top)) {
    expect_s3_class(suppressWarnings(markers_dot_plot(o, top, "seurat_clusters")), "ggplot")
    expect_s3_class(ggplot2::ggplot_build(markers_bar_plot(top)), "ggplot_built")
  }
  o$celltype <- paste0("type", as.character(SeuratObject::Idents(o)))
  tab <- annotation_cluster_table(SeuratObject::Idents(o), o$celltype)
  expect_equal(sum(tab$n), ncol(o))
  expect_true(all(tab$share == 1))
  o <- suppressWarnings(embed_obj(o, "umap", "pca", dims = 15))
  expect_false(is.null(annotate_plot(o, tab, "manual")))
  expect_s3_class(embed_plot(o, "umap", "sample"), "ggplot")
})

test_that("embedding parameters are capped by the data", {
  pr <- embed_params(n_cells = 40, n_dims = 10, dims = 30, n_neighbors = 50, perplexity = 30)
  expect_equal(pr$dims, 10L)
  expect_equal(pr$n_neighbors, 39L)
  expect_lt(pr$perplexity * 3, 40 - 1)
})

test_that("the upstream modules run end to end on a real object and log runnable code", {
  rv <- shiny::reactiveValues(obj = NULL, source = NULL, status = list(), ckpt = new.env(),
                              clinical = NULL)
  log_rv <- shiny::reactiveVal(list())
  bad <- character(0)
  drive <- function(server, id, steps, outs) {
    suppressMessages(suppressWarnings(shiny::testServer(
      server, args = list(rv = rv, log_rv = log_rv), {
        session$flushReact()
        for (st in steps) {
          do.call(session$setInputs, st)
          session$flushReact()
        }
        for (nm in outs) {
          err <- tryCatch({
            output[[nm]]
            NULL
          }, error = function(e) e)
          if (!is.null(err) && !inherits(err, "shiny.silent.error")) {
            bad <<- c(bad, sprintf("%s$%s: %s", id, nm, conditionMessage(err)))
          }
        }
      })))
  }
  drive(mod_import_server, "import", list(list(demo_id = "bundled"), list(load_demo = 1)),
        c("summary", "ov_plot"))
  expect_equal(obj_dims(shiny::isolate(rv$obj))$cells, 300L)
  o <- shiny::isolate(rv$obj)
  o$sample <- rep(c("A", "B"), length.out = ncol(o))
  rv$obj <- o
  drive(mod_qc_server, "qc",
        list(list(species = "human", method = "mad", nmad_lib = 5, nmad_mt = 3,
                  max_mt_cap = NA, min_genes = 200, max_genes = 6000, max_mt = 15),
             list(batch = "sample"), list(run = 1)),
        c("summary", "insight", "preview", "batch_ui"))
  drive(mod_normalize_server, "normalize",
        list(list(method = "LogNormalize", scale_factor = 1e4), list(run = 1)),
        c("summary", "insight", "preview", "regress_ui"))
  drive(mod_reduce_server, "reduce",
        list(list(hvg_method = "vst", n_hvg = 300, npcs = 20.5), list(run = 1)),
        c("summary", "insight", "preview", "hvg_ui"))
  drive(mod_integrate_server, "integrate", list(list(method = "none", dims = 15), list(run = 1)),
        c("summary", "insight", "preview", "batch_ui"))
  drive(mod_cluster_server, "cluster",
        list(list(method = "louvain", dims = 15, neighbors = 20, resolutions = "0.5, 1, x"),
             list(reduction = "pca"), list(run = 1), list(active = "RNA_snn_res.1")),
        c("summary", "insight", "preview", "reduction_ui", "active_ui"))
  expect_equal(obj_misc(shiny::isolate(rv$obj), "omicone_cluster_col"), "RNA_snn_res.1")
  drive(mod_embed_server, "embed",
        list(list(method = "umap", dims = 15, n_neighbors = 30, min_dist = 0.3,
                  perplexity = 30, mask = FALSE),
             list(reduction = "pca"), list(run = 1), list(color_by = "sample")),
        c("summary", "preview", "reduction_ui", "color_ui"))
  drive(mod_markers_server, "markers",
        list(list(test = "wilcox", logfc = 0.25, min_pct = 0.1, only_pos = TRUE, top_n = ""),
             list(group_by = "seurat_clusters"), list(run = 1), list(top_n = 3)),
        c("summary", "insight", "preview", "bars", "tbl", "group_ui"))
  drive(mod_annotate_server, "annotate",
        list(list(method = "manual"), list(lab_1 = "T cell", lab_2 = "B/cell (x)"),
             list(apply_manual = 1)),
        c("summary", "insight", "preview", "manual_inputs", "tbl"))
  expect_equal(bad, character(0))
  st <- shiny::isolate(rv$status)
  for (k in c("import", "qc", "normalize", "reduce", "integrate", "cluster", "embed",
              "markers", "annotate")) {
    expect_true(isTRUE(st[[k]]), info = k)
  }
  expect_true("B/cell (x)" %in% shiny::isolate(rv$obj)$celltype)
  entries <- shiny::isolate(log_rv())
  expect_gte(length(entries), 9)
  for (e in entries) expect_silent(parse(text = e$code))
  # no OmicOne internals in the exported code
  all_code <- unlist(lapply(entries, function(e) e$code))
  expect_false(any(grepl("(run_doublets|integrate_obj|cluster_obj|embed_obj|keep_cells)\\(",
                         all_code)))
})
