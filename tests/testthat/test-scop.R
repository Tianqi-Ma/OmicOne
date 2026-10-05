# The scop engine layer (R/fct_scop.R) and the downstream single-cell steps.
#
# 1. API guard: every `pkg::f(arg = )` call in fct_scop.R must name an
#    exported function and an argument it has (or reach a `...`). scop is
#    only on the user's machine, so these tests skip without it; they are the
#    check that catches the next signature drift (group_by vs group.by).
# 2. Pure helpers: label maps, coverage checks, column detection, input
#    checks -- no Seurat / scop needed.
# 3. Reproducibility code: every *_log_code() output parses and calls only
#    base R or namespaced functions, never OmicOne internals.

# Recurse into the parts of a call, skipping empty arguments (x[i, ]).
walk_args <- function(e, walk, skip_head = FALSE) {
  parts <- as.list(e)
  idx <- seq_along(parts)
  if (skip_head) idx <- idx[-1]
  for (i in idx) {
    if (is.null(parts[[i]])) next
    if (is.symbol(parts[[i]]) && !nzchar(as.character(parts[[i]]))) next
    walk(parts[[i]])
  }
}

scop_src <- function() {
  f <- testthat::test_path("..", "..", "R", "fct_scop.R")
  if (!file.exists(f)) testthat::skip("R/fct_scop.R not available (installed package)")
  f
}

# Every pkg::fun(...) call in a file: data.frame(pkg, fn, args).
ns_calls <- function(path, pkgs) {
  rows <- list()
  walk <- function(e) {
    if (is.call(e)) {
      h <- e[[1]]
      if (is.call(h) && identical(h[[1]], as.name("::")) &&
          as.character(h[[2]]) %in% pkgs) {
        an <- names(as.list(e))[-1]
        an <- an[!is.na(an) & nzchar(an)]
        rows[[length(rows) + 1]] <<- data.frame(
          pkg = as.character(h[[2]]), fn = as.character(h[[3]]),
          args = paste(an, collapse = ","), stringsAsFactors = FALSE)
      }
      walk_args(e, walk)
    }
  }
  for (e in parse(path, keep.source = FALSE)) walk(e)
  if (!length(rows)) return(data.frame(pkg = character(0), fn = character(0), args = character(0)))
  unique(do.call(rbind, rows))
}

# Formals a call can match: the function's, or its Seurat S3 method's.
call_formals <- function(pkg, fn) {
  ns <- asNamespace(pkg)
  f <- get(fn, envir = ns)
  fm <- names(formals(f))
  m <- utils::getS3method(fn, "Seurat", optional = TRUE, envir = ns)
  if (!is.null(m)) fm <- names(formals(m))
  fm
}

check_pkg_calls <- function(pkg) {
  testthat::skip_if_not_installed(pkg)
  calls <- ns_calls(scop_src(), pkg)
  testthat::skip_if(!nrow(calls), paste("no", pkg, "calls"))
  exports <- getNamespaceExports(pkg)
  for (i in seq_len(nrow(calls))) {
    fn <- calls$fn[i]
    testthat::expect_true(fn %in% exports, info = sprintf("%s::%s is not exported", pkg, fn))
    if (!fn %in% exports) next
    fm <- call_formals(pkg, fn)
    args <- strsplit(calls$args[i], ",", fixed = TRUE)[[1]]
    for (a in args) {
      testthat::expect_true(a %in% fm || "..." %in% fm,
                            info = sprintf("%s::%s has no argument '%s'", pkg, fn, a))
    }
  }
}

test_that("scop calls in fct_scop.R match the installed scop", {
  check_pkg_calls("scop")
})

test_that("UCell / CellChat / liana / copykat / mascarade calls match", {
  for (p in c("UCell", "CellChat", "liana", "copykat", "mascarade")) {
    if (requireNamespace(p, quietly = TRUE)) check_pkg_calls(p)
  }
  succeed()
})

test_that("fct_scop.R only uses scop functions that the 0.9.2 API has", {
  calls <- ns_calls(scop_src(), "scop")
  expect_false(any(calls$fn %in% c("palette_scp", "Standard_SCP", "Integration_SCP",
                                   "RunDimReduction", "VolcanoPlot")))
  expect_false(any(grepl("(^|,)group_by(,|$)", calls$args)))
})

# ---- copykat ----------------------------------------------------------------

test_that("copykat labels map to malignant / normal / NA with a confidence", {
  pred <- c("aneuploid", "diploid", "c1:diploid:low.conf", "c2:aneuploid:low.conf",
            "not.defined", NA)
  m <- copykat_label_map(pred)
  expect_equal(m$call, c("malignant", "normal", "normal", "malignant", NA, NA))
  expect_equal(m$confidence, c("high", "high", "low", "low", NA, NA))
  expect_equal(copykat_genome("mouse"), "mm10")
  expect_equal(copykat_genome("human"), "hg20")
})

test_that("reference cells called malignant are measured", {
  calls <- c(a = "normal", b = "malignant", c = "normal", d = NA, e = "malignant")
  expect_equal(copykat_ref_flagged(calls, c("a", "b", "c", "d")), 1 / 3)
  expect_true(is.na(copykat_ref_flagged(calls, character(0))))
})

# ---- mascarade --------------------------------------------------------------

test_that("mascarade outlines come back with x / y columns whatever the embedding names", {
  skip_if_not_installed("mascarade")
  set.seed(1)
  emb <- rbind(matrix(stats::rnorm(400, 0, 0.5), ncol = 2),
               matrix(stats::rnorm(400, 4, 0.5), ncol = 2))
  colnames(emb) <- c("UMAP_1", "UMAP_2")
  rownames(emb) <- paste0("c", seq_len(nrow(emb)))
  labels <- rep(c("A", "B"), each = 200)
  labels[c(1, 250)] <- NA
  mask <- mascarade_mask(emb, labels)
  expect_true(all(c("x", "y", "group", "cluster") %in% colnames(mask)))
  expect_setequal(unique(as.character(mask$cluster)), c("A", "B"))
})

# ---- gene sets / scoring ------------------------------------------------------

test_that("gene-set coverage is counted and too-small sets are refused", {
  cov <- geneset_coverage(list(a = c("X", "Y", "Z", "Q"), b = c("X", "W")), c("X", "Y", "Z"))
  expect_equal(cov$n_input, c(4L, 2L))
  expect_equal(cov$n_found, c(3L, 1L))
  expect_equal(cov$found[[1]], c("X", "Y", "Z"))
  expect_error(geneset_require(cov, 3), "b \\(1 of 2\\)")
  expect_silent(geneset_require(cov[1, ], 3))
})

test_that("mouse case conversion and safe set names", {
  expect_equal(mouse_case(c("MKI67", "TOP2A")), c("Mki67", "Top2a"))
  expect_equal(safe_set_names(c("my set", "2nd-set", "my set")),
               c("my_set", "S2nd_set", "my_set_1"))
})

test_that("custom gene sets are parsed with parse_genes()", {
  sets <- parse_custom_sets("T cell: CD3D, CD3E;TRAC\n\nbad line\nB: MS4A1 CD79A")
  expect_equal(names(sets), c("T_cell", "B"))
  expect_equal(sets$T_cell, c("CD3D", "CD3E", "TRAC"))
  expect_equal(sets$B, c("MS4A1", "CD79A"))
  expect_length(parse_custom_sets(""), 0)
})

test_that("AddModuleScore columns are renamed to <set>_AMS", {
  map <- ams_col_map(c("EMT", "hypoxia"))
  expect_equal(names(map), c("EMT_AMS", "hypoxia_AMS"))
  expect_equal(unname(map), c("omicone_ams1", "omicone_ams2"))
  expect_equal(modulescore_cols("EMT", "UCell"), "EMT_UCell")
  expect_equal(modulescore_cols("EMT", "AddModuleScore"), "EMT_AMS")
})

# ---- enrichment ------------------------------------------------------------------

test_that("ORA background is restricted to tested genes", {
  t2g <- data.frame(Term = c("t1", "t1", "t2", "t2", "t2"),
                    symbol = c("A", "B", "B", "C", NA),
                    entrez_id = 1:5, stringsAsFactors = FALSE)
  out <- enrich_filter_term2gene(t2g, universe = c("A", "B"))
  expect_equal(colnames(out), c("Term", "symbol"))
  expect_equal(nrow(out), 3)
  expect_true(all(out$symbol %in% c("A", "B")))
})

test_that("enrichment database keys follow scop 0.9.2 PrepareDB names", {
  expect_true(all(c("GO_BP", "GO_CC", "GO_MF", "KEGG", "Reactome", "WikiPathway",
                    "MSigDB_H") %in% enrich_db_choices()))
  expect_equal(enrich_db_key("MSigDB_H", "Mus_musculus"), "MSigDB_MH")
  expect_equal(enrich_db_key("MSigDB_H", "Homo_sapiens"), "MSigDB_H")
  expect_equal(enrich_db_key("KEGG", "Mus_musculus"), "KEGG")
  expect_equal(enrich_threshold(0.25, 0.05), "avg_log2FC > 0.25 & p_val_adj < 0.05")
})

test_that("enrichment summary counts BH-significant sets and groups", {
  tab <- data.frame(Groups = c("1", "1", "2", "3"), p.adjust = c(0.01, 0.2, 0.04, NA))
  s <- enrichment_summary(tab)
  expect_equal(s$n_tested, 4L)
  expect_equal(s$n_sig, 2L)
  expect_equal(s$groups, 3L)
  expect_equal(s$groups_hit, 2L)
  expect_equal(enrichment_summary(data.frame())$n_sig, 0L)
})

# ---- trajectory / dynamics / velocity ----------------------------------------------

test_that("trajectory start is validated per method", {
  lv <- c("HSC", "Mono")
  expect_null(traj_check_start("slingshot", "", lv))
  expect_null(traj_check_start("slingshot", "HSC", lv))
  expect_match(traj_check_start("slingshot", "Bcell", lv)$en, "not a level")
  for (m in c("monocle3", "paga", "palantir")) {
    expect_match(traj_check_start(m, "", lv)$en, "needs a start group")
    expect_null(traj_check_start(m, "HSC", lv, reductions = c("pca", "umap")))
    expect_match(traj_check_start(m, "HSC", lv, reductions = "pca")$en, "UMAP")
  }
  expect_match(traj_check_start("wot", "", lv, time_col = "day", time_values = rep(1, 5))$en,
               "two time points")
  expect_null(traj_check_start("wot", "", lv, time_col = "day", time_values = c(0, 3, 7)))
  expect_match(traj_check_start("wot", "", lv, time_col = "day",
                                time_values = c("d0", "d3"))$en, "numeric")
})

test_that("only this run's pseudotime columns are picked up", {
  new <- c("Slingshot_Lineage1", "Slingshot_Lineage2", "Slingshot_BranchID")
  expect_equal(traj_pt_cols(new, "slingshot"), c("Slingshot_Lineage1", "Slingshot_Lineage2"))
  expect_equal(traj_pt_cols(c("Monocle3_clusters", "Monocle3_Pseudotime"), "monocle3"),
               "Monocle3_Pseudotime")
  expect_equal(traj_pt_cols(c("dpt_pseudotime"), "paga"), "dpt_pseudotime")
  expect_length(traj_pt_cols(c("trajectory_0"), "wot"), 0)
  expect_length(traj_pt_cols(character(0), "palantir"), 0)
  expect_true(grepl(traj_output_pattern("slingshot"), "Slingshot_Lineage1"))
  expect_false(grepl(traj_output_pattern("slingshot"), "Lineage1"))
  expect_error(traj_output_pattern("tscan"), "Unknown")
})

test_that("time columns for WOT are numeric with a few distinct values", {
  md <- data.frame(day = rep(c(0, 3, 7), 40), score = seq(0, 1, length.out = 120),
                   label = rep(c("a", "b"), 60), stringsAsFactors = FALSE)
  expect_equal(traj_time_cols(md), "day")
})

test_that("dynamic features find lineage columns or report none", {
  md <- data.frame(Slingshot_Lineage1 = 1:3, Lineage2 = c(1, NA, 2),
                   seurat_clusters = c("0", "1", "1"),
                   Monocle3_Pseudotime = c(0.1, 0.2, 0.3), nCount_RNA = 1:3,
                   stringsAsFactors = FALSE)
  expect_setequal(dynamic_lineage_cols(md),
                  c("Slingshot_Lineage1", "Lineage2", "Monocle3_Pseudotime"))
  expect_equal(dynamic_lineage_cols(md, recorded = "Monocle3_Pseudotime"), "Monocle3_Pseudotime")
  expect_length(dynamic_lineage_cols(md["nCount_RNA"]), 0)
  expect_error(sc_dynamic(NULL, character(0)), "Trajectory step first")
})

test_that("velocity prerequisites are reported", {
  expect_setequal(velocity_missing("RNA", c("pca")), c("spliced", "unspliced", "umap"))
  expect_length(velocity_missing(c("RNA", "spliced", "unspliced"), c("pca", "umap")), 0)
})

# ---- grouping helpers / communication ----------------------------------------------

test_that("default grouping prefers annotation, then the active clustering", {
  expect_equal(default_group_col(c("orig.ident", "seurat_clusters", "celltype")), "celltype")
  expect_equal(default_group_col(c("orig.ident", "RNA_snn_res.1", "seurat_clusters"),
                                 cluster_col = "RNA_snn_res.1"), "RNA_snn_res.1")
  expect_equal(default_group_col(c("orig.ident")), "orig.ident")
  expect_null(default_group_col(character(0)))
  md <- data.frame(g = factor(c("b", "a", NA), levels = c("b", "a", "z")))
  expect_equal(group_levels(md, "g"), c("b", "a"))
  expect_length(group_levels(md, "missing"), 0)
})

test_that("CellChat labels never contain the reserved '0'", {
  expect_equal(cellchat_labels(c("0", "1", "10")), c("C0", "C1", "C10"))
  expect_equal(cellchat_labels(c("T", "B")), c("T", "B"))
  expect_equal(liana_resource("mouse"), "MouseConsensus")
  expect_equal(liana_resource("human"), "Consensus")
})

test_that("LIANA significant pairs use aggregate_rank < 0.05", {
  res <- data.frame(source = c("A", "A", "B"), target = c("B", "A", "A"),
                    ligand.complex = c("L1", "L2", "L3"), receptor.complex = c("R1", "R2", "R3"),
                    aggregate_rank = c(0.01, 0.2, 0.04), stringsAsFactors = FALSE)
  tab <- cellcomm_sig_table(res, "liana")
  expect_equal(nrow(tab), 2)
  expect_equal(tab$ligand, c("L1", "L3"))
  p <- cellcomm_count_plot(tab)
  expect_s3_class(p, "ggplot")
  expect_s3_class(cellcomm_count_plot(tab[0, ]), "ggplot")
})

# ---- reproducibility code -------------------------------------------------------------

# Function names a code snippet calls without a namespace.
bare_calls <- function(code) {
  out <- character(0)
  walk <- function(e) {
    if (is.call(e)) {
      h <- e[[1]]
      if (is.name(h)) out <<- c(out, as.character(h))
      walk_args(e, walk, skip_head = is.call(h) && identical(h[[1]], as.name("::")))
    }
  }
  for (e in parse(text = code, keep.source = FALSE)) walk(e)
  unique(out)
}

base_fun <- function(f) {
  any(vapply(c("base", "stats", "utils", "methods"), function(p)
    exists(f, envir = asNamespace(p), inherits = FALSE), logical(1)))
}

test_that("every log-code builder emits runnable code with no OmicOne internals", {
  codes <- list(
    detest_log_code("seurat_clusters"),
    enrichment_log_code("ora", "seurat_clusters", "GO_BP", "Homo_sapiens", lfc = 0.5),
    enrichment_log_code("gsea", "celltype", "MSigDB_H", "Mus_musculus"),
    velocity_log_code("celltype", "stochastic"),
    dynamic_log_code(c("Slingshot_Lineage1", "Slingshot_Lineage2"), 1000),
    cellcycle_log_code("human"),
    cellcycle_log_code("mouse"),
    modulescore_log_code(list(EMT = c("VIM", "ZEB1", "FN1")), "UCell"),
    modulescore_log_code(list(EMT = c("VIM", "ZEB1", "FN1"), hyp = c("CA9", "VEGFA", "LDHA")),
                         "AddModuleScore"),
    stemness_log_code(c("SOX2", "MYC", "NANOG")),
    cellchat_log_code("seurat_clusters", "mouse", prefixed = TRUE),
    liana_log_code("celltype", "human"),
    copykat_log_code("celltype", c("T cells", "B cells"), "orig.ident", "mouse"),
    copykat_log_code(NULL, NULL, NULL, "human")
  )
  for (m in c("slingshot", "monocle2", "monocle3", "paga", "palantir", "wot")) {
    codes[[length(codes) + 1]] <- trajectory_log_code(m, "celltype", start = "HSC",
                                                      reduction = "harmony", nd = 20,
                                                      time_col = "day")
  }
  codes[[length(codes) + 1]] <- trajectory_log_code("slingshot", "celltype", start = "")
  for (code in codes) {
    expect_true(is.character(code) && length(code) > 0)
    expect_silent(parse(text = code))
    bare <- bare_calls(code)
    ok <- vapply(bare, base_fun, logical(1))
    expect_true(all(ok), info = paste("not base R:", paste(bare[!ok], collapse = ", ")))
  }
  expect_match(paste(trajectory_log_code("slingshot", "g", start = "", nd = 20), collapse = "\n"),
               "start = NULL")
  expect_match(paste(enrichment_log_code("ora", "g", "KEGG", "Homo_sapiens"), collapse = "\n"),
               "TERM2GENE = term2gene")
})

# ---- report / export ----------------------------------------------------------------------

test_that("report sections match log entries by step key, not by step name", {
  entries <- list(
    list(step = "Features / PCA", key = "reduce"),
    list(step = "Dynamic features", key = "dynamic"),
    list(step = "Markers", key = "markers"),
    list(step = "Something new", key = "new_step"),
    list(step = "Legacy entry"))
  kept <- report_filter_entries(entries, c("trajectory"))
  expect_equal(vapply(kept, `[[`, "", "step"),
               c("Dynamic features", "Something new", "Legacy entry"))
  kept <- report_filter_entries(entries, c("reduce", "markers"))
  expect_equal(vapply(kept, `[[`, "", "step"),
               c("Features / PCA", "Markers", "Something new", "Legacy entry"))
  keys <- unlist(lapply(report_sections(), `[[`, "keys"))
  sc_keys <- vapply(steps_sc(), function(s) s$v, character(1))
  analysis <- setdiff(sc_keys, c("import", "viz", "report", "export"))
  expect_true(all(analysis %in% keys), info = paste(setdiff(analysis, keys), collapse = ", "))
})

test_that("the exported script loads the packages it uses and records the session", {
  entries <- list(
    list(step = "Normalization", key = "normalize", params = list(method = "LogNormalize"),
         code = "obj <- Seurat::NormalizeData(obj)"),
    list(step = "Trajectory", key = "trajectory", params = list(),
         code = c("obj <- scop::RunSlingshot(obj, group.by = \"g\", show_plot = FALSE)",
                  "m <- Matrix::rowSums(x); y <- utils::head(m)")))
  expect_equal(log_code_packages(entries), c("Seurat", "Matrix", "scop"))
  txt <- export_script_text(entries, session = c("R version 4.5.1", "scop: 0.9.2"))
  expect_match(txt, "library(Seurat)", fixed = TRUE)
  expect_match(txt, "library(scop)", fixed = TRUE)
  expect_false(grepl("library(utils)", txt, fixed = TRUE))
  expect_match(txt, "#   scop: 0.9.2", fixed = TRUE)
  expect_silent(parse(text = txt))
  info <- session_info_lines("Matrix")
  expect_true(any(grepl("^scop: ", info)))
  expect_true(any(grepl("^Matrix: ", info)))
  expect_true(any(grepl("sessionInfo|R version", info)))
})
