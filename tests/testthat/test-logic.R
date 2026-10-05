# Unit tests for the pure-logic helpers (no Seurat / no network required).
# These run anywhere, including CI without the heavy Suggests packages.

test_that("cell-type dictionary explains common labels (and tolerates variants)", {
  expect_match(explain_celltype("T cell"), "immune", ignore.case = TRUE)
  # tolerant of markers / plurals / case
  expect_match(explain_celltype("CD8+ T cells"), "[Cc]ytotoxic")
  expect_equal(explain_celltype("totally-unknown-type-xyz"), "")
})

test_that("guess_batch_col finds common sample columns", {
  expect_equal(guess_batch_col(data.frame(orig.ident = 1, x = 2)), "orig.ident")
  expect_equal(guess_batch_col(data.frame(sample = 1)), "sample")
  expect_null(guess_batch_col(data.frame(a = 1, b = 2)))
  expect_null(guess_batch_col(data.frame()))
})

test_that("memory_advice warns only for large datasets", {
  expect_equal(memory_advice(1e4), "")
  expect_gt(nchar(memory_advice(2e5)), 0)
  expect_gt(nchar(memory_advice(2e6)), 0)
})

test_that("demo dataset catalogue is well-formed", {
  ds <- demo_datasets()
  expect_s3_class(ds, "data.frame")
  expect_true(all(c("id", "name", "format", "source", "description") %in% names(ds)))
  expect_true(all(c("pbmc3k", "pancreas_sub", "bundled") %in% ds$id))
  expect_true(all(nzchar(ds$description)))
  # ids are unique
  expect_equal(anyDuplicated(ds$id), 0L)
})

test_that("MAD-based QC keeps normal cells and flags outliers", {
  set.seed(1)
  md <- data.frame(
    nCount_RNA   = c(rpois(60, 1000), 50000, 30),   # last two: huge / tiny
    nFeature_RNA = c(rpois(60, 500),  20,    15),
    percent.mt   = c(runif(60, 1, 8), 60,    55)     # last two: high mito
  )
  keep <- qc_mad_keep_from_meta(md, nmads_lib = 5, nmads_mt = 3)
  expect_length(keep, 62)
  expect_false(keep[61])          # extreme high-count / high-mito cell
  expect_false(keep[62])          # extreme low-count cell
  expect_gt(mean(keep[1:60]), 0.9)  # the bulk of normal cells survive
})

test_that("manual QC thresholds exclude out-of-range cells", {
  md <- data.frame(
    nCount_RNA   = c(1000, 1000, 1000),
    nFeature_RNA = c(500,  50,   500),   # 2nd below min
    percent.mt   = c(5,    5,    40)      # 3rd above max
  )
  keep <- qc_manual_keep_from_meta(md, min_genes = 100, max_genes = 2000, max_mt = 15)
  expect_equal(as.vector(keep), c(TRUE, FALSE, FALSE))
  expect_s3_class(attr(keep, "thresholds"), "data.frame")
})

test_that("read_counts_table parses a simple matrix", {
  tf <- tempfile(fileext = ".tsv")
  on.exit(unlink(tf))
  writeLines(c("gene\tc1\tc2", "GENE1\t3\t0", "GENE2\t1\t5"), tf)
  m <- read_counts_table(tf, sep = "\t")
  expect_equal(dim(m), c(2, 2))
  expect_equal(rownames(m), c("GENE1", "GENE2"))
  expect_equal(unname(m["GENE2", "c2"]), 5)
})

# ---- QC: per-sample MAD, MAD = 0, NA handling ---------------------------------

test_that("MAD thresholds can be computed per sample", {
  set.seed(2)
  # a small shallow sample (median ~600 UMIs) next to a large deep one (~4000)
  lib <- c(exp(rnorm(40, log(600), 0.3)), exp(rnorm(160, log(4000), 0.3)))
  md <- data.frame(
    nCount_RNA   = round(lib),
    nFeature_RNA = round(lib * 0.4),
    percent.mt   = runif(200, 1, 6),
    sample       = rep(c("shallow", "deep"), c(40, 160))
  )
  pooled <- qc_mad_keep_from_meta(md, 3, 3)
  per <- qc_mad_keep_from_meta(md, 3, 3, batch = "sample")
  thr <- attr(per, "thresholds")
  expect_setequal(thr$batch, c("shallow", "deep"))
  expect_equal(sum(thr$n), 200L)
  # pooled thresholds throw away the shallow sample; per-sample ones keep both
  expect_lt(mean(pooled[md$sample == "shallow"]), 0.5)
  expect_gt(mean(per[md$sample == "shallow"]), 0.9)
  expect_gt(mean(per[md$sample == "deep"]), 0.9)
})

test_that("a zero mitochondrial MAD disables the adaptive mito cut; the hard cap still applies", {
  md <- data.frame(nCount_RNA = rep(1000, 50), nFeature_RNA = rep(500, 50),
                   percent.mt = c(rep(0, 40), seq(1, 10, length.out = 10)))
  md$nCount_RNA <- md$nCount_RNA + seq_len(50)
  md$nFeature_RNA <- md$nFeature_RNA + seq_len(50)
  keep <- qc_mad_keep_from_meta(md, 5, 3)
  thr <- attr(keep, "thresholds")
  expect_true(grepl("percent.mt", thr$mad_zero))
  expect_equal(thr$mt_hi, Inf)
  expect_true(all(keep))                      # median 0, MAD 0 would otherwise flag every mt > 0
  capped <- qc_mad_keep_from_meta(md, 5, 3, max_mt = 5)
  expect_equal(sum(!capped), sum(md$percent.mt > 5))
})

test_that("cells with a missing metric or an unknown sample are never kept", {
  md <- data.frame(nCount_RNA = c(1000, NA, 1000, 1000), nFeature_RNA = c(500, 500, 500, 500),
                   percent.mt = c(2, 2, NA, 2), sample = c("a", "a", "a", NA))
  keep <- qc_manual_keep_from_meta(md, 100, 2000, 10)
  expect_equal(as.vector(keep), c(TRUE, FALSE, FALSE, TRUE))
  thr <- qc_mad_thresholds(md[c(1, 1, 1, 1), ], batch = "sample")
  per <- qc_apply_thresholds(md, thr, batch = "sample")
  expect_false(per[4])
})

test_that("the QC log code is valid R and reproduces the app's keep vector", {
  set.seed(3)
  md <- data.frame(nCount_RNA = rpois(80, 1500), nFeature_RNA = rpois(80, 700),
                   percent.mt = runif(80, 0, 12), sample = rep(c("x", "y"), 40))
  keep <- qc_mad_keep_from_meta(md, 3, 2, batch = "sample")
  code <- qc_log_code("human", attr(keep, "thresholds"), "sample", "test rule")
  expect_silent(parse(text = code))
  # replay the thresholding part of the script on the same metadata
  thr_lines <- code[grepl("^qc_thr|^i <-|^keep|^  md|^if \\(!is.null|^keep\\[", code)]
  e <- new.env()
  e$md <- md
  eval(parse(text = thr_lines), envir = e)
  expect_equal(e$keep, as.vector(keep))
})

# ---- QC gene patterns / ids --------------------------------------------------

test_that("hemoglobin patterns match globins only", {
  hs <- qc_gene_patterns("human")$hb
  expect_true(all(grepl(hs, c("HBB", "HBA1", "HBA2", "HBD", "HBG1", "HBG2", "HBE1",
                              "HBZ", "HBM", "HBQ1"))))
  expect_false(any(grepl(hs, c("HBEGF", "HBS1L", "HBP1", "HBCBP"))))
  mm <- qc_gene_patterns("mouse")$hb
  expect_true(all(grepl(mm, c("Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt", "Hbb-y", "Hbq1a"))))
  expect_false(any(grepl(mm, c("Hbegf", "Hbs1l", "Hbp1"))))
})

test_that("Ensembl row names are recognised (and their species read from the prefix)", {
  expect_equal(feature_id_type(c("ENSG00000141510", "ENSG00000146648", "TP53")), "ensembl")
  expect_equal(feature_id_type(c("TP53", "EGFR")), "symbol")
  m <- matrix(0, 3, 2, dimnames = list(c("ENSMUSG00000024747", "ENSMUSG00000059552",
                                         "ENSMUSG00000027490"), c("a", "b")))
  expect_equal(guess_species(m), "mouse")
})

# ---- Doublets ----------------------------------------------------------------

test_that("expected doublet rate follows scDblFinder's 0.8% per 1,000 cells", {
  expect_equal(doublet_expected_rate(10000), 0.08)
  expect_equal(doublet_call(c(0.1, 0.9), factor(c("singlet", "singlet")), 0.5),
               c("singlet", "doublet"))
  expect_equal(doublet_call(c(0.1, 0.9), factor(c("singlet", "doublet")), NA),
               c("singlet", "doublet"))
})

# ---- Integration helpers -----------------------------------------------------

test_that("batch columns are checked before integration", {
  md <- data.frame(b = c("a", "a", "b", NA), one = "x", stringsAsFactors = FALSE)
  expect_error(check_batch_col(md, "b"), "missing")
  expect_error(check_batch_col(md, "one"), "single level")
  expect_error(check_batch_col(md, "nope"), "not found")
  expect_true(check_batch_col(md[1:3, ], "b"))
})

test_that("categorical_cols never offers continuous metrics or single-level columns as a batch", {
  md <- data.frame(sample = rep(c("s1", "s2"), 10), nCount_RNA = runif(20, 500, 5000),
                   percent.mt = runif(20), orig.ident = "one", lane = rep(1:2, 10),
                   partly = c(NA, rep("a", 9), rep("b", 10)), stringsAsFactors = FALSE)
  expect_setequal(categorical_cols(md), c("sample", "lane"))
})

test_that("QC metrics are never offered as a sample column, even when integer-valued", {
  md <- data.frame(nFeature_RNA = rep(c(200L, 210L, 220L), 10), nCount_RNA = rep(1:3, 10),
                   percent.mt = rep(c(1, 2, 3), 10), S.Score = rep(c(0, 1), 15),
                   lane = rep(c("L1", "L2", "L3"), 10), stringsAsFactors = FALSE)
  expect_false("nFeature_RNA" %in% categorical_cols(md))  # metrics are excluded by name
  expect_equal(sc_group_cols(md), "lane")
})

test_that("a rebuilt reduction selector keeps the user's choice unless the default moved", {
  expect_equal(reduction_selection("pca", c("pca", "harmony"), "harmony", "pca"), "harmony")
  expect_equal(reduction_selection("pca", c("pca", "harmony"), "harmony", "harmony"), "pca")
  expect_equal(reduction_selection("harmony", "pca", "pca", "pca"), "pca")
})

# ---- Markers / annotation helpers --------------------------------------------

test_that("top-N markers are filtered on adjusted p, ranked by fold change, uniquely labelled", {
  df <- data.frame(gene = c("A", "B", "C", "A", "D"),
                   cluster = factor(c("0", "0", "0", "1", "1")),
                   avg_log2FC = c(3, 5, 1, 2, 1),
                   p_val_adj = c(0.001, 0.2, 0.01, 0.001, 0.04))
  top <- marker_top_n(df, 2)
  expect_equal(top$gene[top$cluster == "0"], c("A", "C"))   # B is big but not significant
  expect_equal(anyDuplicated(top$label), 0L)
  expect_true(all(c("A___0", "A___1") %in% top$label))
  roc <- data.frame(gene = c("x", "y"), cluster = "0", power = c(0.2, 0.8))
  expect_equal(marker_top_n(roc, 1)$gene, "y")
})

test_that("the Azimuth label column follows the requested level", {
  cols <- c("predicted.celltype.l1", "predicted.celltype.l1.score",
            "predicted.celltype.l2", "predicted.celltype.l2.score", "mapping.score")
  expect_equal(azimuth_label_col(cols, "2"), "predicted.celltype.l2")
  expect_equal(azimuth_label_col(cols, "3"), "predicted.celltype.l1")
  expect_equal(azimuth_label_col(c("predicted.ann_level_1", "predicted.ann_level_3"), "3"),
               "predicted.ann_level_3")
  expect_null(azimuth_label_col("nCount_RNA"))
  expect_equal(singler_ref_species("ImmGenData"), "mouse")
  expect_equal(singler_ref_species("MonacoImmuneData"), "human")
})

test_that("cell-type explanations match whole words and prefer the longest key", {
  expect_match(explain_celltype("Mast cells"), "histamine")       # not "T cell"
  expect_match(explain_celltype("CD14 Mono"), "Classical")
  expect_match(explain_celltype("CD16 Mono"), "Non-classical")
  expect_match(explain_celltype("Treg"), "Treg")
  expect_match(explain_celltype("CD8 TEM"), "Effector-memory cytotoxic")
  expect_match(explain_celltype("NK"), "Natural killer")
  expect_match(explain_celltype("pDC"), "Plasmacytoid")
  expect_match(explain_celltype("Plasma"), "antibody")
  expect_match(explain_celltype("Fibroblasts"), "extracellular")
  expect_match(explain_celltype("Endothelial_cells"), "vessels")
  expect_match(explain_celltype("Epithelial"), "lining")
  expect_match(explain_celltype("Tcm/Naive helper T cells"), "Helper T cell")
  expect_equal(explain_celltype(c(NA, "")), c("", ""))
})

# ---- Inputs, tables, log helpers ---------------------------------------------

test_that("numeric inputs are coerced before formatting", {
  expect_identical(int_input(30.5), 30L)
  expect_identical(int_input(NA), NA_integer_)
  expect_identical(int_input(NULL), NA_integer_)
  expect_identical(int_input(-3, lo = 1), NA_integer_)
  expect_identical(num_input("abc"), NA_real_)
  expect_equal(num_input(0.4, 0, 1), 0.4)
})

test_that("table separators follow the file extension and non-numeric columns are refused", {
  expect_equal(guess_table_sep("counts.csv"), ",")
  expect_equal(guess_table_sep("counts.csv.gz"), ",")
  expect_equal(guess_table_sep("counts.tsv"), "\t")
  expect_equal(guess_table_sep("counts.txt"), "\t")
  tf <- tempfile(fileext = ".csv")
  on.exit(unlink(tf))
  writeLines(c("gene,c1,c2", "GENE1,3,0", "GENE2,1,5"), tf)
  expect_equal(dim(read_counts_table(tf)), c(2, 2))
  writeLines(c("gene,symbol,c1", "G1,TP53,3", "G2,EGFR,1"), tf)
  expect_error(read_counts_table(tf), "non-numeric")
})

test_that("log-code builders emit valid R", {
  thr <- qc_manual_thresholds(200, 6000, 15)
  codes <- list(
    qc_log_code("mouse", thr, NULL, "manual"),
    doublet_log_code("sample", 42, 0.4, "remove"),
    doublet_log_code(NULL, 42, NULL, "flag"),
    normalize_log_code("LogNormalize", 1e4, NULL, TRUE),
    normalize_log_code("SCT", 1e4, c("percent.mt", "S.Score")),
    reduce_log_code(FALSE, "vst", 2000, 50, 42, c("pca", "harmony", "umap")),
    reduce_log_code(FALSE, "mvp", 2000, 50),
    reduce_log_code(TRUE, "SCT", 2000, 30),
    integrate_log_code("none"),
    integrate_log_code("harmony", "sample", 30),
    integrate_log_code("CCA", "sample", 30, 42, 50),
    cluster_log_code("harmony", 30, c(0.5, 1), 4, 20, 0, "RNA", "RNA_snn_res.0.5", TRUE),
    embed_log_code("umap", "pca", list(dims = 30L, n_neighbors = 30L), 0.3),
    embed_log_code("tsne", "pca", list(dims = 30L, perplexity = 30), 0.3),
    markers_log_code("celltype", "wilcox", 0.25, 0.1, TRUE),
    markers_log_code(NULL, "roc", 0.25, 0.1, FALSE),
    annotate_singler_log_code("HumanPrimaryCellAtlasData", "fine"),
    annotate_azimuth_log_code("pbmcref", "predicted.celltype.l2"),
    annotate_manual_log_code(c("0" = "T cell", "1" = "B cell's")),
    import_log_code("seurat", c("default_assay:integrated", "join_layers")),
    import_log_code("sce"),
    import_log_code("h5", file = "x.h5"),
    import_log_code("table", sep = ","),
    import_log_code("matrix", c("rebuild_rna:originalexp", "add_ncount"),
                    demo_file = "pbmc3k.rds"),
    import_log_code("scop_demo")
  )
  for (cd in codes) expect_silent(parse(text = cd))
})
