# WES layer: our own logic (guesses, id joins, survival analysis set, BED
# arithmetic, logged code) is tested without maftools; the wrappers around
# maftools are tested on its bundled TCGA LAML cohort (193 samples, 200
# clinical rows) and skipped when maftools is not installed.

laml <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      p <- wes_demo_paths()
      clin <- wes_read_clinical(p$clinical, p$clinical)
      maf <- suppressMessages(wes_read_maf(p$maf, clinical = clin))
      cache <<- list(maf = maf, clin = clin)
    }
    cache
  }
})

test_that("the WES step registry is internally consistent", {
  steps <- steps_wes()
  expect_equal(length(steps), 12)
  expect_equal(vapply(steps, function(s) s$n, numeric(1)), 1:12)
  expect_true(all(vapply(steps, function(s) is.function(s$ui), logical(1))))
  # keys are prefixed so they cannot collide with another omics' namespaces
  expect_true(all(grepl("^wes_", vapply(steps, function(s) s$v, character(1)))))
  # every phase key resolves to a label
  phases <- app_phases()
  for (s in steps) expect_true(!is.null(phases[[s$phase]]), info = s$v)
})

test_that("wes_missing_api reports cleanly when maftools is absent", {
  # On a machine without maftools this must be NA (not an error), so the check
  # is safe to call from anywhere.
  out <- wes_missing_api()
  if (!has_pkg("maftools")) {
    expect_true(is.na(out))
  } else {
    expect_equal(out, character(0))
  }
})

test_that("the VAF and protein-change column guesses handle real MAF namings", {
  tcga <- c("Hugo_Symbol", "Variant_Classification", "i_TumorVAF_WU",
            "HGVSp_Short", "Tumor_Sample_Barcode")
  expect_equal(wes_pick_vaf_col(tcga), "i_TumorVAF_WU")
  expect_equal(wes_pick_aa_col(tcga), "HGVSp_Short")

  # vcf2maf / annovar style
  expect_equal(wes_pick_vaf_col(c("Hugo_Symbol", "VAF")), "VAF")
  expect_equal(wes_pick_aa_col(c("Hugo_Symbol", "AAChange.refGene")),
               "AAChange.refGene")

  # unusual naming still gets found by the fuzzy pass
  expect_equal(wes_pick_vaf_col(c("Hugo_Symbol", "tumor_allele_freq")),
               "tumor_allele_freq")
  expect_equal(wes_pick_aa_col(c("Hugo_Symbol", "my_hgvsp_col")), "my_hgvsp_col")

  # nothing plausible -> NULL, so the module can say so instead of guessing wrong
  expect_null(wes_pick_vaf_col(c("Hugo_Symbol", "Chromosome")))
  expect_null(wes_pick_aa_col(c("Hugo_Symbol", "Chromosome")))
  expect_null(wes_pick_vaf_col(character(0)))
  expect_null(wes_pick_aa_col(character(0)))
})

test_that("population allele frequencies are never mistaken for the tumour VAF", {
  vep <- c("Hugo_Symbol", "ExAC_AF", "gnomAD_AF", "gnomAD_NFE_AF", "AF",
           "1000G_AF", "ESP_allele_freq", "Tumor_Sample_Barcode")
  expect_null(wes_pick_vaf_col(vep))
  expect_equal(wes_pick_vaf_col(c(vep, "tumor_f")), "tumor_f")
  # read counts let maftools compute the VAF itself
  expect_true(wes_vaf_auto(c("t_ref_count", "t_alt_count")))
  expect_false(wes_vaf_auto(c("t_alt_count")))
  vc <- wes_vaf_choice(c("t_ref_count", "t_alt_count"))
  expect_true(vc$auto)
  expect_equal(vc$selected, "")
})

test_that("ids are normalised for joining, with optional TCGA truncation", {
  expect_equal(wes_norm_id(c(" tcga-ab-2802 ", "S01", NA)), c("TCGA-AB-2802", "S01", NA))
  expect_equal(wes_norm_id("TCGA-AB-2802-03B-01W", tcga12 = TRUE), "TCGA-AB-2802")
  # non-TCGA ids are never truncated
  expect_equal(wes_norm_id("PATIENT-0001-TUMOUR", tcga12 = TRUE), "PATIENT-0001-TUMOUR")
  expect_equal(wes_norm_id(NULL), character(0))
})

test_that("the survival analysis set keeps unmatched rows as WT and collapses patients", {
  # five samples in the MAF; S6 and S7 are in the clinical table only
  st <- data.frame(.id = paste0("S", 1:5), in_maf = TRUE,
                   mutated = c(TRUE, TRUE, FALSE, FALSE, FALSE),
                   status = c("Mutant", "Mutant", "WT", "WT", "WT"),
                   stringsAsFactors = FALSE)
  clin <- data.frame(sample = c("s1", "S2 ", "S3", "S4", "S5", "S6", "S7"),
                     t = c(5, 8, 30, 40, Inf, 50, NA),
                     e = c(1, 1, 0, 1, 1, 0, 1), stringsAsFactors = FALSE)
  d <- wes_surv_data(clin, st, "sample", "t", "e", time_unit = "months")
  fl <- attr(d, "flow")
  expect_equal(fl$clin_rows, 7)
  expect_equal(fl$matched, 5)          # case and spaces ignored
  expect_equal(fl$unmatched, 2)
  expect_equal(fl$dropped_na, 2)       # Inf and NA follow-up
  expect_equal(nrow(d), 5)
  expect_equal(as.integer(table(d$.group)), c(3L, 2L))
  expect_true("S6" %in% d$.id[d$.group == "WT"])

  # without the WT assumption the unmatched patients leave the analysis
  d2 <- wes_surv_data(clin, st, "sample", "t", "e", time_unit = "months",
                      unmatched_wt = FALSE)
  expect_false("S6" %in% d2$.id)

  # two samples of one patient: mutant if either is
  clin_p <- data.frame(sample = c("S1", "S3", "S4"), patient = c("P1", "P1", "P2"),
                       t = c(5, 5, 9), e = c(1, 1, 0), stringsAsFactors = FALSE)
  dp <- wes_surv_data(clin_p, st, "sample", "t", "e", time_unit = "months",
                      patient_col = "patient")
  expect_equal(nrow(dp), 2)
  expect_equal(as.character(dp$.group[dp$.id == "P1"]), "Mutant")
  expect_equal(attr(dp, "flow")$dup_rows, 1)

  # repeated ids with different follow-up are a data error, not a silent pick
  clin_bad <- data.frame(sample = c("S1", "s1"), t = c(5, 9), e = c(1, 0),
                         stringsAsFactors = FALSE)
  expect_error(wes_surv_data(clin_bad, st, "sample", "t", "e"), "different follow-up")
})

test_that("the WES survival bridge produces a two-level grouping", {
  # wes_mutation_status() is the join between a MAF and fct_survival.R. Its
  # contract: every sample in the cohort comes back, labelled, so that samples
  # without a mutation are WT rather than dropped.
  st <- data.frame(
    .id     = paste0("S", 1:12),
    mutated = c(rep(TRUE, 5), rep(FALSE, 7)),
    status  = c(rep("Mutant", 5), rep("WT", 7)),
    stringsAsFactors = FALSE)

  clin <- normalise_clinical(
    data.frame(sample = paste0("S", 1:12),
               t = c(5, 8, 3, 6, 9,  40, 35, 50, 44, 60, 38, 47),
               e = c(1, 1, 1, 1, 1,   1,  1,  1,  1,  1,  1,  0),
               stringsAsFactors = FALSE),
    "sample", "t", "e")

  clin$.group <- factor(st$status[match(clin$.id, st$.id)], levels = c("WT", "Mutant"))
  expect_equal(as.integer(table(clin$.group)), c(7L, 5L))

  fit <- km_fit(clin)
  expect_s3_class(fit, "survfit")
  lr <- logrank_test(clin)
  expect_true(is.list(lr) && lr$df == 1)
  med <- km_medians(fit)
  expect_setequal(med$group, c("WT", "Mutant"))
  # the mutant arm dies early here, so it must have the shorter median
  expect_lt(med$median[med$group == "Mutant"], med$median[med$group == "WT"])
})

test_that("chromosome prefix and reference build are read off the data", {
  expect_false(wes_chr_has_prefix(c("1", "17", "X")))
  expect_true(wes_chr_has_prefix(c("chr1", "chr17", "chrX")))
  expect_false(wes_chr_has_prefix(character(0)))
  expect_equal(wes_harmonise_chr(c("chr1", "2"), prefix = TRUE), c("chr1", "chr2"))
  expect_equal(wes_harmonise_chr(c("chr1", "2"), prefix = FALSE), c("1", "2"))
  expect_equal(wes_build_from_ncbi(c("37", "37")), "hg19")
  expect_equal(wes_build_from_ncbi("GRCh38"), "hg38")
  expect_equal(wes_build_from_ncbi("hg19"), "hg19")
  expect_true(is.na(wes_build_from_ncbi(c(NA, ""))))
})

test_that("BED territory merges overlaps and counts 0-based half-open lengths", {
  f <- tempfile(fileext = ".bed")
  writeLines(c("track name=kit", "chr1\t100\t200", "chr1\t150\t300", "chr1\t300\t310",
               "chr2\t0\t50\textra\tcols"), f)
  bed <- wes_read_bed(f)
  expect_equal(nrow(bed), 2)
  # chr1: 100-310 after merging (abutting 300 joins), chr2: 0-50
  expect_equal(wes_bed_size_mb(bed), (210 + 50) / 1e6)
})

test_that("TCGA cohort labels never collide with a TCGA code", {
  codes <- c("ESCA", "LAML")
  expect_equal(wes_tcga_label("ESCA.maf", codes), "ESCA_input")
  expect_equal(wes_tcga_label("esca.maf.gz", codes), "esca_input")
  expect_equal(wes_tcga_label("Demo: TCGA LAML", codes), "TCGA LAML")
  expect_equal(wes_tcga_label(NULL, codes), "This cohort")
})

test_that("protein positions are parsed like lollipopPlot() does", {
  expect_equal(wes_aa_position(c("p.R882H", "p.Arg882His", "p.C229Lfs*18",
                                 "p.Asn1986GlnfsTer13", "p.761_762del")),
               c(882, 882, 229, 1986, 761))
  expect_true(all(is.na(wes_aa_position(c("", NA)))))
})

test_that("logged code is built from the arguments that ran", {
  expect_equal(wes_code("maftools::tmb", list(maf = quote(maf), captureSize = 35.8,
                                              logScale = TRUE), assign = "x"),
               "x <- maftools::tmb(maf = maf, captureSize = 35.8, logScale = TRUE)")
  # NULL = default and is omitted; wes_raw() passes literal code through
  expect_equal(wes_code("f", list(a = NULL, b = c("x", "y"), p = wes_raw("NULL"))),
               "f(b = c(\"x\", \"y\"), p = NULL)")
  code <- wes_surv_code(list(genes = "TP53", source = "maf", id_col = "id",
                             time_col = "t", event_col = "e", time_unit = "days",
                             patient_col = NULL, tcga12 = FALSE, unmatched_wt = FALSE,
                             raw_event = c(0, 1)))
  expect_true(any(grepl("includeSyn = FALSE", code)))
  expect_true(any(grepl("clin\\$.in_maf, \\]", code)))   # unmatched rows dropped
  expect_true(any(grepl("/ 30.4375", code)))
  expect_false(any(grepl(";", code, fixed = TRUE)))
})

test_that("small helpers survive empty or invalid numeric inputs", {
  expect_equal(wes_int(NA, 5, lo = 2), 5L)
  expect_equal(wes_int(NULL, 5), 5L)
  expect_equal(wes_int(100, 5, hi = 60), 60L)
  expect_equal(wes_prob(NA, 0.1), 0.1)
  expect_equal(wes_prob(1.5, 0.1), 0.1)
  expect_equal(wes_fmt(NA), "-")
  expect_equal(wes_fmt(1234.4), "1,234")
  expect_equal(wes_sep_for("clin.csv"), ",")
  expect_equal(wes_sep_for("clin.tsv"), "\t")
  expect_equal(wes_sep_for(NULL), "\t")
})

test_that("clinical ids keep leading zeros", {
  f <- tempfile(fileext = ".csv")
  writeLines(c("sample_id,age,os", "001,60,12.5", "002,71,3"), f)
  d <- wes_read_clinical(f, f, id_col = "sample_id")
  expect_identical(d$sample_id, c("001", "002"))
  expect_true(is.numeric(d$age))
  expect_equal(wes_find_col(names(d), "SAMPLE_ID"), "sample_id")
  expect_equal(unname(wes_maf_renames(c("hugo_symbol", "Chromosome", "tumor_sample_barcode"))),
               c("Hugo_Symbol", "Tumor_Sample_Barcode"))
})

test_that("the demo MAF path is resolvable in shape even without maftools", {
  p <- wes_demo_paths()
  expect_true(all(c("maf", "clinical") %in% names(p)))
  expect_true(is.character(p$maf) && is.character(p$clinical))
})

# ---- on maftools' TCGA LAML ---------------------------------------------------

test_that("mutation status counts non-synonymous mutants and keeps the universe", {
  skip_if_not_installed("maftools")
  d <- laml()
  st <- wes_mutation_status(d$maf, "TP53")
  expect_equal(nrow(st), 193)
  expect_equal(sum(st$mutated), 15)
  # the 7 clinical patients without any MAF record come back as WT, not missing
  st2 <- wes_mutation_status(d$maf, "TP53", universe = d$clin$Tumor_Sample_Barcode)
  expect_equal(nrow(st2), 200)
  expect_equal(sum(!st2$in_maf), 7)
  expect_equal(sum(st2$mutated), 15)
  expect_false(any(st2$mutated[!st2$in_maf]))
})

test_that("the LAML survival set has 188 patients, as maftools::mafSurvival", {
  skip_if_not_installed("maftools")
  d <- laml()
  st <- wes_mutation_status(d$maf, "DNMT3A", universe = d$clin$Tumor_Sample_Barcode)
  s <- wes_surv_data(d$clin, st, "Tumor_Sample_Barcode", "days_to_last_followup",
                     "Overall_Survival_Status", time_unit = "days")
  fl <- attr(s, "flow")
  expect_equal(nrow(s), 188)
  expect_equal(fl$unmatched, 7)
  expect_equal(fl$dropped_na, 12)       # non-finite follow-up
  expect_equal(sum(s$.group == "Mutant"), 45)
  expect_lt(logrank_test(s)$p, 0.01)
  # the old MAF-only set loses the 6 WT patients that read.maf dropped
  s2 <- wes_surv_data(d$clin, st, "Tumor_Sample_Barcode", "days_to_last_followup",
                      "Overall_Survival_Status", time_unit = "days", unmatched_wt = FALSE)
  expect_equal(nrow(s2), 182)

  # the logged code reproduces the same analysis set
  code <- wes_surv_code(list(genes = "DNMT3A", source = "maf",
                             id_col = "Tumor_Sample_Barcode",
                             time_col = "days_to_last_followup",
                             event_col = "Overall_Survival_Status", time_unit = "days",
                             patient_col = NULL, tcga12 = FALSE, unmatched_wt = TRUE,
                             raw_event = d$clin$Overall_Survival_Status))
  env <- new.env()
  env$maf <- d$maf
  env$clin_raw <- d$clin
  utils::capture.output(eval(parse(text = paste(code, collapse = "\n")), envir = env))
  expect_equal(nrow(env$clin), 188)
  expect_equal(sum(env$clin$.group == "Mutant"), 45)
})

test_that("TMB returns every sample, with the documented column names", {
  skip_if_not_installed("maftools")
  d <- laml()
  tm <- wes_tmb(d$maf, capture_size = 35.8)
  expect_true(all(c("Tumor_Sample_Barcode", "total", "total_perMB") %in% names(tm$df)))
  expect_equal(nrow(tm$df), 193)
  expect_equal(tm$n_zero, 1)
  st <- wes_tmb_stats(tm$df)
  expect_equal(st$median, 9 / 35.8)
  expect_equal(st$n_high, 0)
  expect_error(wes_tmb(d$maf, capture_size = NA), "positive")
  # a BED restricts the variants, sizes the territory, and keeps zero samples
  f <- tempfile(fileext = ".bed")
  writeLines(c("chr17\t7571719\t7590868", "chr2\t25457000\t25470000"), f)
  tb <- wes_tmb(d$maf, bed = wes_read_bed(f))
  expect_equal(nrow(tb$df), 193)
  expect_equal(tb$capture, (19149 + 13000) / 1e6)
  expect_gt(tb$n_zero, 100)
})

test_that("clinical enrichment counts genes, uses FDR and drops missing levels", {
  skip_if_not_installed("maftools")
  d <- laml()
  cd <- as.data.frame(maftools::getClinicalData(d$maf))
  expect_true("FAB_classification" %in% wes_enrichment_features(cd))
  enr <- suppressMessages(wes_clinical_enrichment(d$maf, "FAB_classification"))
  s <- wes_enrichment_sig(enr, 0.05)
  expect_equal(s$tested_genes, 20)
  expect_setequal(s$genes, c("IDH1", "TP53"))
  expect_equal(s$sig_genes, length(unique(s$sig$Hugo_Symbol)))
  # rows are gene x level; the gene count must not equal the row count by accident
  expect_gt(nrow(wes_enrichment_table(enr)), s$tested_genes)

  cd2 <- cd
  cd2$FAB_classification[1:10] <- NA
  m2 <- suppressMessages(maftools::read.maf(d$maf@data, clinicalData = cd2,
                                            verbose = FALSE))
  enr2 <- suppressMessages(wes_clinical_enrichment(m2, "FAB_classification"))
  cd_m2 <- as.data.frame(maftools::getClinicalData(m2))
  expect_gte(attr(enr2, "n_missing"), 10)
  expect_equal(attr(enr2, "n_missing"), sum(is.na(cd_m2$FAB_classification)))
  expect_equal(attr(enr2, "n_used") + attr(enr2, "n_missing"), nrow(cd_m2))
  # the missing samples are in no group: each level's "rest" excludes them
  tab <- wes_enrichment_table(enr2)
  rest_n <- as.numeric(sub(".* of ", "", tab$n_mutated_group2))
  grp_n <- as.numeric(sub(".* of ", "", tab$n_mutated_group1))
  expect_true(all(rest_n + grp_n == attr(enr2, "n_used")))
  expect_false(anyNA(enr2$groupwise_comparision$Group1))
})

test_that("cohort comparison reports the n tested and BH q-values", {
  skip_if_not_installed("maftools")
  d <- laml()
  out <- suppressMessages(wes_compare_cohorts(d$maf, "FAB_classification", "M2", "M3"))
  expect_equal(c(out$n1, out$n2), c(44, 21))
  fd <- wes_forest_data(out$res, 0.25, out$n1, out$n2)
  expect_true(all(fd$adjPval < 0.25))
  # NPM1 is 8 vs 0: infinite OR, drawn with a pseudo-count, flagged as such
  expect_true(fd$pseudo[fd$Hugo_Symbol == "NPM1"])
  expect_true(is.finite(fd$or_plot[fd$Hugo_Symbol == "NPM1"]))
  expect_s3_class(wes_forest_plot(fd, "M2", "M3", out$n1, out$n2, 0.25), "ggplot")
  expect_error(wes_forest_plot(fd[0, ], "M2", "M3", 44, 21, 0.05), "No gene passes")
})

test_that("driver, interaction, heterogeneity and lollipop helpers run on LAML", {
  skip_if_not_installed("maftools")
  skip_if_not_installed("mclust")
  d <- laml()
  od <- suppressWarnings(suppressMessages(
    wes_oncodrive(d$maf, aa_col = "Protein_Change", min_mut = 5)))
  expect_true("fdr" %in% names(od$res))
  # LAML has too few synonymous variants: the preset background is reported
  expect_true(length(od$bg_note) >= 1)
  it <- wes_interactions(d$maf, top = 10)
  expect_true(all(c("pAdj", "Event") %in% names(it)))

  het <- suppressMessages(wes_heterogeneity(d$maf, "TCGA-AB-2972", "i_TumorVAF_WU"))
  expect_equal(wes_clone_count(het$clusterData$cluster), 2)
  expect_equal(wes_clone_count(c("1", "2", "outlier", "CN_altered", "2")), 2)
  pc <- wes_math_percentile(wes_math_scores(d$maf, "i_TumorVAF_WU"), "TCGA-AB-2972")
  expect_true(pc$percentile > 0 && pc$percentile <= 100)

  rp <- wes_recurrent_positions(d$maf, "DNMT3A", "Protein_Change")
  expect_true(882 %in% rp$positions)
  expect_true(all(table(rp$all)[as.character(rp$positions)] >= 2))

  pooled <- wes_titv_pooled(wes_titv(d$maf))
  expect_equal(sum(pooled[c("C>A", "C>G", "C>T", "T>C", "T>A", "T>G")]), 100)
  expect_equal(wes_guess_build(d$maf), "hg19")
  expect_false(wes_chr_has_prefix(d$maf@data$Chromosome))
})

test_that("the logged import code rebuilds the same MAF", {
  skip_if_not_installed("maftools")
  code <- wes_import_code(maf_file = NULL, clin_file = "tcga_laml_annot.tsv", demo = TRUE,
                          id_from = "Tumor_Sample_Barcode",
                          text_cols = "Tumor_Sample_Barcode",
                          read_args = list(isTCGA = FALSE, rmFlags = FALSE))
  env <- new.env()
  utils::capture.output(suppressMessages(
    eval(parse(text = paste(code, collapse = "\n")), envir = env)))
  expect_equal(nrow(maftools::getSampleSummary(env$maf)), 193)
  expect_equal(nrow(env$clin_raw), 200)
})

test_that("WES steps are marked done only after a successful run, and never write rv$clinical", {
  skip_if_not_installed("maftools")
  d <- laml()
  shared <- data.frame(.id = "X1", .time = 1, .event = 1)
  run_step <- function(server, inputs) {
    rv <- shiny::reactiveValues(omics = "wes", maf = d$maf, maf_source = "Demo: TCGA LAML",
                                status = list(), clinical = shared, wes_clin_raw = d$clin,
                                epoch_wes = 1L, ckpt = new.env())
    log_rv <- shiny::reactiveVal(list())
    suppressWarnings(suppressMessages(shiny::testServer(server,
      args = list(rv = rv, log_rv = log_rv), {
        session$flushReact()
        do.call(session$setInputs, inputs)
        session$setInputs(run = 1)
        session$flushReact()
      })))
    list(status = shiny::isolate(rv$status), clinical = shiny::isolate(rv$clinical),
         log = shiny::isolate(log_rv()))
  }
  # no capture size and no BED: refused, not done
  r <- run_step(mod_wes_tmb_server, list(capture = NA, log = TRUE))
  expect_null(r$status$wes_tmb)
  r <- run_step(mod_wes_tmb_server, list(capture = 35.8, log = TRUE))
  expect_true(isTRUE(r$status$wes_tmb))
  # a gene list matching nothing: refused, not done
  r <- run_step(mod_wes_onco_server, list(gene_mode = "list", genes = "NOTAGENE"))
  expect_null(r$status$wes_onco)
  r <- run_step(mod_wes_onco_server, list(gene_mode = "list", genes = "tp53, FLT3"))
  expect_true(isTRUE(r$status$wes_onco))

  skip_if_not(exists("cox_hr") && exists("survival_warnings"))
  r <- run_step(mod_wes_surv_server,
                list(genes = "DNMT3A", endpoint = "Overall survival", src = "maf",
                     id_col = "Tumor_Sample_Barcode", time_col = "days_to_last_followup",
                     time_unit = "days", event_col = "Overall_Survival_Status",
                     patient_col = "", unmatched_wt = TRUE))
  expect_true(isTRUE(r$status$wes_surv))
  # the shared cohort is read, never overwritten by the WES analysis set
  expect_identical(r$clinical, shared)
  expect_true(any(grepl("survdiff", r$log[[1]]$code)))
})

test_that("signature extraction recovers known COSMIC signatures from a synthetic matrix", {
  # trinucleotideMatrix() needs a BSgenome; the NMF / COSMIC half is tested on
  # a 96-channel matrix simulated from SBS1, SBS4 and SBS13.
  skip_if_not_installed("maftools")
  skip_if_not_installed("NMF")
  db <- readRDS(system.file("extdata", "SBS_signatures.RDs", package = "maftools"))
  set.seed(1)
  w <- as.matrix(db$db[, c("SBS1", "SBS4", "SBS13")])
  h <- matrix(stats::rgamma(3 * 40, shape = 2, rate = 0.02), nrow = 3)
  cnt <- matrix(stats::rpois(96 * 40, lambda = as.vector(w %*% h)), nrow = 96)
  cnt[1:3, ] <- 0                          # empty channels, as in low-burden exomes
  tnm <- list(nmf_matrix = t(cnt))
  dimnames(tnm$nmf_matrix) <- list(paste0("S", 1:40), rownames(db$db))
  out <- suppressMessages(wes_signatures(tnm, n = 3))
  # NMF refuses all-zero channels: the retry adds pConstant = 1e-4, no more
  expect_equal(out$pconstant, 1e-4)
  m <- wes_sig_matches(out$cmp)
  expect_setequal(m$best_match, c("SBS1", "SBS4", "SBS13"))
  expect_true(all(m$cosine > 0.85))
  expect_equal(out$db, wes_sig_db())
})
