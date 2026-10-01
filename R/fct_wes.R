#' WES / somatic mutation layer — wrappers around maftools
#'
#' Every function here is gated by `require_pkgs("maftools")` at the call site
#' and wrapped in tryCatch, so a missing package or a signature mismatch shows a
#' friendly message instead of crashing the app. maftools installs as a plain
#' Bioconductor binary, so unlike the single-cell engine it needs no source
#' build and runs on a locked-down Windows box.
#'
#' NOTE: maftools renamed a few entry points across versions (most notably
#' `OncogenicPathways()`/`PlotOncogenicPathways()` -> `pathways()`/
#' `plotPathways()` in 2.12). Where that happened the wrapper tries the new name
#' first and falls back to the old one, so both generations work.
#'
#' The maftools entry points this pipeline depends on, for checking against your
#' installed version on a first real run:
#' `read.maf`, `getSampleSummary`, `getGeneSummary`, `getClinicalData`,
#' `getFields`, `subsetMaf`, `plotmafSummary`, `oncoplot`, `titv`, `plotTiTv`,
#' `plotVaf`, `rainfallPlot`, `tmb`, `lollipopPlot`, `oncodrive`,
#' `plotOncodrive`, `somaticInteractions`, `pathways`, `trinucleotideMatrix`,
#' `extractSignatures`, `compareSignatures`, `plotSignatures`, `plotApobecDiff`,
#' `clinicalEnrichment`, `plotEnrichmentResults`, `drugInteractions`,
#' `mafCompare`, `forestPlot`, `coBarplot`, `inferHeterogeneity`, `plotClusters`.
#'
#' @name fct_wes
#' @keywords internal
NULL

#' Which of the maftools entry points this pipeline needs are actually present
#'
#' Used by the tests and useful at the console: a fast way to see whether an
#' installed maftools has drifted away from what the modules call.
#'
#' @return Character vector of missing function names (empty when all present);
#'   `NA_character_` when maftools itself is not installed.
#' @keywords internal
wes_missing_api <- function() {
  if (!has_pkg("maftools")) return(NA_character_)
  needed <- c("read.maf", "getSampleSummary", "getGeneSummary", "getClinicalData",
              "getFields", "subsetMaf", "plotmafSummary", "oncoplot", "titv",
              "plotTiTv", "plotVaf", "rainfallPlot", "tmb", "lollipopPlot",
              "oncodrive", "plotOncodrive", "somaticInteractions",
              "trinucleotideMatrix", "extractSignatures", "compareSignatures",
              "plotSignatures", "plotApobecDiff", "clinicalEnrichment",
              "plotEnrichmentResults", "drugInteractions", "mafCompare",
              "forestPlot", "coBarplot", "inferHeterogeneity", "plotClusters")
  ns <- asNamespace("maftools")
  missing <- needed[!vapply(needed, exists, logical(1), where = ns, inherits = FALSE)]
  # pathways() vs the pre-2.12 OncogenicPathways(): either is fine
  if (!exists("pathways", where = ns, inherits = FALSE) &&
      !exists("OncogenicPathways", where = ns, inherits = FALSE)) {
    missing <- c(missing, "pathways / OncogenicPathways")
  }
  missing
}

# ---- shared UI fragments ----------------------------------------------------

#' The "no MAF yet" placeholder every WES step shows before Import has run
#' @keywords internal
wes_no_maf <- function() {
  shiny::div(class = "omicone-placeholder",
             i18n("Load a MAF file on the <b>Import MAF</b> step first.",
                  "请先在<b>导入 MAF</b> 步骤加载 MAF 文件。"))
}

#' The "run this step" prompt shown before a WES step has produced anything
#' @param en,zh Prompt text.
#' @keywords internal
wes_prompt <- function(en, zh) {
  shiny::div(class = "omicone-placeholder", i18n(en, zh))
}

# ---- input ------------------------------------------------------------------

#' Read a MAF file into a maftools MAF object
#'
#' @param path Path to a `.maf` / `.maf.gz` / tab-delimited mutation table.
#' @param clinical Optional per-sample clinical data.frame; its first column
#'   must be `Tumor_Sample_Barcode`.
#' @param vc_nonSyn Optional character vector overriding which
#'   `Variant_Classification` values count as non-synonymous.
#' @return A maftools MAF object.
#' @keywords internal
wes_read_maf <- function(path, clinical = NULL, vc_nonSyn = NULL) {
  args <- list(maf = path)
  if (!is.null(clinical)) args$clinicalData <- clinical
  if (!is.null(vc_nonSyn)) args$vc_nonSyn <- vc_nonSyn
  do.call(maftools::read.maf, args)
}

#' Path to maftools' bundled TCGA LAML example (MAF + clinical annotation)
#'
#' Ships with maftools, so the WES pipeline has an instant offline demo in the
#' same spirit as the bundled pbmc3k on the single-cell side.
#'
#' @return list(maf=, clinical=); components are "" when maftools is absent.
#' @keywords internal
wes_demo_paths <- function() {
  list(
    maf = system.file("extdata", "tcga_laml.maf.gz", package = "maftools"),
    clinical = system.file("extdata", "tcga_laml_annot.tsv", package = "maftools")
  )
}

#' Sample ids in a MAF
#' @param maf A MAF object.
#' @keywords internal
wes_samples <- function(maf) {
  s <- tryCatch(as.character(maftools::getSampleSummary(maf)$Tumor_Sample_Barcode),
                error = function(e) character(0))
  s[!is.na(s)]
}

#' Gene names in a MAF, most-mutated first
#' @param maf A MAF object. @param n How many to return (Inf for all).
#' @keywords internal
wes_genes <- function(maf, n = Inf) {
  g <- tryCatch(as.character(maftools::getGeneSummary(maf)$Hugo_Symbol),
                error = function(e) character(0))
  if (is.finite(n)) utils::head(g, n) else g
}

#' Clinical columns carried inside the MAF object
#' @param maf A MAF object.
#' @keywords internal
wes_clinical_cols <- function(maf) {
  cd <- tryCatch(maftools::getClinicalData(maf), error = function(e) NULL)
  if (is.null(cd)) return(character(0))
  setdiff(colnames(cd), "Tumor_Sample_Barcode")
}

#' Columns of the MAF itself (used to find the VAF column)
#' @param maf A MAF object.
#' @keywords internal
wes_fields <- function(maf) {
  tryCatch(as.character(maftools::getFields(maf)), error = function(e) character(0))
}

#' Pick the variant-allele-frequency column out of a set of MAF field names
#'
#' MAFs disagree on this completely: TCGA uses `i_TumorVAF_WU`, others use
#' `VAF`, `tumor_vaf`, or only the read counts. Several steps (VAF plot,
#' heterogeneity) are useless without it, so guess and let the user override.
#' Split from [wes_guess_vaf_col()] so the guessing rule is testable on its own.
#'
#' @param fields Character vector of column names.
#' @return A column name, or NULL.
#' @keywords internal
wes_pick_vaf_col <- function(fields) {
  if (!length(fields)) return(NULL)
  exact <- c("i_TumorVAF_WU", "VAF", "vaf", "tumor_vaf", "TumorVAF", "t_vaf")
  hit <- fields[fields %in% exact]
  if (length(hit)) return(hit[1])
  hit <- grep("vaf|allele_freq", fields, ignore.case = TRUE, value = TRUE)
  if (length(hit)) return(hit[1])
  NULL
}

#' Pick the protein-change column out of a set of MAF field names
#' @param fields Character vector of column names.
#' @keywords internal
wes_pick_aa_col <- function(fields) {
  if (!length(fields)) return(NULL)
  exact <- c("HGVSp_Short", "AAChange", "Protein_Change", "amino_acid_change",
             "AAChange.refGene")
  hit <- fields[fields %in% exact]
  if (length(hit)) return(hit[1])
  hit <- grep("aachange|hgvsp|protein_change", fields, ignore.case = TRUE, value = TRUE)
  if (length(hit)) return(hit[1])
  NULL
}

#' Guess the VAF column of a MAF object
#' @param maf A MAF object.
#' @keywords internal
wes_guess_vaf_col <- function(maf) wes_pick_vaf_col(wes_fields(maf))

#' Guess the protein-change column of a MAF object
#' @param maf A MAF object.
#' @keywords internal
wes_guess_aa_col <- function(maf) wes_pick_aa_col(wes_fields(maf))

# ---- cohort summaries -------------------------------------------------------

#' Headline numbers for the summary tiles
#' @param maf A MAF object.
#' @return list(samples, genes, variants, median_per_sample, top_gene, top_pct).
#' @keywords internal
wes_overview <- function(maf) {
  ss <- tryCatch(maftools::getSampleSummary(maf), error = function(e) NULL)
  gs <- tryCatch(maftools::getGeneSummary(maf), error = function(e) NULL)
  n_samples <- if (!is.null(ss)) nrow(ss) else NA_integer_
  per_sample <- if (!is.null(ss) && "total" %in% colnames(ss)) ss$total else NA_real_
  list(
    samples  = n_samples,
    genes    = if (!is.null(gs)) nrow(gs) else NA_integer_,
    variants = if (!is.null(per_sample)) sum(per_sample, na.rm = TRUE) else NA_real_,
    median_per_sample = if (!is.null(per_sample)) stats::median(per_sample, na.rm = TRUE) else NA_real_,
    top_gene = if (!is.null(gs) && nrow(gs)) as.character(gs$Hugo_Symbol[1]) else NA_character_,
    top_pct  = if (!is.null(gs) && nrow(gs) && !is.na(n_samples) && n_samples > 0)
                 100 * gs$MutatedSamples[1] / n_samples else NA_real_
  )
}

#' Per-gene mutation frequency table
#' @param maf A MAF object. @param n Number of genes.
#' @return data.frame: Hugo_Symbol, MutatedSamples, pct, total.
#' @keywords internal
wes_gene_table <- function(maf, n = 50) {
  gs <- as.data.frame(maftools::getGeneSummary(maf))
  n_samples <- nrow(maftools::getSampleSummary(maf))
  gs <- utils::head(gs, n)
  keep <- intersect(c("Hugo_Symbol", "MutatedSamples", "total"), colnames(gs))
  out <- gs[, keep, drop = FALSE]
  if ("MutatedSamples" %in% keep && n_samples > 0) {
    out$pct <- round(100 * out$MutatedSamples / n_samples, 1)
  }
  out
}

# ---- burden -----------------------------------------------------------------

#' Tumour mutational burden per sample
#' @param maf A MAF object.
#' @param capture_size Exome capture size in Mb (used as the denominator).
#' @param log_scale Draw the plot on a log scale.
#' @return data.frame of per-sample TMB (maftools draws its own plot).
#' @keywords internal
wes_tmb <- function(maf, capture_size = 50, log_scale = TRUE) {
  as.data.frame(maftools::tmb(maf = maf, captureSize = capture_size,
                              logScale = log_scale))
}

#' Transition / transversion summary
#' @param maf A MAF object. @param use_syn Include synonymous variants.
#' @keywords internal
wes_titv <- function(maf, use_syn = TRUE) {
  maftools::titv(maf = maf, plot = FALSE, useSyn = use_syn)
}

# ---- drivers ----------------------------------------------------------------

#' Positional clustering driver detection (oncodrive)
#' @param maf A MAF object. @param aa_col Protein-change column.
#' @param min_mut Minimum mutations per gene to test.
#' @keywords internal
wes_oncodrive <- function(maf, aa_col = NULL, min_mut = 5) {
  maftools::oncodrive(maf = maf, AACol = aa_col, minMut = min_mut,
                      pvalMethod = "zscore")
}

#' Mutually exclusive / co-occurring gene pairs
#' @param maf A MAF object. @param top Number of top genes to test.
#' @keywords internal
wes_interactions <- function(maf, top = 25) {
  maftools::somaticInteractions(maf = maf, top = top, pvalue = c(0.05, 0.1))
}

#' Oncogenic pathway summary, across the maftools rename
#'
#' `OncogenicPathways()` became `pathways()` in maftools 2.12.
#' @param maf A MAF object.
#' @keywords internal
wes_pathways <- function(maf) {
  if (exists("pathways", where = asNamespace("maftools"), inherits = FALSE)) {
    return(maftools::pathways(maf = maf, plotType = "treemap"))
  }
  fn <- utils::getFromNamespace("OncogenicPathways", "maftools")
  fn(maf = maf)
}

# ---- signatures -------------------------------------------------------------

#' The BSgenome package a reference build needs
#' @param build "hg19" or "hg38".
#' @keywords internal
wes_bsgenome_pkg <- function(build = c("hg19", "hg38")) {
  build <- match.arg(build)
  if (build == "hg19") "BSgenome.Hsapiens.UCSC.hg19" else "BSgenome.Hsapiens.UCSC.hg38"
}

#' Trinucleotide context matrix (the input to signature extraction)
#' @param maf A MAF object. @param build "hg19"/"hg38".
#' @keywords internal
wes_trinuc <- function(maf, build = "hg19") {
  maftools::trinucleotideMatrix(maf = maf, ref_genome = wes_bsgenome_pkg(build),
                                prefix = "chr", add = TRUE)
}

#' Extract de-novo mutational signatures and match them to COSMIC
#' @param tnm A [wes_trinuc()] result. @param n Number of signatures.
#' @return list(sig=, cmp=) — the extracted signatures and the COSMIC comparison.
#' @keywords internal
wes_signatures <- function(tnm, n = 3) {
  # NMF (0.28, current CRAN) resolves its internal `seed` S4 generic with a
  # bare getGeneric("seed"), which returns NULL once the Bioconductor stack —
  # attached by trinucleotideMatrix's BSgenome load — is on the search path
  # while NMF itself is not. Attaching NMF restores the lookup.
  if (!"package:NMF" %in% search()) suppressMessages(attachNamespace("NMF"))
  # parallel = NULL: NMF's parallel mode (.opt = "P4") dies inside Shiny where
  # no foreach backend is registered; sequential NMF is seconds on this matrix.
  sig <- maftools::extractSignatures(mat = tnm, n = n, pConstant = 0.1,
                                     parallel = NULL)
  cmp <- tryCatch(maftools::compareSignatures(nmfRes = sig, sig_db = "SBS"),
                  error = function(e) NULL)
  list(sig = sig, cmp = cmp)
}

# ---- clinical / comparison --------------------------------------------------

#' Genes enriched in one level of a clinical feature
#' @param maf A MAF object. @param feature Clinical column name.
#' @keywords internal
wes_clinical_enrichment <- function(maf, feature) {
  maftools::clinicalEnrichment(maf = maf, clinicalFeature = feature)
}

#' Compare two cohorts defined by a clinical column's levels
#'
#' @param maf A MAF object.
#' @param feature Clinical column to split on.
#' @param level1,level2 The two levels to compare.
#' @param min_mut Minimum mutations for a gene to be tested.
#' @return list(m1=, m2=, res=, n1=, n2=) — the two sub-MAFs and the test.
#' @keywords internal
wes_compare_cohorts <- function(maf, feature, level1, level2, min_mut = 5) {
  cd <- as.data.frame(maftools::getClinicalData(maf))
  s1 <- as.character(cd$Tumor_Sample_Barcode[which(as.character(cd[[feature]]) == level1)])
  s2 <- as.character(cd$Tumor_Sample_Barcode[which(as.character(cd[[feature]]) == level2)])
  if (length(s1) < 2 || length(s2) < 2) {
    stop("Each group needs at least 2 samples (got ", length(s1), " and ",
         length(s2), ").")
  }
  m1 <- maftools::subsetMaf(maf = maf, tsb = s1, mafObj = TRUE)
  m2 <- maftools::subsetMaf(maf = maf, tsb = s2, mafObj = TRUE)
  res <- maftools::mafCompare(m1 = m1, m2 = m2, m1Name = level1, m2Name = level2,
                              minMut = min_mut)
  list(m1 = m1, m2 = m2, res = res, n1 = length(s1), n2 = length(s2))
}

# ---- survival ---------------------------------------------------------------

#' Per-sample mutation status of a gene set
#'
#' The bridge from a MAF to the shared survival layer: which samples carry a
#' mutation in any of `genes`. Every sample in the MAF is returned, so samples
#' with no mutation in the set are correctly "WT" rather than missing.
#'
#' @param maf A MAF object. @param genes Character vector of gene symbols.
#' @return data.frame: `.id`, `mutated` (logical), `status` ("Mutant"/"WT").
#' @keywords internal
wes_mutation_status <- function(maf, genes) {
  all_s <- wes_samples(maf)
  if (!length(all_s)) stop("No samples found in the MAF.")
  dat <- as.data.frame(maftools::subsetMaf(maf = maf, genes = genes,
                                           mafObj = FALSE))
  hit <- if (nrow(dat)) unique(as.character(dat$Tumor_Sample_Barcode)) else character(0)
  data.frame(
    .id     = all_s,
    mutated = all_s %in% hit,
    status  = ifelse(all_s %in% hit, "Mutant", "WT"),
    stringsAsFactors = FALSE
  )
}

# ---- heterogeneity ----------------------------------------------------------

#' Infer clonal structure of one sample from its VAF distribution
#' @param maf A MAF object. @param sample Tumor_Sample_Barcode.
#' @param vaf_col VAF column name.
#' @keywords internal
wes_heterogeneity <- function(maf, sample, vaf_col) {
  maftools::inferHeterogeneity(maf = maf, tsb = sample, vafCol = vaf_col)
}
