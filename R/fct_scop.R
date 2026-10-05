#' scop engine: compute and plotting wrappers
#'
#' Everything that goes through the scop package (mengxu98/scop) -- or through
#' the packages the advanced single-cell steps call directly (UCell, CellChat,
#' LIANA, copykat, mascarade) -- lives here, so the API surface that drifts
#' between releases is in one file. Every `scop::f(arg = )` call below was
#' checked against the scop 0.9.2 sources (GitHub main, 2026-09-12); the
#' installed version's formals are checked again by tests/testthat/test-scop.R.
#'
#' Next to each compute wrapper sits the builder of its reproducibility code
#' (`*_log_code()`): native scop / Seurat / CellChat / liana / copykat calls with
#' the parameters actually used, never OmicOne internals, so the exported
#' script runs without this package.
#'
#' @name fct_scop
#' @keywords internal
NULL

# ---- Shared helpers ----------------------------------------------------------

#' scop's species label for the working object
#'
#' scop's database, cell-cycle and CNV helpers take `"Homo_sapiens"` /
#' `"Mus_musculus"`; the object only tells us the gene-symbol casing.
#' @param obj A Seurat object.
#' @return "Homo_sapiens" or "Mus_musculus".
#' @keywords internal
scop_species <- function(obj) {
  if (identical(guess_species(obj), "mouse")) "Mus_musculus" else "Homo_sapiens"
}

#' A per-run scratch directory under the session's tempdir()
#'
#' scop's Python bridges and copykat write plots, h5ad files and transport maps
#' into a directory (default: the working directory). Pointing them at the
#' session tempdir keeps the user's folders clean and is removed with the
#' R session.
#' @param name Sub-directory name.
#' @return The directory path (created).
#' @keywords internal
scop_outdir <- function(name) {
  d <- file.path(tempdir(), "omicone_scop", name)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

#' Default grouping column for a step that groups cells
#'
#' Prefers an annotation column, then the active clustering column recorded by
#' the Cluster step, then `seurat_clusters`.
#' @param cols Candidate column names (already filtered to categorical ones).
#' @param cluster_col The active clustering column (`obj@misc$omicone_cluster_col`).
#' @return One column name, or NULL when `cols` is empty.
#' @keywords internal
default_group_col <- function(cols, cluster_col = NULL) {
  pref <- c("celltype", "cell_type", "CellType", "SingleR", cluster_col, "seurat_clusters")
  hit <- intersect(pref, cols)
  if (length(hit)) return(hit[1])
  if (length(cols)) cols[1] else NULL
}

#' Levels of a metadata column, in the order figures show them
#' @param md Cell metadata. @param col Column name.
#' @return Character vector (empty when the column is missing); see [level_order()].
#' @keywords internal
group_levels <- function(md, col) {
  if (is.null(col) || !col %in% colnames(md)) return(character(0))
  level_order(md[[col]])
}

#' Coverage of gene sets in the genes of the object
#' @param sets Named list of character vectors.
#' @param genes Genes present in the object (rownames).
#' @return data.frame: set, n_input, n_found, and a list column `found`.
#' @keywords internal
geneset_coverage <- function(sets, genes) {
  found <- lapply(sets, function(g) intersect(unique(g), genes))
  out <- data.frame(set = names(sets),
                    n_input = vapply(sets, function(g) length(unique(g)), integer(1)),
                    n_found = vapply(found, length, integer(1)),
                    stringsAsFactors = FALSE)
  out$found <- unname(found)
  rownames(out) <- NULL
  out
}

#' Stop when a gene set has too few genes in the object
#'
#' A module score on one or two genes is a gene-expression value with extra
#' noise, not a programme score.
#' @param cov Output of [geneset_coverage()].
#' @param min_genes Minimum genes found per set.
#' @return `cov`, invisibly, when every set passes.
#' @keywords internal
geneset_require <- function(cov, min_genes = 3) {
  bad <- cov[cov$n_found < min_genes, , drop = FALSE]
  if (nrow(bad)) {
    stop(sprintf("Gene set(s) with fewer than %d genes found in the data: %s.",
                 min_genes,
                 paste(sprintf("%s (%d of %d)", bad$set, bad$n_found, bad$n_input),
                       collapse = ", ")),
         call. = FALSE)
  }
  invisible(cov)
}

#' Human gene symbols in mouse casing (MKI67 -> Mki67)
#'
#' A case conversion, not an orthology mapping: it is right for most 1:1
#' orthologues and wrong for genes whose mouse symbol differs. The coverage
#' check that follows reports how many genes were found.
#' @param genes Character vector of human symbols.
#' @keywords internal
mouse_case <- function(genes) {
  paste0(substr(genes, 1, 1), tolower(substring(genes, 2)))
}

#' Safe metadata names for gene sets
#'
#' Scores are written to metadata columns named after the set; spaces and
#' punctuation would break FetchData() later.
#' @param x Character vector of set names.
#' @keywords internal
safe_set_names <- function(x) {
  x <- gsub("[^A-Za-z0-9_.]+", "_", trimws(x))
  x <- sub("^([0-9_.])", "S\\1", x)
  make.unique(x, sep = "_")
}

# ---- Differential expression feeding enrichment / GSEA ---------------------

#' Parameters of the DE table enrichment and GSEA read
#'
#' One-vs-rest Wilcoxon per group, both directions, no fold-change pre-filter
#' (`fc.threshold = 1`), BH within each group. Without both directions and the
#' full tested set, GSEA would rank a truncated list and ORA would have no
#' honest background.
#' @param group_by Grouping column.
#' @return A string recorded on the DE table, so a later run can tell whether
#'   the stored table was made with these parameters for this grouping.
#' @keywords internal
de_fingerprint <- function(group_by) {
  paste0("omicone:", group_by, ":wilcox;only.pos=FALSE;fc.threshold=1;BH")
}

#' Make sure the object carries a DE table enrichment / GSEA can read
#'
#' scop's RunEnrichment() / RunGSEA() read
#' `obj@tools[["DEtest_<group.by>"]][["AllMarkers_wilcox"]]`, which only
#' scop::RunDEtest() writes (Seurat::FindAllMarkers() does not).
#' @param srt Seurat object. @param group_by Grouping column.
#' @param reuse Reuse a stored table made by this function (the caller passes
#'   TRUE only when nothing upstream changed since).
#' @return The object with the DE table.
#' @keywords internal
sc_ensure_detest <- function(srt, group_by, reuse = FALSE) {
  key <- paste0("DEtest_", group_by)
  de <- srt@tools[[key]][["AllMarkers_wilcox"]]
  if (isTRUE(reuse) && is.data.frame(de) && nrow(de) &&
      identical(attr(de, "omicone"), de_fingerprint(group_by))) {
    return(srt)
  }
  srt <- scop::RunDEtest(srt, group.by = group_by, test.use = "wilcox",
                         only.pos = FALSE, fc.threshold = 1,
                         p.adjust.method = "BH", verbose = FALSE)
  de <- srt@tools[[key]][["AllMarkers_wilcox"]]
  if (!is.data.frame(de) || !nrow(de)) {
    stop("scop::RunDEtest() returned no DE results for '", group_by, "'.", call. = FALSE)
  }
  attr(de, "omicone") <- de_fingerprint(group_by)
  srt@tools[[key]][["AllMarkers_wilcox"]] <- de
  srt
}

#' R code reproducing [sc_ensure_detest()]
#' @param group_by Grouping column.
#' @keywords internal
detest_log_code <- function(group_by) {
  sprintf(paste0("obj <- scop::RunDEtest(obj, group.by = %s, test.use = \"wilcox\", ",
                 "only.pos = FALSE, fc.threshold = 1, p.adjust.method = \"BH\")"),
          r_lit(group_by))
}

# ---- Enrichment (ORA) and GSEA ---------------------------------------------

#' Gene-set databases offered by the Enrichment step
#'
#' Names as scop 0.9.2's PrepareDB() expects them (R/PrepareDB.R: `db`).
#' @return Named character vector (label -> db key).
#' @keywords internal
enrich_db_choices <- function() {
  c("GO biological process" = "GO_BP",
    "GO cellular component" = "GO_CC",
    "GO molecular function" = "GO_MF",
    "KEGG" = "KEGG",
    "Reactome" = "Reactome",
    "WikiPathways" = "WikiPathway",
    "MSigDB Hallmark" = "MSigDB_H")
}

#' The PrepareDB() key of a database for one species
#'
#' MSigDB collections are keyed by their own collection code; the mouse
#' Hallmark collection is "MH", the human one "H".
#' @param db A value of [enrich_db_choices()]. @param species scop species.
#' @keywords internal
enrich_db_key <- function(db, species) {
  if (identical(db, "MSigDB_H") && identical(species, "Mus_musculus")) return("MSigDB_MH")
  db
}

#' Restrict a TERM2GENE table to the background genes
#'
#' scop's ORA treats every gene annotated in the database as the background,
#' which inflates enrichment for any list drawn from the few thousand genes a
#' single-cell assay detects. Filtering the gene-set table to the genes that
#' were tested makes the hypergeometric universe = annotated tested genes,
#' the clusterProfiler `universe =` convention.
#' @param term2gene data.frame with a `Term` column and an `idtype` column.
#' @param universe Character vector of background genes.
#' @param idtype Gene-ID column. @return Two-column data.frame (Term, idtype).
#' @keywords internal
enrich_filter_term2gene <- function(term2gene, universe, idtype = "symbol") {
  t2g <- term2gene[, c("Term", idtype), drop = FALSE]
  t2g <- t2g[!is.na(t2g[[idtype]]) & t2g[[idtype]] %in% universe, , drop = FALSE]
  t2g <- unique(t2g)
  rownames(t2g) <- NULL
  t2g
}

#' Over-representation analysis per group with a tested-gene background
#'
#' @param srt Seurat object. @param group_by Grouping column.
#' @param db A value of [enrich_db_choices()]. @param species scop species.
#' @param lfc Minimum avg log2 fold-change for a gene to enter a group's list.
#' @param padj BH-adjusted p-value cut-off for that list.
#' @param reuse_de Passed to [sc_ensure_detest()].
#' @return The object; the result sits in
#'   `srt@tools[["Enrichment_<group_by>_wilcox"]]`, with the database name in
#'   its `Database` column and the background size in `$universe_size`.
#' @keywords internal
sc_enrichment <- function(srt, group_by, db = "GO_BP", species = "Homo_sapiens",
                          lfc = 0.25, padj = 0.05, reuse_de = FALSE) {
  srt <- sc_ensure_detest(srt, group_by, reuse = reuse_de)
  de <- srt@tools[[paste0("DEtest_", group_by)]][["AllMarkers_wilcox"]]
  universe <- unique(as.character(de$gene))
  key <- enrich_db_key(db, species)
  db_list <- scop::PrepareDB(species = species, db = key, db_IDtypes = "symbol")
  entry <- db_list[[species]][[key]]
  if (is.null(entry$TERM2GENE)) {
    stop("scop::PrepareDB() returned no gene sets for ", key, " (", species, ").",
         call. = FALSE)
  }
  t2g <- enrich_filter_term2gene(entry$TERM2GENE, universe, "symbol")
  if (!nrow(t2g)) {
    stop("None of the tested genes is annotated in ", key, " for ", species,
         "; check the species.", call. = FALSE)
  }
  t2n <- entry$TERM2NAME
  res_key <- paste("Enrichment", group_by, "wilcox", sep = "_")
  # RunEnrichment() merges into an existing result of the same key; start clean
  # so the stored table is exactly this run.
  srt@tools[[res_key]] <- NULL
  srt <- scop::RunEnrichment(srt, group.by = group_by, test.use = "wilcox",
                             DE_threshold = enrich_threshold(lfc, padj),
                             species = species, TERM2GENE = t2g, TERM2NAME = t2n,
                             minGSSize = 10, maxGSSize = 500, verbose = FALSE)
  out <- srt@tools[[res_key]]
  if (is.data.frame(out$enrichment) && nrow(out$enrichment)) {
    out$enrichment$Database[out$enrichment$Database == "custom"] <- db
  }
  out$universe_size <- length(intersect(universe, t2g$symbol))
  out$universe_tested <- length(universe)
  srt@tools[[res_key]] <- out
  srt
}

#' DE_threshold expression for the ORA gene lists
#' @param lfc,padj See [sc_enrichment()].
#' @keywords internal
enrich_threshold <- function(lfc = 0.25, padj = 0.05) {
  sprintf("avg_log2FC > %s & p_val_adj < %s", format(lfc), format(padj))
}

#' GSEA per group on every tested gene ranked by avg log2 fold-change
#' @inheritParams sc_enrichment
#' @return The object; result in `srt@tools[["GSEA_<group_by>_wilcox"]]`.
#' @keywords internal
sc_gsea <- function(srt, group_by, db = "GO_BP", species = "Homo_sapiens",
                    reuse_de = FALSE) {
  srt <- sc_ensure_detest(srt, group_by, reuse = reuse_de)
  key <- enrich_db_key(db, species)
  res_key <- paste("GSEA", group_by, "wilcox", sep = "_")
  srt@tools[[res_key]] <- NULL
  srt <- scop::RunGSEA(srt, group.by = group_by, test.use = "wilcox",
                       DE_threshold = "p_val <= 1", species = species, db = key,
                       minGSSize = 10, maxGSSize = 500, verbose = FALSE)
  out <- srt@tools[[res_key]]
  if (is.data.frame(out$enrichment) && nrow(out$enrichment) && !identical(key, db)) {
    out$enrichment$Database[out$enrichment$Database == key] <- db
    srt@tools[[res_key]] <- out
  }
  srt
}

#' The stored enrichment / GSEA table of the last run
#' @param srt Seurat object. @param analysis "ora" or "gsea".
#' @param group_by Grouping column.
#' @return data.frame (possibly empty).
#' @keywords internal
enrichment_table <- function(srt, analysis, group_by) {
  key <- paste(if (identical(analysis, "ora")) "Enrichment" else "GSEA",
               group_by, "wilcox", sep = "_")
  tab <- tryCatch(srt@tools[[key]][["enrichment"]], error = function(e) NULL)
  if (is.data.frame(tab)) tab else data.frame()
}

#' Headline numbers of an enrichment / GSEA table
#' @param tab Output of [enrichment_table()]. @param q BH cut-off.
#' @return list(n_tested, n_sig, groups, groups_hit).
#' @keywords internal
enrichment_summary <- function(tab, q = 0.05) {
  if (!is.data.frame(tab) || !nrow(tab) || !"p.adjust" %in% colnames(tab)) {
    return(list(n_tested = 0L, n_sig = 0L, groups = 0L, groups_hit = 0L))
  }
  sig <- !is.na(tab$p.adjust) & tab$p.adjust < q
  list(n_tested = nrow(tab),
       n_sig = sum(sig),
       groups = length(unique(tab$Groups)),
       groups_hit = length(unique(tab$Groups[sig])))
}

#' R code reproducing [sc_enrichment()] / [sc_gsea()]
#' @param analysis "ora" or "gsea". @inheritParams sc_enrichment
#' @keywords internal
enrichment_log_code <- function(analysis, group_by, db, species, lfc = 0.25,
                                padj = 0.05) {
  key <- enrich_db_key(db, species)
  de <- detest_log_code(group_by)
  if (identical(analysis, "gsea")) {
    res_key <- paste("GSEA", group_by, "wilcox", sep = "_")
    relabel <- if (!identical(key, db)) {
      c(sprintf("enr <- obj@tools[[%s]]$enrichment", r_lit(res_key)),
        sprintf("enr$Database[enr$Database == %s] <- %s", r_lit(key), r_lit(db)),
        sprintf("obj@tools[[%s]]$enrichment <- enr", r_lit(res_key)))
    }
    return(c(
      de,
      sprintf("obj@tools[[%s]] <- NULL", r_lit(res_key)),
      sprintf(paste0("obj <- scop::RunGSEA(obj, group.by = %s, test.use = \"wilcox\", ",
                     "DE_threshold = \"p_val <= 1\", species = %s, db = %s, ",
                     "minGSSize = 10, maxGSSize = 500)"),
              r_lit(group_by), r_lit(species), r_lit(key)),
      relabel))
  }
  res_key <- paste("Enrichment", group_by, "wilcox", sep = "_")
  c(
    de,
    "# Background = genes tested in the DE step (not the whole annotation database)",
    sprintf("universe <- unique(as.character(obj@tools[[%s]][[\"AllMarkers_wilcox\"]]$gene))",
            r_lit(paste0("DEtest_", group_by))),
    sprintf("db_list <- scop::PrepareDB(species = %s, db = %s, db_IDtypes = \"symbol\")",
            r_lit(species), r_lit(key)),
    sprintf("term2gene <- db_list[[%s]][[%s]][[\"TERM2GENE\"]][, c(\"Term\", \"symbol\")]",
            r_lit(species), r_lit(key)),
    "term2gene <- unique(term2gene[term2gene$symbol %in% universe, , drop = FALSE])",
    sprintf("term2name <- db_list[[%s]][[%s]][[\"TERM2NAME\"]]", r_lit(species), r_lit(key)),
    sprintf("obj@tools[[%s]] <- NULL", r_lit(res_key)),
    sprintf(paste0("obj <- scop::RunEnrichment(obj, group.by = %s, test.use = \"wilcox\", ",
                   "DE_threshold = %s, species = %s, TERM2GENE = term2gene, ",
                   "TERM2NAME = term2name, minGSSize = 10, maxGSSize = 500)"),
            r_lit(group_by), r_lit(enrich_threshold(lfc, padj)), r_lit(species)),
    sprintf("enr <- obj@tools[[%s]]$enrichment", r_lit(res_key)),
    sprintf("enr$Database[enr$Database == \"custom\"] <- %s", r_lit(db)),
    sprintf("obj@tools[[%s]]$enrichment <- enr", r_lit(res_key)))
}

#' Preview of an enrichment / GSEA run (scop dot-matrix "comparison" plot)
#'
#' Rows are the top gene sets per group (BH q < 0.05), columns the groups.
#' @param srt Seurat object. @param analysis "ora" or "gsea".
#' @param group_by Grouping column. @param db Database label stored in the table.
#' @keywords internal
sc_enrichplot <- function(srt, analysis, group_by, db) {
  if (!require_pkgs("scop", "Enrichment plot")) return(NULL)
  s <- enrichment_summary(enrichment_table(srt, analysis, group_by))
  if (s$n_sig == 0) {
    return(empty_plot(paste0("No gene set reached BH q < 0.05 in any group.\n",
                             "Nothing to draw; the full table is in obj@tools.")))
  }
  if (identical(analysis, "ora")) {
    scop::EnrichmentPlot(srt, db = db, group.by = group_by, test.use = "wilcox",
                         plot_type = "comparison", padjustCutoff = 0.05, topTerm = 5)
  } else {
    scop::GSEAPlot(srt, db = db, group.by = group_by, test.use = "wilcox",
                   plot_type = "comparison", padjustCutoff = 0.05, topTerm = 5)
  }
}

#' A blank ggplot carrying a message (used when there is nothing to draw)
#' @param msg Message text.
#' @keywords internal
empty_plot <- function(msg) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = msg, size = 5) +
    ggplot2::theme_void()
}

# ---- Trajectory --------------------------------------------------------------

#' Metadata columns a trajectory method writes (regular expression)
#' @param method Trajectory method key.
#' @keywords internal
traj_output_pattern <- function(method) {
  switch(method,
         slingshot = "^Slingshot_",
         monocle2  = "^Monocle2_",
         monocle3  = "^Monocle3_",
         paga      = "^dpt_pseudotime$",
         palantir  = "^palantir_",
         wot       = "^(trajectory|fates|transition|coupling)_",
         stop("Unknown trajectory method: ", method, call. = FALSE))
}

#' Pseudotime columns among the columns a run added
#'
#' WOT returns transport-map trajectory / fate scores, not a pseudotime, so it
#' yields none here.
#' @param new_cols Columns added by the run. @param method Method key.
#' @return Character vector (possibly empty).
#' @keywords internal
traj_pt_cols <- function(new_cols, method) {
  pat <- switch(method,
                slingshot = "^Slingshot_Lineage[0-9]+$",
                monocle2  = "^Monocle2_Pseudotime$",
                monocle3  = "^Monocle3_Pseudotime$",
                paga      = "^dpt_pseudotime$",
                palantir  = "^palantir_pseudotime$",
                wot       = "^$",
                stop("Unknown trajectory method: ", method, call. = FALSE))
  new_cols[grepl(pat, new_cols)]
}

#' Check the trajectory inputs before a run
#'
#' Monocle3 without a root opens an interactive gadget; PAGA and Palantir need
#' a root group; WOT needs real sampling times.
#' @param method Method key. @param start Start group ("" or NULL = none).
#' @param levels Levels of the grouping column.
#' @param reductions Reduction names of the object.
#' @param time_col,time_values For WOT: the sampling-time column and its values.
#' @return NULL when valid, otherwise list(en =, zh =) with the reason.
#' @keywords internal
traj_check_start <- function(method, start, levels, reductions = "umap",
                             time_col = NULL, time_values = NULL) {
  has_start <- length(start) == 1 && !is.na(start) && nzchar(start)
  if (has_start && !start %in% levels) {
    return(list(en = sprintf("Start group '%s' is not a level of the chosen grouping.", start),
                zh = sprintf("起始分组“%s”不是所选分组列中的取值。", start)))
  }
  if (method %in% c("monocle3", "paga", "palantir") && !has_start) {
    return(list(en = sprintf("%s needs a start group to root the pseudotime.", method),
                zh = sprintf("%s 需要指定起始分组作为拟时序的根。", method)))
  }
  if (method %in% c("monocle3", "paga", "palantir") && !"umap" %in% reductions) {
    return(list(en = sprintf("%s needs a UMAP; run the Embed step first.", method),
                zh = sprintf("%s 需要 UMAP；请先运行降维图步骤。", method)))
  }
  if (identical(method, "wot")) {
    tv <- if (is.numeric(time_values)) time_values else numeric(0)
    if (is.null(time_col) || !nzchar(time_col) || length(unique(tv[!is.na(tv)])) < 2) {
      return(list(en = "WOT needs a numeric sampling-time column with at least two time points.",
                  zh = "WOT 需要一个数值型的采样时间列，且至少包含两个时间点。"))
    }
  }
  NULL
}

#' Linear space trajectories are fitted in
#'
#' The Integrate step records the corrected space; a 2-D UMAP distorts
#' distances and is never used for fitting.
#' @param obj A Seurat object.
#' @keywords internal
traj_default_reduction <- function(obj) {
  red <- NULL
  if (exists("default_graph_reduction", mode = "function")) red <- default_graph_reduction(obj)
  red %||% "pca"
}

#' Linear spaces a trajectory can be fitted in
#' @param obj A Seurat object.
#' @return Reduction names other than the 2-D display maps.
#' @keywords internal
traj_reduction_choices <- function(obj) {
  if (exists("graph_reduction_choices", mode = "function")) return(graph_reduction_choices(obj))
  reds <- obj_reductions(obj)
  reds[!grepl("umap|tsne|pacmap", tolower(reds))]
}

#' Numeric metadata columns that look like sampling time points (for WOT)
#' @param md Cell metadata.
#' @return Names of numeric columns with 2-50 distinct values.
#' @keywords internal
traj_time_cols <- function(md) {
  if (!ncol(md)) return(character(0))
  ok <- vapply(md, function(v) {
    is.numeric(v) && length(unique(v[!is.na(v)])) %in% 2:50
  }, logical(1))
  names(md)[ok]
}

#' Trajectory / pseudotime with one of scop's engines
#'
#' Columns the chosen method wrote in an earlier run are removed first, so the
#' pseudotime found afterwards is from this run only.
#' @param srt Seurat object. @param method slingshot / monocle2 / monocle3 /
#'   paga / palantir / wot.
#' @param group_by Grouping column. @param start Start group or NULL.
#' @param reduction Linear space for Slingshot / PAGA / Palantir.
#' @param dims Number of `reduction` dimensions Slingshot uses.
#' @param time_col WOT sampling-time column.
#' @return list(obj, pt_cols, fate_cols, new_cols, nd = Slingshot dims used).
#' @keywords internal
sc_trajectory <- function(srt, method = "slingshot", group_by, start = NULL,
                          reduction = "pca", dims = 20, time_col = NULL) {
  start <- if (length(start) == 1 && nzchar(start)) start else NULL
  stale <- grep(traj_output_pattern(method), colnames(srt@meta.data), value = TRUE)
  for (col in stale) srt[[col]] <- NULL
  before <- colnames(srt@meta.data)
  md <- srt@meta.data
  nd <- NA_integer_
  if (method == "slingshot") {
    nd <- min(dims, ncol(SeuratObject::Embeddings(srt, reduction = reduction)))
    srt <- scop::RunSlingshot(srt, group.by = group_by, reduction = reduction,
                              dims = seq_len(nd), start = start, prefix = "Slingshot",
                              show_plot = FALSE, verbose = FALSE)
  } else if (method == "monocle2") {
    srt <- scop::RunMonocle2(srt, group.by = group_by, show_plot = FALSE, verbose = FALSE)
    if (!is.null(start)) srt <- monocle2_reroot(srt, md[[group_by]] %in% start)
  } else if (method == "monocle3") {
    root <- rownames(md)[as.character(md[[group_by]]) %in% start]
    srt <- scop::RunMonocle3(srt, group.by = group_by, reduction = "umap",
                             root_cells = root, show_plot = FALSE, verbose = FALSE)
  } else if (method == "paga") {
    srt <- scop::RunPAGA(srt, group.by = group_by, linear_reduction = reduction,
                         nonlinear_reduction = "umap", infer_pseudotime = TRUE,
                         root_group = start, show_plot = FALSE, save_plot = FALSE,
                         dirpath = scop_outdir("paga"), return_seurat = TRUE,
                         verbose = FALSE)
  } else if (method == "palantir") {
    srt <- scop::RunPalantir(srt, group.by = group_by, linear_reduction = reduction,
                             nonlinear_reduction = "umap", early_group = start,
                             show_plot = FALSE, save_plot = FALSE,
                             dirpath = scop_outdir("palantir"), return_seurat = TRUE,
                             verbose = FALSE)
  } else if (method == "wot") {
    tv <- md[[time_col]]
    srt <- scop::RunWOT(srt, group.by = group_by, time_field = time_col,
                        time_from = min(tv, na.rm = TRUE), time_to = max(tv, na.rm = TRUE),
                        tmap_out = file.path(scop_outdir("wot"), "tmap"),
                        show_plot = FALSE, save_plot = FALSE,
                        dirpath = scop_outdir("wot"), return_seurat = TRUE,
                        verbose = FALSE)
  } else {
    stop("Unknown trajectory method: ", method, call. = FALSE)
  }
  new_cols <- setdiff(colnames(srt@meta.data), before)
  pt_cols <- traj_pt_cols(new_cols, method)
  fate_cols <- if (method == "wot") grep("^trajectory_", new_cols, value = TRUE) else character(0)
  srt@misc$omicone_trajectory <- list(method = method, group_by = group_by,
                                      start = start, pt_cols = pt_cols)
  list(obj = srt, pt_cols = pt_cols, fate_cols = fate_cols, new_cols = new_cols, nd = nd)
}

#' Re-root a Monocle2 ordering at the State holding most start-group cells
#'
#' RunMonocle2() roots at the first State; its `root_state` is a Monocle State,
#' not a cell group, so the State is picked from the cells the user named.
#' @param srt Seurat object after scop::RunMonocle2(). @param in_start Logical
#'   vector over cells: is the cell in the start group?
#' @keywords internal
monocle2_reroot <- function(srt, in_start) {
  st <- as.character(srt$Monocle2_State)
  root <- names(sort(table(st[in_start]), decreasing = TRUE))[1]
  cds <- monocle::orderCells(srt@tools$Monocle2$cds, root_state = root)
  srt$Monocle2_Pseudotime <- cds$Pseudotime
  srt@tools$Monocle2$cds <- cds
  srt
}

#' R code reproducing [sc_trajectory()]
#' @inheritParams sc_trajectory
#' @param nd Slingshot dimensions actually used.
#' @keywords internal
trajectory_log_code <- function(method, group_by, start = NULL, reduction = "pca",
                                nd = 20, time_col = NULL) {
  has_start <- length(start) == 1 && nzchar(start)
  g <- r_lit(group_by)
  drop <- sprintf("for (col in grep(%s, colnames(obj@meta.data), value = TRUE)) obj[[col]] <- NULL",
                  r_lit(traj_output_pattern(method)))
  run <- switch(
    method,
    slingshot = sprintf(paste0("obj <- scop::RunSlingshot(obj, group.by = %s, reduction = %s, ",
                               "dims = 1:%d, start = %s, prefix = \"Slingshot\", ",
                               "show_plot = FALSE)"),
                        g, r_lit(reduction), nd, if (has_start) r_lit(start) else "NULL"),
    monocle2 = c(
      sprintf("obj <- scop::RunMonocle2(obj, group.by = %s, show_plot = FALSE)", g),
      if (has_start) c(
        sprintf("in_start <- obj@meta.data[[%s]] %%in%% %s", g, r_lit(start)),
        "root <- names(sort(table(obj$Monocle2_State[in_start]), decreasing = TRUE))[1]",
        "cds <- monocle::orderCells(obj@tools$Monocle2$cds, root_state = root)",
        "obj$Monocle2_Pseudotime <- cds$Pseudotime",
        "obj@tools$Monocle2$cds <- cds")),
    monocle3 = c(
      sprintf("root <- colnames(obj)[as.character(obj@meta.data[[%s]]) %%in%% %s]",
              g, r_lit(start)),
      sprintf(paste0("obj <- scop::RunMonocle3(obj, group.by = %s, reduction = \"umap\", ",
                     "root_cells = root, show_plot = FALSE)"), g)),
    paga = sprintf(paste0("obj <- scop::RunPAGA(obj, group.by = %s, linear_reduction = %s, ",
                          "nonlinear_reduction = \"umap\", infer_pseudotime = TRUE, ",
                          "root_group = %s, show_plot = FALSE, save_plot = FALSE, ",
                          "dirpath = file.path(tempdir(), \"paga\"), return_seurat = TRUE)"),
                   g, r_lit(reduction), r_lit(start)),
    palantir = sprintf(paste0("obj <- scop::RunPalantir(obj, group.by = %s, ",
                              "linear_reduction = %s, ",
                              "nonlinear_reduction = \"umap\", early_group = %s, ",
                              "show_plot = FALSE, save_plot = FALSE, ",
                              "dirpath = file.path(tempdir(), \"palantir\"), ",
                              "return_seurat = TRUE)"),
                       g, r_lit(reduction), r_lit(start)),
    wot = c(
      sprintf("tv <- obj@meta.data[[%s]]", r_lit(time_col)),
      sprintf(paste0("obj <- scop::RunWOT(obj, group.by = %s, time_field = %s, ",
                     "time_from = min(tv, na.rm = TRUE), time_to = max(tv, na.rm = TRUE), ",
                     "tmap_out = file.path(tempdir(), \"wot\", \"tmap\"), show_plot = FALSE, ",
                     "save_plot = FALSE, dirpath = file.path(tempdir(), \"wot\"), ",
                     "return_seurat = TRUE)"),
              g, r_lit(time_col))),
    stop("Unknown trajectory method: ", method, call. = FALSE))
  c(drop, run)
}

#' The 2-D map trajectories are drawn on
#' @param srt Seurat object.
#' @return "umap" when present, else scop's default reduction (NULL lets scop pick).
#' @keywords internal
traj_display_reduction <- function(srt) {
  if (has_reduction(srt, "umap")) "umap" else NULL
}

#' Preview of a trajectory run
#'
#' Slingshot: cells coloured by group with the lineage curves (pseudotime fitted
#' in the linear space, curves smoothed onto the map). Monocle: pseudotime
#' colour with the principal graph. PAGA / Palantir: pseudotime colour, no
#' curve. WOT: the first transport-map trajectory score.
#' @param srt Seurat object. @param method Method key. @param group_by Grouping.
#' @param pt_cols,fate_cols From [sc_trajectory()].
#' @keywords internal
sc_trajectory_plot <- function(srt, method, group_by, pt_cols, fate_cols = character(0)) {
  if (!require_pkgs("scop", "Trajectory plot")) return(NULL)
  red <- traj_display_reduction(srt)
  if (method == "slingshot") {
    pp <- scop_group_prep(srt, group_by)
    return(scop::CellDimPlot(pp$srt, group.by = group_by, reduction = red, palcolor = pp$palcolor,
                             lineages = pt_cols, lineages_span = 0.1))
  }
  if (method %in% c("monocle2", "monocle3")) {
    tool <- srt@tools[[if (method == "monocle2") "Monocle2" else "Monocle3"]]
    red_m <- if (method == "monocle2") "DDRTree" else red
    p <- scop::FeatureDimPlot(srt, features = pt_cols[1], reduction = red_m)
    if (!is.null(tool$trajectory)) p <- p + tool$trajectory
    return(p)
  }
  feat <- if (length(pt_cols)) pt_cols[1] else fate_cols[1]
  scop::FeatureDimPlot(srt, features = feat, reduction = red)
}

# ---- RNA velocity ------------------------------------------------------------

#' Missing prerequisites for RNA velocity
#' @param assays Assay names of the object. @param reductions Reduction names.
#' @return Character vector of what is missing (empty = ready).
#' @keywords internal
velocity_missing <- function(assays, reductions) {
  miss <- setdiff(c("spliced", "unspliced"), assays)
  miss <- c(miss, setdiff(c("pca", "umap"), reductions))
  miss
}

#' RNA velocity with scVelo through scop::RunSCVELO()
#' @param srt Seurat object with spliced / unspliced assays.
#' @param group_by Grouping column. @param mode stochastic / deterministic /
#'   dynamical.
#' @keywords internal
sc_velocity <- function(srt, group_by, mode = "stochastic") {
  miss <- velocity_missing(obj_assays(srt), obj_reductions(srt))
  if (length(miss)) {
    stop("RNA velocity needs: ", paste(miss, collapse = ", "),
         " (spliced/unspliced assays from velocyto or kallisto|bustools; pca and umap reductions).",
         call. = FALSE)
  }
  scop::RunSCVELO(srt, group.by = group_by, linear_reduction = "pca",
                  nonlinear_reduction = "umap", mode = mode, show_plot = FALSE,
                  save_plot = FALSE, dirpath = scop_outdir("scvelo"),
                  return_seurat = TRUE, verbose = FALSE)
}

#' R code reproducing [sc_velocity()]
#' @inheritParams sc_velocity
#' @keywords internal
velocity_log_code <- function(group_by, mode) {
  sprintf(paste0("obj <- scop::RunSCVELO(obj, group.by = %s, linear_reduction = \"pca\", ",
                 "nonlinear_reduction = \"umap\", mode = %s, show_plot = FALSE, ",
                 "save_plot = FALSE, dirpath = file.path(tempdir(), \"scvelo\"), ",
                 "return_seurat = TRUE)"),
          r_lit(group_by), r_lit(mode))
}

#' Velocity streamlines on the UMAP over cells coloured by group
#'
#' scop::VelocityPlot() colours only its "raw" arrows by group; CellDimPlot()
#' draws the cells by group and the stream layer on top.
#' @param srt Seurat object after [sc_velocity()]. @param mode Velocity mode.
#' @param group_by Grouping column.
#' @keywords internal
sc_velocityplot <- function(srt, mode, group_by) {
  if (!require_pkgs("scop", "Velocity plot")) return(NULL)
  pp <- scop_group_prep(srt, group_by)
  scop::CellDimPlot(pp$srt, group.by = group_by, reduction = "umap", velocity = mode,
                    velocity_plot_type = "stream", palcolor = pp$palcolor)
}

# ---- Dynamic features ----------------------------------------------------------

#' Pseudotime columns usable as lineages for dynamic features
#'
#' The Trajectory step records what it produced in
#' `obj@misc$omicone_trajectory$pt_cols`; otherwise numeric columns with the
#' names scop's trajectory engines write are taken.
#' @param md Cell metadata. @param recorded Recorded pseudotime columns or NULL.
#' @return Character vector (possibly empty).
#' @keywords internal
dynamic_lineage_cols <- function(md, recorded = NULL) {
  num <- names(md)[vapply(md, is.numeric, logical(1))]
  rec <- intersect(recorded, num)
  if (length(rec)) return(rec)
  pat <- paste0("^(.+_)?Lineage[0-9]+$|^Monocle[23]_Pseudotime$|",
                "^dpt_pseudotime$|^palantir_pseudotime$")
  num[grepl(pat, num)]
}

#' Genes varying along pseudotime, scop::RunDynamicFeatures()
#' @param srt Seurat object. @param lineages Pseudotime column(s).
#' @param n_candidates Candidate genes (highly variable) to test.
#' @keywords internal
sc_dynamic <- function(srt, lineages, n_candidates = 1000) {
  if (!length(lineages)) {
    stop("No pseudotime column found. Run the Trajectory step first.", call. = FALSE)
  }
  scop::RunDynamicFeatures(srt, lineages = lineages, n_candidates = n_candidates,
                           verbose = FALSE)
}

#' R code reproducing [sc_dynamic()]
#' @inheritParams sc_dynamic
#' @keywords internal
dynamic_log_code <- function(lineages, n_candidates) {
  sprintf("obj <- scop::RunDynamicFeatures(obj, lineages = %s, n_candidates = %d)",
          r_lit(lineages), as.integer(n_candidates))
}

#' Number of dynamic genes per lineage at DynamicHeatmap()'s default filters
#'
#' r.sq > 0.2, deviance explained > 0.2 and BH-type adjusted p < 0.05, the
#' thresholds scop::DynamicHeatmap() applies by default.
#' @param srt Seurat object. @param lineages Lineage columns.
#' @return Named integer vector.
#' @keywords internal
dynamic_counts <- function(srt, lineages) {
  vapply(lineages, function(l) {
    df <- tryCatch(srt@tools[[paste0("DynamicFeatures_", l)]][["DynamicFeatures"]],
                   error = function(e) NULL)
    if (!is.data.frame(df) || !nrow(df)) return(0L)
    sum(df$r.sq > 0.2 & df$dev.expl > 0.2 & df$padjust < 0.05, na.rm = TRUE)
  }, integer(1))
}

#' Dynamic (pseudotime) heatmap, scop::DynamicHeatmap()
#' @param srt Seurat object. @param lineages Lineage columns.
#' @keywords internal
sc_dynamicheatmap <- function(srt, lineages) {
  if (!require_pkgs("scop", "DynamicHeatmap")) return(NULL)
  ht <- scop::DynamicHeatmap(srt, lineages = lineages, r.sq = 0.2, dev.expl = 0.2,
                             padjust = 0.05, verbose = FALSE)
  ht$plot %||% ht
}

# ---- Cell cycle & signatures ---------------------------------------------------

#' Seurat's 2019 S / G2M gene lists, in the object's species casing
#' @param species "human" or "mouse".
#' @return list(S = , G2M = ).
#' @keywords internal
cellcycle_genes <- function(species = "human") {
  cc <- Seurat::cc.genes.updated.2019
  out <- list(S = cc$s.genes, G2M = cc$g2m.genes)
  if (identical(species, "mouse")) out <- lapply(out, mouse_case)
  out
}

#' Cell-cycle scoring, Seurat::CellCycleScoring()
#' @param srt Seurat object (normalised). @param species "human" or "mouse".
#' @return list(obj, coverage).
#' @keywords internal
sc_cellcycle <- function(srt, species = guess_species(srt)) {
  cov <- geneset_coverage(cellcycle_genes(species), rownames(srt))
  geneset_require(cov, min_genes = 3)
  obj <- Seurat::CellCycleScoring(srt, s.features = cov$found[[1]],
                                  g2m.features = cov$found[[2]], set.ident = FALSE)
  list(obj = obj, coverage = cov)
}

#' R code reproducing [sc_cellcycle()]
#' @param species "human" or "mouse".
#' @keywords internal
cellcycle_log_code <- function(species = "human") {
  conv <- if (identical(species, "mouse")) {
    c("# mouse data: human symbols converted by case (MKI67 -> Mki67)",
      "cc <- lapply(cc, function(x) paste0(substr(x, 1, 1), tolower(substring(x, 2))))")
  }
  c("cc <- Seurat::cc.genes.updated.2019",
    conv,
    "s_genes <- intersect(cc$s.genes, rownames(obj))",
    "g2m_genes <- intersect(cc$g2m.genes, rownames(obj))",
    paste0("obj <- Seurat::CellCycleScoring(obj, s.features = s_genes, ",
           "g2m.features = g2m_genes, set.ident = FALSE)"))
}

#' Metadata names AddModuleScore() writes, and the names they are given
#'
#' Seurat names the scores `<name><index>`; with the set name as `name` the
#' columns become e.g. "EMT2", which nothing downstream recognises. A fixed
#' prefix is scored and renamed to `<set>_AMS`.
#' @param set_names Gene-set names. @param tmp Temporary prefix.
#' @return Named character vector: new name -> Seurat's name.
#' @keywords internal
ams_col_map <- function(set_names, tmp = "omicone_ams") {
  stats::setNames(paste0(tmp, seq_along(set_names)), paste0(set_names, "_AMS"))
}

#' Score column names a method writes for these sets
#' @param set_names Gene-set names. @param method "UCell" or "AddModuleScore".
#' @keywords internal
modulescore_cols <- function(set_names, method = "UCell") {
  if (identical(method, "UCell")) paste0(set_names, "_UCell") else paste0(set_names, "_AMS")
}

#' Signature scoring with UCell or Seurat::AddModuleScore()
#' @param srt Seurat object. @param features Named list of gene sets (already
#'   restricted to genes in the object). @param method "UCell" / "AddModuleScore".
#' @return The object with one score column per set.
#' @keywords internal
sc_modulescore <- function(srt, features, method = "UCell") {
  if (identical(method, "UCell")) {
    return(UCell::AddModuleScore_UCell(srt, features = features, name = "_UCell"))
  }
  map <- ams_col_map(names(features))
  srt <- Seurat::AddModuleScore(srt, features = unname(features), name = "omicone_ams",
                                seed = 1)
  for (i in seq_along(map)) {
    srt@meta.data[[names(map)[i]]] <- srt@meta.data[[map[[i]]]]
    srt@meta.data[[map[[i]]]] <- NULL
  }
  srt
}

#' R code reproducing [sc_modulescore()]
#' @inheritParams sc_modulescore
#' @keywords internal
modulescore_log_code <- function(features, method = "UCell") {
  sets <- sprintf("sets <- %s", r_lit(features))
  if (identical(method, "UCell")) {
    return(c(sets, "obj <- UCell::AddModuleScore_UCell(obj, features = sets, name = \"_UCell\")"))
  }
  c(sets,
    "obj <- Seurat::AddModuleScore(obj, features = unname(sets), name = \"omicone_ams\", seed = 1)",
    "for (i in seq_along(sets)) {",
    paste0("  obj@meta.data[[paste0(names(sets)[i], \"_AMS\")]] <- ",
           "obj@meta.data[[paste0(\"omicone_ams\", i)]]"),
    "  obj@meta.data[[paste0(\"omicone_ams\", i)]] <- NULL",
    "}")
}

#' User-defined stemness gene-set score (UCell)
#'
#' A rank-based score of the genes the user supplies. It is not mRNAsi
#' (Malta et al. 2018, an OCLR model trained on stem-cell transcriptomes) and
#' not a validated potency estimate such as CytoTRACE 2.
#' @param srt Seurat object. @param genes Gene symbols.
#' @return list(obj, coverage); score column "stemness_UCell".
#' @keywords internal
sc_stemness <- function(srt, genes) {
  cov <- geneset_coverage(list(stemness = genes), rownames(srt))
  geneset_require(cov, min_genes = 3)
  obj <- UCell::AddModuleScore_UCell(srt, features = list(stemness = cov$found[[1]]),
                                     name = "_UCell")
  list(obj = obj, coverage = cov)
}

#' R code reproducing [sc_stemness()]
#' @param genes Genes found in the object.
#' @keywords internal
stemness_log_code <- function(genes) {
  sprintf(paste0("obj <- UCell::AddModuleScore_UCell(obj, features = list(stemness = %s), ",
                 "name = \"_UCell\")"),
          r_lit(genes))
}

# ---- Cell-cell communication ---------------------------------------------------

#' Cell-group labels CellChat accepts
#'
#' CellChat rejects the label "0", which every Seurat clustering produces;
#' all labels get a "C" prefix then, so the groups stay recognisable.
#' @param labels Character vector of group labels.
#' @keywords internal
cellchat_labels <- function(labels) {
  labels <- as.character(labels)
  if (any(labels == "0", na.rm = TRUE)) labels <- paste0("C", labels)
  labels
}

#' CellChat inference, the full standard workflow
#'
#' createCellChat -> DB by species -> subsetData -> identifyOverExpressedGenes
#' -> identifyOverExpressedInteractions -> computeCommunProb (triMean,
#' 100 permutations) -> filterCommunication -> computeCommunProbPathway ->
#' aggregateNet. Cells of all samples are pooled.
#' @param srt Seurat object (log-normalised `data` layer).
#' @param group_by Grouping column. @param species "human" or "mouse".
#' @param min_cells Minimum cells per group kept by filterCommunication().
#' @param do_fast Use presto's Wilcoxon in identifyOverExpressedGenes() (CellChat
#'   2's default, which errors when presto is not installed).
#' @return A CellChat object.
#' @keywords internal
sc_cellchat <- function(srt, group_by, species = "human", min_cells = 10,
                        do_fast = has_pkg("presto")) {
  data <- obj_layer(srt, "data")
  if (is.null(data)) stop("No normalised 'data' layer; run Normalize first.", call. = FALSE)
  meta <- data.frame(group = cellchat_labels(srt@meta.data[[group_by]]),
                     samples = factor("sample1"), row.names = colnames(srt))
  cc <- CellChat::createCellChat(object = data, meta = meta, group.by = "group")
  cc@DB <- if (identical(species, "mouse")) {
    CellChat::CellChatDB.mouse
  } else {
    CellChat::CellChatDB.human
  }
  cc <- CellChat::subsetData(cc)
  cc <- CellChat::identifyOverExpressedGenes(cc, do.fast = do_fast)
  cc <- CellChat::identifyOverExpressedInteractions(cc)
  cc <- CellChat::computeCommunProb(cc)
  cc <- CellChat::filterCommunication(cc, min.cells = min_cells)
  cc <- CellChat::computeCommunProbPathway(cc)
  CellChat::aggregateNet(cc)
}

#' R code reproducing [sc_cellchat()]
#' @inheritParams sc_cellchat
#' @param prefixed Were the labels prefixed by [cellchat_labels()]?
#' @keywords internal
cellchat_log_code <- function(group_by, species = "human", min_cells = 10, prefixed = FALSE,
                              do_fast = FALSE) {
  lab <- sprintf("labels <- as.character(obj@meta.data[[%s]])", r_lit(group_by))
  c(lab,
    if (prefixed) "labels <- paste0(\"C\", labels)  # CellChat rejects the label \"0\"",
    "meta <- data.frame(group = labels, samples = factor(\"sample1\"), row.names = colnames(obj))",
    paste0("cc <- CellChat::createCellChat(object = ",
           "SeuratObject::LayerData(obj, layer = \"data\"), meta = meta, group.by = \"group\")"),
    sprintf("cc@DB <- CellChat::%s",
            if (identical(species, "mouse")) "CellChatDB.mouse" else "CellChatDB.human"),
    "cc <- CellChat::subsetData(cc)",
    sprintf("cc <- CellChat::identifyOverExpressedGenes(cc, do.fast = %s)", r_lit(isTRUE(do_fast))),
    "cc <- CellChat::identifyOverExpressedInteractions(cc)",
    "cc <- CellChat::computeCommunProb(cc)",
    sprintf("cc <- CellChat::filterCommunication(cc, min.cells = %d)", as.integer(min_cells)),
    "cc <- CellChat::computeCommunProbPathway(cc)",
    "cc <- CellChat::aggregateNet(cc)",
    "lr <- CellChat::subsetCommunication(cc, thresh = 0.05)")
}

#' LIANA consensus over its default methods, with liana_aggregate()
#' @param srt Seurat object. @param group_by Grouping column.
#' @param species "human" or "mouse" (mouse uses the "MouseConsensus" resource).
#' @return The aggregated LIANA table (one row per source-target-LR).
#' @keywords internal
sc_liana <- function(srt, group_by, species = "human") {
  res <- liana::liana_wrap(srt, idents_col = group_by, resource = liana_resource(species))
  liana::liana_aggregate(res)
}

#' LIANA resource for a species
#' @param species "human" or "mouse".
#' @keywords internal
liana_resource <- function(species = "human") {
  if (identical(species, "mouse")) "MouseConsensus" else "Consensus"
}

#' R code reproducing [sc_liana()]
#' @inheritParams sc_liana
#' @keywords internal
liana_log_code <- function(group_by, species = "human") {
  c(sprintf("liana_res <- liana::liana_wrap(obj, idents_col = %s, resource = %s)",
            r_lit(group_by), r_lit(liana_resource(species))),
    "liana_agg <- liana::liana_aggregate(liana_res)",
    "lr <- liana_agg[liana_agg$aggregate_rank < 0.05, ]")
}

#' The significant ligand-receptor pairs of a communication result
#'
#' LIANA: aggregate_rank < 0.05 (a rank-aggregation p-value, not FDR-adjusted).
#' CellChat: permutation p < 0.05 from subsetCommunication().
#' @param result LIANA aggregate table or CellChat object. @param method Method.
#' @return data.frame with columns source, target, ligand, receptor.
#' @keywords internal
cellcomm_sig_table <- function(result, method) {
  if (identical(method, "liana")) {
    df <- as.data.frame(result)
    df <- df[!is.na(df$aggregate_rank) & df$aggregate_rank < 0.05, , drop = FALSE]
    return(data.frame(source = as.character(df$source), target = as.character(df$target),
                      ligand = as.character(df$ligand.complex),
                      receptor = as.character(df$receptor.complex),
                      stringsAsFactors = FALSE))
  }
  df <- CellChat::subsetCommunication(result, thresh = 0.05)
  data.frame(source = as.character(df$source), target = as.character(df$target),
             ligand = as.character(df$ligand), receptor = as.character(df$receptor),
             stringsAsFactors = FALSE)
}

#' Sender x receiver heatmap of the number of significant LR pairs
#' @param tab Output of [cellcomm_sig_table()].
#' @keywords internal
cellcomm_count_plot <- function(tab) {
  if (!nrow(tab)) return(empty_plot("No significant ligand-receptor pair."))
  counts <- as.data.frame(table(source = tab$source, target = tab$target),
                          stringsAsFactors = FALSE)
  counts$Freq <- as.numeric(counts$Freq)
  ggplot2::ggplot(counts, ggplot2::aes(x = .data$target, y = .data$source, fill = .data$Freq)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = .data$Freq), size = 3) +
    ggplot2::scale_fill_gradient(low = "#eef3f8", high = style_tokens()$kept,
                                 name = "Significant\nLR pairs") +
    ggplot2::labs(x = "Receiver (target) group", y = "Sender (source) group",
                  title = "Number of significant ligand-receptor pairs") +
    omicone_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

# ---- Malignant cells: copykat --------------------------------------------------

#' Map copykat predictions to a call and a confidence
#'
#' copykat labels cells "aneuploid", "diploid", "c1:diploid:low.conf",
#' "c2:aneuploid:low.conf" or "not.defined" (filtered cells). Exact matching
#' on "aneuploid" sent the low-confidence aneuploid cells to "normal".
#' @param pred Character vector of copykat.pred values.
#' @return data.frame: call ("malignant" / "normal" / NA), confidence
#'   ("high" / "low" / NA).
#' @keywords internal
copykat_label_map <- function(pred) {
  pred <- as.character(pred)
  aneu <- grepl("aneuploid", pred, fixed = TRUE)
  dipl <- !aneu & grepl("diploid", pred, fixed = TRUE)
  call <- ifelse(aneu, "malignant", ifelse(dipl, "normal", NA_character_))
  conf <- ifelse(is.na(call), NA_character_,
                 ifelse(grepl("low.conf", pred, fixed = TRUE), "low", "high"))
  data.frame(call = call, confidence = conf, stringsAsFactors = FALSE)
}

#' copykat genome label for a species
#' @param species "human" or "mouse".
#' @keywords internal
copykat_genome <- function(species = "human") {
  if (identical(species, "mouse")) "mm10" else "hg20"
}

#' copykat run per sample
#'
#' One copykat run per sample, so the diploid baseline and clustering never
#' mix samples; each sample's counts are made dense on their own (copykat
#' needs a dense matrix), after dropping genes with no count in that sample,
#' which copykat's own filter would drop anyway. copykat writes its tables
#' and heatmaps into the working directory, so it runs in a temporary one.
#' @param srt Seurat object. @param ref_cells Known-normal cell names or NULL.
#' @param sample_col Sample column or NULL (one run). @param species "human"/"mouse".
#' @param min_cells Samples with fewer cells are skipped.
#' @return list(pred = data.frame(cell, sample, copykat.pred), skipped = names).
#' @keywords internal
sc_copykat <- function(srt, ref_cells = NULL, sample_col = NULL, species = "human",
                       min_cells = 50) {
  counts <- obj_counts(srt)
  samples <- if (is.null(sample_col)) rep("all", ncol(srt)) else
    as.character(srt@meta.data[[sample_col]])
  run_dir <- tempfile("omicone_copykat_")
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  old_wd <- setwd(run_dir)
  on.exit(setwd(old_wd), add = TRUE)
  preds <- list()
  skipped <- character(0)
  for (s in unique(samples[!is.na(samples)])) {
    idx <- which(samples == s)
    if (length(idx) < min_cells) {
      skipped <- c(skipped, s)
      next
    }
    m <- counts[, idx, drop = FALSE]
    m <- as.matrix(m[Matrix::rowSums(m) > 0, , drop = FALSE])
    ref <- intersect(ref_cells, colnames(m))
    ck <- copykat::copykat(rawmat = m, id.type = "S", sam.name = make.names(s),
                           genome = copykat_genome(species),
                           norm.cell.names = if (length(ref) > 1) ref else "",
                           plot.genes = "FALSE", output.seg = "FALSE", n.cores = 1)
    p <- as.data.frame(ck$prediction, stringsAsFactors = FALSE)
    preds[[s]] <- data.frame(cell = as.character(p$cell.names), sample = s,
                             copykat.pred = as.character(p$copykat.pred),
                             stringsAsFactors = FALSE)
  }
  if (!length(preds)) {
    stop(sprintf("No sample has at least %d cells; copykat was not run.", min_cells),
         call. = FALSE)
  }
  pred <- do.call(rbind, preds)
  rownames(pred) <- NULL
  list(pred = pred, skipped = skipped)
}

#' Write copykat calls into the object's metadata
#' @param srt Seurat object. @param pred `pred` from [sc_copykat()].
#' @return The object with `copykat_pred`, `malignant` and `malignant_confidence`.
#' @keywords internal
copykat_add_calls <- function(srt, pred) {
  i <- match(colnames(srt), pred$cell)
  raw <- pred$copykat.pred[i]
  lab <- copykat_label_map(raw)
  srt$copykat_pred <- raw
  srt$malignant <- lab$call
  srt$malignant_confidence <- lab$confidence
  srt
}

#' Share of reference (normal) cells copykat called malignant
#' @param calls Named call vector (names = cells). @param ref_cells Reference cells.
#' @return Fraction, or NA without reference cells.
#' @keywords internal
copykat_ref_flagged <- function(calls, ref_cells) {
  r <- calls[intersect(ref_cells, names(calls))]
  r <- r[!is.na(r)]
  if (!length(r)) return(NA_real_)
  mean(r == "malignant")
}

#' R code reproducing [sc_copykat()] + [copykat_add_calls()]
#' @inheritParams sc_copykat
#' @param ref_col,ref_groups Reference column and its normal groups.
#' @keywords internal
copykat_log_code <- function(ref_col = NULL, ref_groups = NULL, sample_col = NULL,
                             species = "human", min_cells = 50) {
  ref <- if (length(ref_groups) && !is.null(ref_col)) {
    sprintf("ref_cells <- colnames(obj)[as.character(obj@meta.data[[%s]]) %%in%% %s]",
            r_lit(ref_col), r_lit(ref_groups))
  } else {
    "ref_cells <- character(0)"
  }
  smp <- if (is.null(sample_col)) "samples <- rep(\"all\", ncol(obj))" else
    sprintf("samples <- as.character(obj@meta.data[[%s]])", r_lit(sample_col))
  c("counts <- SeuratObject::LayerData(obj, assay = \"RNA\", layer = \"counts\")",
    ref, smp,
    "old_wd <- setwd(tempdir())  # copykat writes tables and heatmaps to the working directory",
    "preds <- list()",
    "for (s in unique(samples[!is.na(samples)])) {",
    "  idx <- which(samples == s)",
    sprintf("  if (length(idx) < %d) next", as.integer(min_cells)),
    "  m <- counts[, idx, drop = FALSE]",
    "  m <- as.matrix(m[Matrix::rowSums(m) > 0, , drop = FALSE])",
    "  ref <- intersect(ref_cells, colnames(m))",
    sprintf(paste0("  ck <- copykat::copykat(rawmat = m, id.type = \"S\", ",
                   "sam.name = make.names(s), ",
                   "genome = %s, norm.cell.names = if (length(ref) > 1) ref else \"\", ",
                   "plot.genes = \"FALSE\", output.seg = \"FALSE\", n.cores = 1)"),
            r_lit(copykat_genome(species))),
    "  preds[[s]] <- as.data.frame(ck$prediction)",
    "}",
    "setwd(old_wd)",
    "pred <- do.call(rbind, preds)",
    "raw <- as.character(pred$copykat.pred)[match(colnames(obj), pred$cell.names)]",
    "obj$copykat_pred <- raw",
    paste0("obj$malignant <- ifelse(grepl(\"aneuploid\", raw), \"malignant\", ",
           "ifelse(grepl(\"diploid\", raw), \"normal\", NA))"),
    paste0("obj$malignant_confidence <- ifelse(is.na(obj$malignant), NA, ",
           "ifelse(grepl(\"low.conf\", raw), \"low\", \"high\"))"))
}

# ---- scop plotting wrappers ----------------------------------------------------

#' A grouping prepared for a scop plot in the app's colours
#'
#' scop assigns `palcolor` to the levels in order, so the column is made a
#' factor in [group_colors()] order (on a local copy) and the colours are
#' passed in that order: a cluster has the same colour in a scop figure as in
#' every other figure.
#' @param srt Seurat object. @param group_by Grouping column.
#' @return list(srt, palcolor).
#' @keywords internal
scop_group_prep <- function(srt, group_by) {
  if (is.null(group_by) || !group_by %in% obj_meta_cols(srt)) return(list(srt = srt, palcolor = NULL))
  cols <- group_colors(srt, group_by)
  srt@meta.data[[group_by]] <- factor(as.character(srt@meta.data[[group_by]]), levels = names(cols))
  list(srt = srt, palcolor = unname(cols))
}

#' Dimensional-reduction scatter (clusters / metadata), scop::CellDimPlot(),
#' with optional mascarade outlines of the same grouping
#' @param srt Seurat object. @param group_by Column to colour and outline by.
#' @param reduction Reduction (NULL = scop's default).
#' @param mask Draw mascarade outlines. @param label Passed to scop.
#' @param ... Further CellDimPlot() arguments.
#' @keywords internal
sc_dimplot <- function(srt, group_by, reduction = NULL, mask = FALSE, label = TRUE, ...) {
  if (!require_pkgs("scop", "Dimension plot")) return(NULL)
  pp <- scop_group_prep(srt, group_by)
  p <- scop::CellDimPlot(pp$srt, group.by = group_by, reduction = reduction,
                         palcolor = pp$palcolor, label = label, ...)
  if (isTRUE(mask) && has_pkg("mascarade")) {
    p <- tryCatch(add_mascarade(p, srt, group_by, reduction),
                  error = function(e) p)
  }
  p
}

#' PAGA graph on an embedding, scop::PAGAPlot()
#' @param srt Seurat object. @param ... PAGAPlot() arguments.
#' @keywords internal
sc_pagaplot <- function(srt, ...) {
  if (!require_pkgs("scop", "PAGA plot")) return(NULL)
  scop::PAGAPlot(srt, ...)
}

# ---- mascarade outlines --------------------------------------------------------

#' Outline polygons of each group on a 2-D embedding
#'
#' mascarade::generateMask() names the coordinate columns after the
#' embedding's own column names (UMAP_1 / umap_1 / ...), so they are fixed to
#' x / y first; cells with a missing label are left out.
#' @param emb Two-column embedding matrix (rownames = cells).
#' @param labels Group label per row.
#' @return data.frame with columns x, y, cluster, group.
#' @keywords internal
mascarade_mask <- function(emb, labels) {
  keep <- !is.na(labels)
  dims <- as.matrix(emb[keep, 1:2, drop = FALSE])
  colnames(dims) <- c("x", "y")
  mask <- mascarade::generateMask(dims = dims, clusters = as.character(labels[keep]))
  as.data.frame(mask)
}

#' Mascarade outlines overlaid on a scop dim plot
#' @param p ggplot dim plot. @param srt Seurat object.
#' @param group_by Grouping to outline (the one the user picked).
#' @param reduction Reduction (NULL = the object's default).
#' @keywords internal
add_mascarade <- function(p, srt, group_by, reduction = NULL) {
  if (!require_pkgs("mascarade", "Cell-type outlines")) return(p)
  red <- reduction %||% SeuratObject::DefaultDimReduc(srt)
  emb <- SeuratObject::Embeddings(srt, reduction = red)
  mask <- mascarade_mask(emb, obj_meta(srt)[[group_by]])
  p + ggplot2::geom_path(
    data = mask,
    ggplot2::aes(x = .data$x, y = .data$y, group = .data$group),
    colour = style_tokens()$ink, linewidth = 0.4, inherit.aes = FALSE
  )
}
