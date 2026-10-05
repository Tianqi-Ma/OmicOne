# Sample-level analyses (fct_sample.R): pseudobulk DE and differential
# abundance on a simulated 8-sample experiment with a known answer, the
# sample-sheet join, the guards against pseudo-replication, and that the logged
# code is runnable R that reproduces the app's numbers.

sim_experiment <- function(seed = 1) {
  set.seed(seed)
  samples <- sprintf("S%d", 1:8)
  cond <- rep(c("ctrl", "trt"), each = 4)
  types <- c("Tcell", "Bcell", "Mono")
  n_genes <- 300
  cells <- list()
  for (i in seq_along(samples)) {
    # Bcell share doubles in trt; the others take up the rest
    share <- if (cond[i] == "trt") c(0.45, 0.30, 0.25) else c(0.55, 0.15, 0.30)
    n <- as.vector(stats::rmultinom(1, 240, share * stats::runif(3, 0.85, 1.15)))
    cells[[i]] <- data.frame(sample = samples[i], condition = cond[i], patient = sprintf("P%d", (i - 1) %% 4 + 1),
                             celltype = rep(types, n), stringsAsFactors = FALSE)
  }
  md <- do.call(rbind, cells)
  rownames(md) <- sprintf("cell%04d", seq_len(nrow(md)))
  base <- stats::rgamma(n_genes, 2, 0.5)
  mu <- outer(base, rep(1, nrow(md)))
  up <- md$celltype == "Tcell" & md$condition == "trt"
  mu[1:15, up] <- mu[1:15, up] * 4                                  # genes 1-15 up in T cells
  sample_effect <- stats::setNames(stats::rnorm(8, 0, 0.25), samples)
  mu <- sweep(mu, 2, exp(sample_effect[md$sample]), "*")
  counts <- matrix(stats::rnbinom(length(mu), mu = mu * 0.15, size = 5), n_genes,
                   dimnames = list(sprintf("G%03d", 1:n_genes), rownames(md)))
  list(counts = Matrix::Matrix(counts, sparse = TRUE), md = md)
}

test_that("sample_table: one row per sample, and a cell-level condition is refused", {
  ex <- sim_experiment()
  st <- sample_table(ex$md, "sample", "condition")
  expect_equal(nrow(st), 8)
  expect_equal(sum(st$n_cells), nrow(ex$md))
  expect_setequal(st$condition, c("ctrl", "trt"))
  # a per-cell column cannot be a condition
  expect_error(sample_table(ex$md, "sample", "celltype"), "not constant within 8 sample")
  expect_error(sample_table(ex$md, "sample", "sample"), "different columns")
  expect_error(sample_table(ex$md, "sample", "nope"), "Column not found")
  st2 <- sample_table(ex$md, "sample", "condition", covariate = "patient")
  expect_equal(st2$covariate, rep(sprintf("P%d", 1:4), 2))
})

test_that("apply_sample_sheet joins sample-level columns and reports mismatches", {
  ex <- sim_experiment()
  sheet <- data.frame(id = c(sprintf("S%d", 1:7), "S99"), response = rep(c("R", "NR"), 4),
                      stringsAsFactors = FALSE)
  j <- apply_sample_sheet(ex$md, sheet, "sample")
  expect_equal(nrow(j$cols), nrow(ex$md))
  expect_equal(rownames(j$cols), rownames(ex$md))
  expect_equal(j$missing, "S8")
  expect_equal(j$unused, "S99")
  expect_true(all(is.na(j$cols$response[ex$md$sample == "S8"])))
  expect_true(all(j$cols$response[ex$md$sample == "S1"] == "R"))
  # a sheet column named like the sample column is used as the key
  sheet2 <- data.frame(response = "R", sample = "S1")
  expect_equal(apply_sample_sheet(ex$md, sheet2, "sample")$matched, "S1")
  expect_error(apply_sample_sheet(ex$md, data.frame(id = c("S1", "S1"), x = 1:2), "sample"), "repeated")
})

test_that("pseudobulk_counts sums each sample x cell type and drops small ones", {
  ex <- sim_experiment()
  pb <- pseudobulk_counts(ex$counts, ex$md, "sample", "celltype", min_cells = 0)
  expect_equal(ncol(pb$counts), 24)
  one <- ex$md$sample == "S3" & ex$md$celltype == "Bcell"
  expect_equal(unname(pb$counts[, "S3|Bcell"]), unname(Matrix::rowSums(ex$counts[, one])))
  expect_equal(pb$meta$n_cells[pb$meta$column == "S3|Bcell"], sum(one))
  big <- pseudobulk_counts(ex$counts, ex$md, "sample", "celltype", min_cells = 1000)
  expect_equal(ncol(big$counts), 0)
  expect_equal(nrow(big$dropped), 24)
})

test_that("pseudobulk DE (edgeR, limma) finds the simulated T-cell genes and nothing in B cells", {
  skip_if_not_installed("edgeR")
  skip_if_not_installed("limma")
  ex <- sim_experiment()
  pb <- pseudobulk_counts(ex$counts, ex$md, "sample", "celltype", min_cells = 10)
  st <- sample_table(ex$md, "sample", "condition")
  for (m in c("edgeR", "limma")) {
    t_res <- pseudobulk_de(pb, st, "Tcell", "ctrl", "trt", method = m)
    sig <- t_res$gene[t_res$fdr < 0.05]
    expect_true(all(sprintf("G%03d", 1:15) %in% sig), info = m)
    expect_lte(length(setdiff(sig, sprintf("G%03d", 1:15))), 3)
    expect_true(all(t_res$logFC[t_res$gene %in% sprintf("G%03d", 1:15)] > 1), info = m)
    expect_equal(as.integer(attr(t_res, "n")), c(4L, 4L))
    b_res <- pseudobulk_de(pb, st, "Bcell", "ctrl", "trt", method = m)
    expect_lte(sum(b_res$fdr < 0.05), 2)
  }
  # paired by patient: still finds them
  st2 <- sample_table(ex$md, "sample", "condition", covariate = "patient")
  pr <- pseudobulk_de(pb, st2, "Tcell", "ctrl", "trt", covariate = TRUE)
  expect_true(all(sprintf("G%03d", 1:15) %in% pr$gene[pr$fdr < 0.05]))
})

test_that("pseudobulk DE refuses designs without replication", {
  skip_if_not_installed("edgeR")
  ex <- sim_experiment()
  ex$md$condition[ex$md$sample %in% c("S2", "S3", "S4")] <- "trt"   # 1 ctrl sample left
  pb <- pseudobulk_counts(ex$counts, ex$md, "sample", "celltype")
  st <- sample_table(ex$md, "sample", "condition")
  expect_error(pseudobulk_de(pb, st, "Tcell", "ctrl", "trt"), "at least 2 samples per condition")
  all <- pseudobulk_de_all(pb, st, "ctrl", "trt")
  expect_null(all$table)
  expect_equal(nrow(all$skipped), 3)
  # a covariate that is the condition itself is refused
  ex <- sim_experiment()
  ex$md$batch <- ex$md$condition
  pb <- pseudobulk_counts(ex$counts, ex$md, "sample", "celltype")
  st <- sample_table(ex$md, "sample", "condition", covariate = "batch")
  expect_error(pseudobulk_de(pb, st, "Tcell", "ctrl", "trt", covariate = TRUE), "confounded")
})

test_that("pseudobulk_de_all reports per cell type, FDR within each", {
  skip_if_not_installed("edgeR")
  ex <- sim_experiment()
  pb <- pseudobulk_counts(ex$counts, ex$md, "sample", "celltype", min_cells = 10)
  st <- sample_table(ex$md, "sample", "condition")
  all <- pseudobulk_de_all(pb, st, "ctrl", "trt")
  expect_setequal(unique(all$table$cell_type), c("Bcell", "Mono", "Tcell"))
  ov <- all$overview
  expect_equal(ov$up[ov$cell_type == "Tcell"] >= 15, TRUE)
  expect_true(all(ov$samples_ref == 4 & ov$samples_alt == 4))
  one <- pseudobulk_de(pb, st, "Mono", "ctrl", "trt")
  expect_equal(all$table$fdr[all$table$cell_type == "Mono"], one$fdr)
  expect_s3_class(pseudobulk_volcano(one), "ggplot")
  expect_s3_class(pseudobulk_overview_plot(ov, ref = "ctrl", alt = "trt"), "ggplot")
})

test_that("abundance_test (propeller) finds the B-cell shift and matches a direct limma fit", {
  skip_if_not_installed("limma")
  ex <- sim_experiment()
  res <- abundance_test(ex$md, "sample", "celltype", "condition")
  expect_equal(res$group[1], "Bcell")
  expect_lt(res$fdr[1], 0.05)
  expect_gt(res$prop_ratio[1], 1.5)
  expect_equal(names(res)[2:3], c("prop_ctrl", "prop_trt"))
  # the textbook propeller computation, written out
  tab <- table(ex$md$sample, ex$md$celltype)
  ps <- (tab + 0.5) / rowSums(tab + 0.5)
  y <- t(log(ps / (1 - ps)))
  cond <- factor(c("ctrl", "trt")[(rownames(tab) %in% sprintf("S%d", 5:8)) + 1])
  fit <- limma::eBayes(limma::lmFit(y, stats::model.matrix(~ cond)), robust = TRUE)
  tt <- limma::topTable(fit, coef = 2, number = Inf)
  expect_equal(res$p, tt[res$group, "P.Value"], tolerance = 1e-10)
  expect_s3_class(abundance_plot(res), "ggplot")
  expect_s3_class(abundance_stack_plot(res), "ggplot")
  # three conditions: F test
  ex$md$condition[ex$md$sample %in% c("S7", "S8")] <- "mid"
  ex$md$condition[ex$md$sample %in% c("S3", "S4")] <- "mid2"
  ex$md$condition[ex$md$sample %in% c("S3", "S4")] <- "mid"
  res3 <- abundance_test(ex$md, "sample", "celltype", "condition")
  expect_equal(names(res3)[2:4], c("prop_ctrl", "prop_mid", "prop_trt"))
  expect_false("prop_ratio" %in% names(res3))
  expect_true(all(res3$p >= 0 & res3$p <= 1))
})

test_that("abundance_test refuses fewer than two samples per condition", {
  skip_if_not_installed("limma")
  ex <- sim_experiment()
  ex$md$condition[ex$md$sample != "S1"] <- "trt"
  expect_error(abundance_test(ex$md, "sample", "celltype", "condition"), "at least 2 samples")
})

test_that("the logged pseudobulk and abundance code reproduces the app's results", {
  skip_if_not_installed("edgeR")
  skip_if_not_installed("limma")
  skip_if_not_installed("Seurat")
  ex <- sim_experiment()
  obj <- suppressWarnings(Seurat::CreateSeuratObject(ex$counts, meta.data = ex$md))
  pb <- pseudobulk_counts(ex$counts, ex$md, "sample", "celltype", min_cells = 10)
  st <- sample_table(ex$md, "sample", "condition")
  app <- pseudobulk_de_all(pb, st, "ctrl", "trt", method = "edgeR")
  code <- pseudobulk_log_code("sample", "celltype", "condition", "ctrl", "trt", "edgeR", 10)
  expect_silent(parse(text = code))
  run <- new.env()
  run$obj <- obj
  txt <- paste(code, collapse = "\n")
  if (!exists("LayerData", envir = asNamespace("SeuratObject"), inherits = FALSE)) {
    # Seurat 4: the log targets Seurat 5's LayerData(); replay with the v4 accessor
    txt <- gsub("SeuratObject::LayerData(obj, assay = \"RNA\", layer = \"counts\")",
                "SeuratObject::GetAssayData(obj, assay = \"RNA\", slot = \"counts\")", txt, fixed = TRUE)
  }
  eval(parse(text = txt), run)
  tcell <- run$res$Tcell
  mine <- app$table[app$table$cell_type == "Tcell", ]
  expect_equal(tcell[mine$gene, "PValue"], mine$p, tolerance = 1e-8)
  for (m in c("limma", "DESeq2")) expect_silent(parse(text = pseudobulk_log_code("s", "g", "c", "a", "b", m, 10, covariate = "p")))

  ab <- abundance_test(ex$md, "sample", "celltype", "condition")
  run2 <- new.env()
  run2$obj <- obj
  out <- eval(parse(text = abundance_log_code("sample", "celltype", "condition")), run2)
  expect_equal(out[ab$group, "P.Value"], ab$p, tolerance = 1e-10)
  expect_silent(parse(text = abundance_log_code("s", "g", "c", n_levels = 3)))
})

test_that("the pseudobulk and abundance modules run end to end, with a sample sheet", {
  skip_if_not_installed("edgeR")
  skip_if_not_installed("limma")
  skip_if_not_installed("Seurat")
  ex <- sim_experiment()
  md <- ex$md[, c("sample", "celltype")]           # no condition yet: it comes from the sheet
  obj <- suppressWarnings(Seurat::CreateSeuratObject(ex$counts, meta.data = md))
  rv <- shiny::reactiveValues(obj = obj, source = NULL, status = list(), ckpt = new.env(), clinical = NULL)
  log_rv <- shiny::reactiveVal(list())
  sheet <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(sample = sprintf("S%d", 1:8), arm = rep(c("ctrl", "trt"), each = 4)),
                   sheet, row.names = FALSE)
  bad <- character(0)
  read_outs <- function(output, id, outs) {
    for (nm in outs) {
      err <- tryCatch({
        output[[nm]]
        NULL
      }, error = function(e) e)
      if (!is.null(err) && !inherits(err, "shiny.silent.error")) {
        bad <<- c(bad, sprintf("%s$%s: %s", id, nm, conditionMessage(err)))
      }
    }
  }
  suppressMessages(suppressWarnings(shiny::testServer(
    mod_pseudobulk_server, args = list(rv = rv, log_rv = log_rv), {
      session$flushReact()
      session$setInputs(sample = "sample", group = "celltype", method = "edgeR", min_cells = 10, fdr = 0.05)
      session$setInputs(sheet = list(name = "sheet.csv", datapath = sheet))
      session$flushReact()
      expect_true("arm" %in% obj_meta_cols(rv$obj))
      session$setInputs(condition = "arm", ref = "ctrl", alt = "trt", covariate = "")
      session$setInputs(run = 1)
      session$flushReact()
      out <- res$out
      expect_false(is.null(out))
      expect_setequal(out$overview$cell_type, c("Bcell", "Mono", "Tcell"))
      expect_gte(out$overview$up[out$overview$cell_type == "Tcell"], 15)
      session$setInputs(volcano_ct = "Tcell")
      read_outs(output, "pseudobulk", c("summary", "insight", "preview", "volcano", "volcano_pick", "smp", "tbl",
                                        "sample_ui", "cond_ui", "levels_ui", "cov_ui", "group_ui"))
    })))
  expect_true(isTRUE(shiny::isolate(rv$status$pseudobulk)))
  entry <- Filter(function(e) identical(e$step, "Pseudobulk DE"), shiny::isolate(log_rv()))   # testServer ids are not the step key
  expect_length(entry, 1)
  expect_silent(parse(text = entry[[1]]$code))

  suppressMessages(suppressWarnings(shiny::testServer(
    mod_abundance_server, args = list(rv = rv, log_rv = log_rv), {
      session$flushReact()
      session$setInputs(sample = "sample", condition = "arm", group = "celltype")
      session$setInputs(run = 1)
      session$flushReact()
      expect_equal(res$df$group[1], "Bcell")
      read_outs(output, "abundance", c("summary", "insight", "preview", "stack", "tbl"))
    })))
  expect_true(isTRUE(shiny::isolate(rv$status$abundance)))
  expect_equal(bad, character(0))

  # a cell-level column as the condition is refused with a message, not a crash
  suppressMessages(suppressWarnings(shiny::testServer(
    mod_abundance_server, args = list(rv = rv, log_rv = log_rv), {
      session$setInputs(sample = "sample", condition = "celltype", group = "celltype")
      session$setInputs(run = 1)
      expect_null(res$df)
    })))
})
