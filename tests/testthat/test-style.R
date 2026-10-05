# Figure style (fct_style.R, fct_sc_plots.R): a level keeps its colour, cell
# types inherit their cluster's colour, composed panels line up, and the
# single-cell figures are computed correctly and draw without Seurat plotting.

draw_ok <- function(p) {
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f, width = 9, height = 7)
  on.exit(grDevices::dev.off(), add = TRUE)
  draw_plot_object(p)
  TRUE
}

test_that("levels are ordered naturally, with unassigned-type levels last", {
  expect_equal(level_order(c("10", "2", "1")), c("1", "2", "10"))
  expect_equal(level_order(c("C10", "C2", "Unassigned", "B")), c("B", "C2", "C10", "Unassigned"))
  expect_equal(level_order(factor(c("b", "a"), levels = c("b", "a"))), c("b", "a"))
})

test_that("a level keeps its colour when others are added or removed", {
  a <- extend_colors(NULL, c("T", "B", "NK"))
  b <- extend_colors(a, c("NK", "B", "Mono", "T"))
  expect_equal(unname(b[c("T", "B", "NK")]), unname(a[c("T", "B", "NK")]))
  expect_false(b[["Mono"]] %in% a)                          # a new level gets an unused colour
  expect_equal(unname(extend_colors(a, c("B", "Unassigned"))[["Unassigned"]]), style_tokens()$unassigned)
  v <- value_colors(c("Mutant", "WT"))
  expect_equal(unname(v[c("WT", "Mutant")]), c(style_tokens()$wt, style_tokens()$mutant))
  expect_equal(unname(km_group_colors(c("Low", "High"))), c(style_tokens()$wt, style_tokens()$mutant))
})

test_that("the colour registry is stored in the object, and cell types inherit cluster colours", {
  skip_if_not_installed("Seurat")
  counts <- readRDS(demo_bundled_path())
  obj <- suppressWarnings(as_seurat(counts))
  obj$cl <- factor(rep(c("0", "1", "2"), length.out = ncol(obj)))
  obj <- set_group_colors(obj, "cl")
  cl <- group_colors(obj, "cl")
  map <- c(`0` = "T cells", `1` = "B cells", `2` = "T cells")
  obj$celltype <- unname(map[as.character(obj$cl)])
  obj@misc$omicone_cluster_col <- "cl"
  obj <- color_celltypes(obj)
  ct <- group_colors(obj, "celltype")
  expect_equal(ct[["B cells"]], cl[["1"]])
  expect_true(ct[["T cells"]] %in% cl[c("0", "2")])
  # stored maps win over the default, and survive a subset
  obj2 <- set_group_colors(obj, "cl", c(`0` = "#000000", `1` = "#111111", `2` = "#222222"))
  sub <- subset(obj2, cells = colnames(obj2)[obj2$cl != "1"])
  expect_equal(unname(group_colors(sub, "cl")), c("#000000", "#222222"))
  # fixed meanings for known groupings
  obj$Phase <- rep(c("G1", "S", "G2M"), length.out = ncol(obj))
  expect_equal(group_colors(obj, "Phase")[["G2M"]], style_known_groups()$Phase[["G2M"]])
})

test_that("compose_grid aligns panels in a column and nests grids", {
  p1 <- ggplot2::ggplot(data.frame(x = 1:3, y = c(1, 100000, 3)), ggplot2::aes(x, y)) + ggplot2::geom_col() +
    omicone_theme()
  p2 <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) + ggplot2::geom_point() +
    ggplot2::labs(y = "a much longer axis title") + omicone_theme()
  g <- compose_grid(list(p1, p2), nrow = 2, ncol = 1, title = "t")
  expect_s3_class(g, "omicone_grid")
  cells <- g$grobs[grepl("^cell-", g$layout$name)]
  left <- vapply(cells, function(x) {
    lay <- x$layout[grepl("^panel", x$layout$name), ]
    unit_cm(x$widths[seq_len(min(lay$l) - 1)], "width")
  }, 0)
  expect_equal(left[1], left[2], tolerance = 1e-6)
  nested <- compose_grid(list(g, p1), nrow = 1, ncol = 2)
  expect_true(draw_ok(nested))
})

test_that("the single-cell figures compute and draw (no Seurat plotting)", {
  skip_if_not_installed("Seurat")
  obj <- suppressWarnings(as_seurat(readRDS(demo_bundled_path())))
  obj <- normalize_obj(obj, "LogNormalize")
  obj <- suppressWarnings(reduce_obj(obj, n_hvg = 300, npcs = 15))
  obj <- suppressWarnings(cluster_obj(obj, dims = 10, resolutions = 0.8, algorithm = 1))
  obj <- suppressWarnings(embed_obj(obj, "umap", dims = 10))
  col <- obj_misc(obj, "omicone_cluster_col")
  expect_false(is.null(obj@misc$omicone_colors[[col]]))     # clustering stores its colours
  genes <- rownames(obj)[1:4]
  expect_true(draw_ok(sc_dim_plot(obj, "umap", col)))
  expect_true(draw_ok(sc_feature_plot(obj, genes, "umap")))
  expect_true(draw_ok(sc_violin_plot(obj, genes[1:2], col)))
  expect_true(draw_ok(sc_dot_plot(obj, genes, col)))
  expect_true(draw_ok(sc_heatmap(obj, genes, col)))
  # dot-plot numbers: % expressing and log of the mean of expm1
  d <- sc_dot_data(obj, genes[1], col)
  ex <- as.numeric(obj_layer(obj, "data")[genes[1], ])
  g <- as.character(obj_meta(obj)[[col]])
  lv <- as.character(d$group[1])
  expect_equal(d$pct[1], 100 * mean(ex[g == lv] > 0))
  expect_equal(d$mean[1], log1p(mean(expm1(ex[g == lv]))))
  # a numeric metadata column works as a "feature" (scores, pseudotime)
  obj$score1 <- stats::rnorm(ncol(obj))
  expect_true(draw_ok(sc_feature_plot(obj, "score1")))
  expect_error(sc_feature_plot(obj, "NOT_A_GENE", "umap"), "None of the genes")
})
