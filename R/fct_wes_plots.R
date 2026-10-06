#' WES figures in the shared style (ggplot2 instead of maftools base graphics)
#'
#' maftools computes; these functions draw. Every WES figure uses the app's
#' theme (fct_style.R) and the same colours for the same things: one colour
#' per variant classification ([wes_vc_colors()]), the COSMIC colours for the
#' six base changes ([wes_sbs_colors()]), and [value_colors()] for clinical
#' groups, so a FAB subtype or a cohort keeps its colour from the oncoplot to
#' the survival curve. Multi-panel figures go through [compose_grid()].
#'
#' @name fct_wes_plots
#' @keywords internal
NULL

#' One colour per variant classification (close to the maftools convention)
#' @keywords internal
wes_vc_colors <- function() {
  c(Missense_Mutation = "#3a9a50", Nonsense_Mutation = "#d7362f", Frame_Shift_Del = "#2b6cb0",
    Frame_Shift_Ins = "#7b4fb0", In_Frame_Del = "#e09b3a", In_Frame_Ins = "#c47bb8",
    Splice_Site = "#f08a3c", Splice_Region = "#f6b47f", Translation_Start_Site = "#8c5a3c",
    Nonstop_Mutation = "#4bb3ad", Multi_Hit = "#2d3640", Silent = "#b8c0c8",
    Amp = "#c1476b", Del = "#3b6ea5")
}

# colours for the classes present, in the canonical order; unknown ones grey
wes_vc_scale_values <- function(classes) {
  vc <- wes_vc_colors()
  classes <- unique(as.character(classes))
  known <- intersect(names(vc), classes)
  other <- setdiff(classes, names(vc))
  c(vc[known], stats::setNames(rep(style_tokens()$faint, length(other)), other))
}

#' COSMIC colours of the six base-change classes
#' @keywords internal
wes_sbs_colors <- function() {
  c(`C>A` = "#1EBFF0", `C>G` = "#2d3640", `C>T` = "#E62725", `T>A` = "#CBCACB",
    `T>C` = "#A1CF64", `T>G` = "#EDC8C5")
}

wes_vc_label <- function(x) gsub("_", " ", x)

#' Variant records of a MAF as one data.frame
#' @param maf A MAF object. @param silent Include the silent / non-coding records.
#' @keywords internal
wes_records <- function(maf, silent = FALSE) {
  parts <- list(as.data.frame(maf@data))
  if (isTRUE(silent) && nrow(maf@maf.silent)) parts[[2]] <- as.data.frame(maf@maf.silent)
  cols <- Reduce(union, lapply(parts, names))
  parts <- lapply(parts, function(x) {
    x[setdiff(cols, names(x))] <- NA
    x[cols]
  })
  d <- do.call(rbind, parts)
  d$Tumor_Sample_Barcode <- as.character(d$Tumor_Sample_Barcode)
  d$Hugo_Symbol <- as.character(d$Hugo_Symbol)
  d$Variant_Classification <- as.character(d$Variant_Classification)
  d
}

#' The six-class base change of each SNV (pyrimidine reference)
#' @param ref,alt Reference and alternate alleles.
#' @return Character vector ("C>T", ...), NA for non-SNVs.
#' @keywords internal
wes_sbs6 <- function(ref, alt) {
  ref <- toupper(as.character(ref))
  alt <- toupper(as.character(alt))
  ok <- nchar(ref) == 1 & nchar(alt) == 1 & ref %in% c("A", "C", "G", "T") & alt %in% c("A", "C", "G", "T") & ref != alt
  comp <- c(A = "T", C = "G", G = "C", T = "A")
  flip <- ok & ref %in% c("G", "A")
  ref[flip] <- comp[ref[flip]]
  alt[flip] <- comp[alt[flip]]
  out <- paste0(ref, ">", alt)
  out[!ok] <- NA
  out
}

# fine grey grid of tile backgrounds, used by the matrix figures
wes_tile_theme <- function() {
  bl <- ggplot2::element_blank()
  omicone_theme(grid = "none") +
    ggplot2::theme(axis.line = bl, axis.ticks = bl, panel.background = bl)
}

# ---- cohort summary ------------------------------------------------------------

#' Cohort summary dashboard (replaces maftools::plotmafSummary)
#'
#' Six panels: variant classifications, variant types, SNV classes, variants
#' per sample (stacked by classification, with the median or mean), the
#' per-sample distribution of each classification, and the most mutated genes.
#' @param maf A MAF object. @param top Genes in the last panel.
#' @param stat "median", "mean" or "none" for the per-sample line.
#' @param rm_outlier Hide boxplot outliers.
#' @keywords internal
wes_summary_gg <- function(maf, top = 10, stat = "median", rm_outlier = TRUE) {
  tk <- style_tokens()
  d <- wes_records(maf)
  vcs <- wes_vc_scale_values(d$Variant_Classification)
  d$vc <- factor(d$Variant_Classification, levels = rev(names(vcs)))
  n_samples <- length(wes_samples(maf))
  bars <- function(df, x, fill = NULL, title, values = NULL, flip = TRUE) {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x]], y = .data$n)) +
      ggplot2::labs(x = NULL, y = NULL, title = title) +
      omicone_theme(base_size = 11, grid = if (flip) "x" else "y")
    p <- if (is.null(fill)) p + ggplot2::geom_col(fill = tk$kept, width = 0.7)
         else p + ggplot2::geom_col(ggplot2::aes(fill = .data[[fill]]), width = 0.7) +
           ggplot2::scale_fill_manual(values = values, guide = "none")
    if (flip) p + ggplot2::coord_flip() else p
  }
  vc_tab <- as.data.frame(table(vc = d$vc), stringsAsFactors = FALSE)
  names(vc_tab)[2] <- "n"
  vc_tab$vc <- factor(vc_tab$vc, levels = rev(names(vcs)))
  p1 <- bars(vc_tab, "vc", "vc", "Variant classification", vcs) +
    ggplot2::scale_x_discrete(labels = wes_vc_label)
  vt <- as.data.frame(table(type = as.character(d$Variant_Type)), stringsAsFactors = FALSE)
  names(vt)[2] <- "n"
  vt$type <- factor(vt$type, levels = vt$type[order(vt$n)])
  p2 <- bars(vt, "type", NULL, "Variant type")
  s6 <- wes_sbs6(d$Reference_Allele, d$Tumor_Seq_Allele2)
  st <- as.data.frame(table(cls = factor(s6, levels = names(wes_sbs_colors()))), stringsAsFactors = FALSE)
  names(st)[2] <- "n"
  st$cls <- factor(st$cls, levels = rev(names(wes_sbs_colors())))
  p3 <- bars(st, "cls", "cls", "SNV class", wes_sbs_colors())
  per <- as.data.frame(table(sample = d$Tumor_Sample_Barcode, vc = d$vc), stringsAsFactors = FALSE)
  names(per)[3] <- "n"
  tot <- tapply(per$n, per$sample, sum)
  per$sample <- factor(per$sample, levels = names(sort(tot, decreasing = TRUE)))
  per$vc <- factor(per$vc, levels = rev(names(vcs)))
  line <- switch(stat, median = stats::median(tot), mean = mean(tot), NA)
  p4 <- ggplot2::ggplot(per, ggplot2::aes(.data$sample, .data$n, fill = .data$vc)) +
    ggplot2::geom_col(width = 1) +
    ggplot2::scale_fill_manual(values = vcs, guide = "none") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::labs(x = sprintf("%d samples, sorted", nlevels(per$sample)), y = NULL,
                  title = "Variants per sample") +
    omicone_theme(base_size = 11, grid = "y") +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
  if (is.finite(line)) {
    p4 <- p4 + ggplot2::geom_hline(yintercept = line, linetype = 2, colour = tk$muted, linewidth = 0.4) +
      ggplot2::annotate("text", x = Inf, y = line, hjust = 1.05, vjust = -0.5, size = 3, colour = tk$muted,
                        label = sprintf("%s %s", stat, format(round(line, 1))))
  }
  p5 <- ggplot2::ggplot(per, ggplot2::aes(.data$vc, .data$n, fill = .data$vc)) +
    ggplot2::geom_boxplot(outlier.shape = if (isTRUE(rm_outlier)) NA else 19, outlier.size = 0.6,
                          colour = tk$text, linewidth = 0.3, width = 0.65) +
    ggplot2::scale_fill_manual(values = vcs, guide = "none") +
    ggplot2::scale_x_discrete(labels = wes_vc_label) +
    ggplot2::coord_flip(ylim = if (isTRUE(rm_outlier)) c(0, max(1, stats::quantile(per$n, 0.98))) else NULL) +
    ggplot2::labs(x = NULL, y = "variants per sample", title = "Classification per sample") +
    omicone_theme(base_size = 11, grid = "x")
  gs <- as.data.frame(maftools::getGeneSummary(maf))
  gs <- utils::head(gs[order(-gs$AlteredSamples), , drop = FALSE], top)
  gtab <- d[d$Hugo_Symbol %in% gs$Hugo_Symbol, , drop = FALSE]
  gg <- as.data.frame(table(gene = gtab$Hugo_Symbol, vc = gtab$vc), stringsAsFactors = FALSE)
  names(gg)[3] <- "n"
  gg$gene <- factor(gg$gene, levels = rev(gs$Hugo_Symbol))
  gg$vc <- factor(gg$vc, levels = rev(names(vcs)))
  lab <- data.frame(gene = factor(gs$Hugo_Symbol, levels = rev(gs$Hugo_Symbol)),
                    total = tapply(gg$n, gg$gene, sum)[gs$Hugo_Symbol],
                    pct = sprintf("%.0f%%", 100 * gs$AlteredSamples / max(1, n_samples)))
  p6 <- ggplot2::ggplot(gg, ggplot2::aes(.data$gene, .data$n, fill = .data$vc)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_text(data = lab, ggplot2::aes(x = .data$gene, y = .data$total, label = .data$pct),
                       inherit.aes = FALSE, hjust = -0.15, size = 3, colour = tk$text) +
    ggplot2::scale_fill_manual(values = vcs, guide = "none") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.18))) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "variants", title = sprintf("Top %d mutated genes", nrow(gs))) +
    omicone_theme(base_size = 11, grid = "x") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(face = "italic"))
  leg <- ggplot2::ggplot(data.frame(vc = factor(names(vcs), levels = names(vcs)), n = 1),
                         ggplot2::aes(.data$vc, .data$n, fill = .data$vc)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = vcs, labels = wes_vc_label, name = NULL) +
    omicone_theme() + ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2))
  compose_grid(list(p1, p2, p3, p4, p5, p6), nrow = 2, ncol = 3, legends = list(leg),
               title = "Cohort summary",
               subtitle = sprintf("%s non-synonymous variants in %s samples", format(nrow(d), big.mark = ","),
                                  format(n_samples, big.mark = ",")))
}

# ---- oncoplot --------------------------------------------------------------------

#' Mutation matrix of some genes: one row per mutated sample x gene
#' @param d [wes_records()] output. @param genes Genes.
#' @return data.frame(sample, gene, vc): several classes in one gene and sample
#'   become "Multi_Hit".
#' @keywords internal
wes_mut_cells <- function(d, genes) {
  m <- d[d$Hugo_Symbol %in% genes, c("Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification"), drop = FALSE]
  if (!nrow(m)) return(data.frame(sample = character(0), gene = character(0), vc = character(0)))
  key <- paste(m$Tumor_Sample_Barcode, m$Hugo_Symbol, sep = "\r")
  vc <- tapply(m$Variant_Classification, key, function(x) if (length(unique(x)) > 1) "Multi_Hit" else x[1])
  parts <- strsplit(names(vc), "\r", fixed = TRUE)
  data.frame(sample = vapply(parts, `[`, "", 1), gene = vapply(parts, `[`, "", 2), vc = unname(vc),
             stringsAsFactors = FALSE)
}

#' Oncoplot (replaces maftools::oncoplot)
#'
#' Waterfall-sorted mutation matrix with the per-sample burden on top, the
#' share of mutated samples per gene on the right, clinical annotation rows
#' underneath and, optionally, each sample's SNV spectrum. The share is out of
#' every sample in the MAF; only samples with a shown gene mutated are drawn.
#' @param maf A MAF object. @param genes Genes, most frequent first.
#' @param clin Clinical columns to annotate. @param sort_anno Sort samples by the
#'   first annotation before the waterfall. @param show_pct Write the shares.
#' @param titv Add the SNV-spectrum row.
#' @keywords internal
wes_oncoplot_gg <- function(maf, genes, clin = NULL, sort_anno = FALSE, show_pct = TRUE, titv = FALSE) {
  tk <- style_tokens()
  d <- wes_records(maf)
  cells <- wes_mut_cells(d, genes)
  if (!nrow(cells)) stop("None of these genes is mutated in this cohort.")
  n_all <- length(wes_samples(maf))
  freq <- table(factor(cells$gene, levels = genes))
  genes <- names(sort(freq, decreasing = TRUE))
  genes <- genes[freq[genes] > 0]
  samples <- unique(cells$sample)
  bin <- matrix(0L, length(samples), length(genes), dimnames = list(samples, genes))
  bin[cbind(cells$sample, cells$gene)] <- 1L
  cd <- tryCatch(as.data.frame(maftools::getClinicalData(maf)), error = function(e) NULL)
  clin <- intersect(clin %||% character(0), names(cd %||% data.frame()))
  keys <- lapply(seq_along(genes), function(j) -bin[, j])
  if (isTRUE(sort_anno) && length(clin)) {
    a <- as.character(cd[[clin[1]]][match(samples, as.character(cd$Tumor_Sample_Barcode))])
    keys <- c(list(match(a, level_order(a))), keys)
  }
  samples <- samples[do.call(order, keys)]
  vcs <- wes_vc_scale_values(c(cells$vc, d$Variant_Classification))
  grid <- expand.grid(sample = samples, gene = genes, stringsAsFactors = FALSE)
  cells$sample <- factor(cells$sample, levels = samples)
  cells$gene <- factor(cells$gene, levels = rev(genes))
  grid$sample <- factor(grid$sample, levels = samples)
  grid$gene <- factor(grid$gene, levels = rev(genes))
  gene_size <- if (length(genes) > 40) 7 else if (length(genes) > 25) 8.5 else 10
  bl <- ggplot2::element_blank()
  p_main <- ggplot2::ggplot(cells, ggplot2::aes(.data$sample, .data$gene)) +
    ggplot2::geom_tile(data = grid, fill = "#eef1f5", colour = "white", linewidth = 0.25) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$vc), colour = "white", linewidth = 0.25) +
    ggplot2::scale_fill_manual(values = vcs, labels = wes_vc_label, name = NULL, drop = TRUE) +
    ggplot2::scale_x_discrete(expand = c(0, 0)) + ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::labs(x = sprintf("%d of %d samples (with a shown gene mutated)", length(samples), n_all), y = NULL) +
    wes_tile_theme() +
    ggplot2::theme(axis.text.x = bl, axis.text.y = ggplot2::element_text(face = "italic", size = gene_size),
                   legend.position = "none")
  burden <- d[d$Tumor_Sample_Barcode %in% samples, , drop = FALSE]
  bt <- as.data.frame(table(sample = factor(burden$Tumor_Sample_Barcode, levels = samples),
                            vc = burden$Variant_Classification), stringsAsFactors = FALSE)
  names(bt)[3] <- "n"
  bt$sample <- factor(bt$sample, levels = samples)
  bt$vc <- factor(bt$vc, levels = rev(names(vcs)))
  p_top <- ggplot2::ggplot(bt, ggplot2::aes(.data$sample, .data$n, fill = .data$vc)) +
    ggplot2::geom_col(width = 0.9) +
    ggplot2::scale_fill_manual(values = vcs, guide = "none") +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05)), breaks = scales::breaks_pretty(3)) +
    ggplot2::labs(x = NULL, y = "variants") +
    omicone_theme(base_size = 10, grid = "y") +
    ggplot2::theme(axis.text.x = bl, axis.ticks.x = bl, axis.line.x = bl, plot.margin = ggplot2::margin(6, 12, 0, 10))
  rt <- as.data.frame(table(gene = factor(cells$gene, levels = rev(genes)), vc = cells$vc), stringsAsFactors = FALSE)
  names(rt)[3] <- "n"
  rt$gene <- factor(rt$gene, levels = rev(genes))
  rt$vc <- factor(rt$vc, levels = rev(names(vcs)))
  rt$pct <- 100 * rt$n / max(1, n_all)
  lab <- data.frame(gene = factor(genes, levels = rev(genes)), pct = 100 * as.numeric(freq[genes]) / max(1, n_all))
  p_right <- ggplot2::ggplot(rt, ggplot2::aes(.data$gene, .data$pct, fill = .data$vc)) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::scale_fill_manual(values = vcs, guide = "none") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, if (isTRUE(show_pct)) 0.45 else 0.05)),
                                breaks = scales::breaks_pretty(2), labels = function(v) paste0(v, "%")) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "of all samples") +
    omicone_theme(base_size = 10, grid = "x") +
    ggplot2::theme(axis.text.y = bl, axis.ticks.y = bl, axis.line.y = bl, plot.margin = ggplot2::margin(0, 12, 8, 0))
  if (isTRUE(show_pct)) {
    p_right <- p_right + ggplot2::geom_text(data = lab, ggplot2::aes(x = .data$gene, y = .data$pct,
                                                                     label = sprintf("%.0f%%", .data$pct)),
                                            inherit.aes = FALSE, hjust = -0.15, size = 2.9, colour = tk$text)
  }
  rows <- list(p_top, NULL, p_main, p_right)
  heights <- c(1.1, max(2.5, 0.32 * length(genes)))
  legends <- list(p_main + ggplot2::theme(legend.position = "bottom") +
                    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2)))
  for (cl in clin) {
    v <- as.character(cd[[cl]][match(samples, as.character(cd$Tumor_Sample_Barcode))])
    cols <- value_colors(v)
    a <- data.frame(sample = factor(samples, levels = samples), feature = cl, value = v)
    pa <- ggplot2::ggplot(a, ggplot2::aes(.data$sample, .data$feature, fill = .data$value)) +
      ggplot2::geom_tile(colour = "white", linewidth = 0.25) +
      scale_group(cols, "fill", name = cl) +
      ggplot2::scale_x_discrete(expand = c(0, 0)) + ggplot2::scale_y_discrete(expand = c(0, 0)) +
      ggplot2::labs(x = NULL, y = NULL) +
      wes_tile_theme() +
      ggplot2::theme(axis.text.x = bl, legend.position = "none", plot.margin = ggplot2::margin(1, 12, 1, 10))
    rows <- c(rows, list(pa, NULL))
    heights <- c(heights, 0.32)
    legends <- c(legends, list(pa + ggplot2::theme(legend.position = "bottom") +
                                 ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2))))
  }
  if (isTRUE(titv)) {
    s6 <- wes_sbs6(d$Reference_Allele, d$Tumor_Seq_Allele2)
    tv <- d[!is.na(s6) & d$Tumor_Sample_Barcode %in% samples, , drop = FALSE]
    tv$cls <- s6[!is.na(s6) & d$Tumor_Sample_Barcode %in% samples]
    tt <- as.data.frame(table(sample = factor(tv$Tumor_Sample_Barcode, levels = samples),
                              cls = factor(tv$cls, levels = names(wes_sbs_colors()))), stringsAsFactors = FALSE)
    names(tt)[3] <- "n"
    tt$sample <- factor(tt$sample, levels = samples)
    pt <- ggplot2::ggplot(tt, ggplot2::aes(.data$sample, .data$n, fill = .data$cls)) +
      ggplot2::geom_col(position = "fill", width = 1) +
      ggplot2::scale_fill_manual(values = wes_sbs_colors(), name = "SNV class") +
      ggplot2::scale_x_discrete(expand = c(0, 0)) +
      ggplot2::scale_y_continuous(labels = function(v) paste0(100 * v, "%"), breaks = c(0, 0.5, 1)) +
      ggplot2::labs(x = NULL, y = "SNVs") +
      omicone_theme(base_size = 10, grid = "none") +
      ggplot2::theme(axis.text.x = bl, axis.ticks.x = bl, legend.position = "none")
    rows <- c(rows, list(pt, NULL))
    heights <- c(heights, 1)
    legends <- c(legends, list(pt + ggplot2::theme(legend.position = "bottom")))
  }
  compose_grid(rows, nrow = length(rows) / 2, ncol = 2, widths = c(6, 1.1), heights = heights,
               legends = legends, title = sprintf("Oncoplot: %d genes", length(genes)))
}

# ---- base-change spectra -------------------------------------------------------------

#' Transitions and transversions (replaces maftools::plotTiTv)
#' @param tv A [wes_titv()] result.
#' @keywords internal
wes_titv_gg <- function(tv) {
  tk <- style_tokens()
  fc <- as.data.frame(tv$fraction.contribution)
  cls <- intersect(names(wes_sbs_colors()), names(fc))
  long <- data.frame(sample = rep(as.character(fc$Tumor_Sample_Barcode), length(cls)),
                     cls = factor(rep(cls, each = nrow(fc)), levels = cls),
                     pct = unlist(fc[cls], use.names = FALSE))
  p1 <- ggplot2::ggplot(long, ggplot2::aes(.data$cls, .data$pct, fill = .data$cls)) +
    ggplot2::geom_boxplot(outlier.shape = NA, colour = tk$text, linewidth = 0.3, width = 0.6) +
    ggplot2::geom_jitter(width = 0.15, height = 0, size = 0.5, alpha = 0.35, colour = tk$ink) +
    ggplot2::scale_fill_manual(values = wes_sbs_colors(), guide = "none") +
    ggplot2::labs(x = NULL, y = "% of the sample's SNVs", title = "Six base-change classes") +
    omicone_theme(grid = "y")
  tt <- as.data.frame(tv$TiTv.fractions)
  tl <- data.frame(kind = factor(rep(c("Ti", "Tv"), each = nrow(tt)), levels = c("Ti", "Tv")),
                   pct = c(tt$Ti, tt$Tv))
  p2 <- ggplot2::ggplot(tl, ggplot2::aes(.data$kind, .data$pct, fill = .data$kind)) +
    ggplot2::geom_boxplot(outlier.shape = NA, colour = tk$text, linewidth = 0.3, width = 0.6) +
    ggplot2::geom_jitter(width = 0.15, height = 0, size = 0.5, alpha = 0.35, colour = tk$ink) +
    ggplot2::scale_fill_manual(values = c(Ti = "#E62725", Tv = "#1EBFF0"), guide = "none") +
    ggplot2::scale_x_discrete(labels = c(Ti = "transitions", Tv = "transversions")) +
    ggplot2::labs(x = NULL, y = "%", title = "Ti / Tv") +
    omicone_theme(grid = "y")
  ord <- fc$Tumor_Sample_Barcode[order(-fc[["C>T"]] %||% 0)]
  long$sample <- factor(long$sample, levels = as.character(ord))
  p3 <- ggplot2::ggplot(long, ggplot2::aes(.data$sample, .data$pct, fill = .data$cls)) +
    ggplot2::geom_col(width = 1) +
    ggplot2::scale_fill_manual(values = wes_sbs_colors(), name = NULL) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::labs(x = sprintf("%d samples, sorted by C>T", nrow(fc)), y = "%", title = "Per sample") +
    omicone_theme(grid = "none") +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
                   legend.position = "bottom") +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1))
  top <- compose_grid(list(p1, p2), nrow = 1, ncol = 2, widths = c(3, 1.3))
  compose_grid(list(top, p3), nrow = 2, ncol = 1, heights = c(1.1, 1), title = "Mutation spectrum")
}

#' VAF of the most mutated genes (replaces maftools::plotVaf)
#' @param maf A MAF object. @param vaf_col VAF column, or NULL to compute
#'   t_alt / (t_ref + t_alt). @param top Genes.
#' @keywords internal
wes_vaf_gg <- function(maf, vaf_col = NULL, top = 10) {
  tk <- style_tokens()
  d <- wes_records(maf)
  v <- if (!is.null(vaf_col) && vaf_col %in% names(d)) suppressWarnings(as.numeric(d[[vaf_col]]))
       else if ("t_vaf" %in% names(d)) suppressWarnings(as.numeric(d$t_vaf))
       else {
         a <- suppressWarnings(as.numeric(d$t_alt_count))
         r <- suppressWarnings(as.numeric(d$t_ref_count))
         a / (a + r)
       }
  if (all(is.na(v))) stop("No allele fractions: pick a VAF column, or provide t_ref_count / t_alt_count in the MAF.")
  if (any(v > 1, na.rm = TRUE)) v <- v / 100
  d$vaf <- v
  gs <- as.data.frame(maftools::getGeneSummary(maf))
  genes <- utils::head(gs$Hugo_Symbol[order(-gs$AlteredSamples)], top)
  x <- d[d$Hugo_Symbol %in% genes & !is.na(d$vaf), , drop = FALSE]
  med <- tapply(x$vaf, x$Hugo_Symbol, stats::median)
  x$gene <- factor(x$Hugo_Symbol, levels = names(sort(med, decreasing = TRUE)))
  vcs <- wes_vc_scale_values(x$Variant_Classification)
  ggplot2::ggplot(x, ggplot2::aes(.data$gene, .data$vaf)) +
    ggplot2::geom_boxplot(outlier.shape = NA, fill = tk$panel, colour = tk$faint, width = 0.6) +
    ggplot2::geom_jitter(ggplot2::aes(colour = .data$Variant_Classification), width = 0.15, height = 0,
                         size = 1.5, alpha = 0.8) +
    ggplot2::scale_colour_manual(values = vcs, labels = wes_vc_label, name = NULL) +
    ggplot2::scale_y_continuous(limits = c(0, 1), labels = function(v) paste0(100 * v, "%")) +
    ggplot2::labs(x = NULL, y = "variant allele fraction",
                  title = sprintf("VAF of the %d most mutated genes", length(genes)),
                  subtitle = "Clonal heterozygous mutations sit near purity / 2; lower values suggest subclones") +
    omicone_theme(grid = "y") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(face = "italic", angle = 45, hjust = 1),
                   legend.position = "bottom")
}

#' Chromosome lengths (bp) of hg19 / hg38, chr1-22, X, Y
#' @keywords internal
wes_chrom_lengths <- function(build = c("hg19", "hg38")) {
  build <- match.arg(build)
  len <- if (build == "hg19") {
    c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022,
      141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753,
      81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566)
  } else {
    c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636,
      138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345,
      83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895, 57227415)
  }
  stats::setNames(len, c(as.character(1:22), "X", "Y"))
}

#' Rainfall plot of one sample (replaces maftools::rainfallPlot)
#'
#' Each SNV at its genomic position against the distance to the previous SNV
#' (log scale); kataegis shows as a run of close, same-coloured points. With
#' `changepoints = TRUE`, maftools' change-point detection marks the regions.
#' @param maf A MAF object. @param tsb Sample. @param build "hg19" or "hg38".
#' @param changepoints Detect kataegis with maftools (needs 'changepoint').
#' @keywords internal
wes_rainfall_gg <- function(maf, tsb, build = "hg19", changepoints = FALSE) {
  tk <- style_tokens()
  d <- wes_records(maf, silent = TRUE)
  d <- d[d$Tumor_Sample_Barcode == tsb, , drop = FALSE]
  d$cls <- wes_sbs6(d$Reference_Allele, d$Tumor_Seq_Allele2)
  d <- d[!is.na(d$cls), , drop = FALSE]
  len <- wes_chrom_lengths(build)
  d$chr <- sub("^chr", "", as.character(d$Chromosome), ignore.case = TRUE)
  d$chr[d$chr == "23"] <- "X"
  d$chr[d$chr == "24"] <- "Y"
  d <- d[d$chr %in% names(len), , drop = FALSE]
  if (nrow(d) < 2) stop("Fewer than two SNVs in this sample: nothing to space out.")
  d <- d[!duplicated(paste(d$chr, d$Start_Position)), , drop = FALSE]     # one point per genomic position
  off <- c(0, cumsum(as.numeric(len)))[seq_along(len)]
  names(off) <- names(len)
  d$pos <- as.numeric(d$Start_Position) + off[d$chr]
  d <- d[order(d$pos), , drop = FALSE]
  d$dist <- c(NA, diff(d$pos))
  d$dist[c(TRUE, d$chr[-1] != d$chr[-nrow(d)])] <- NA
  d <- d[!is.na(d$dist) & d$dist > 0, , drop = FALSE]
  shade <- data.frame(xmin = off, xmax = off + len, chr = names(len))[seq(2, length(len), 2), ]
  p <- ggplot2::ggplot(d, ggplot2::aes(.data$pos, .data$dist)) +
    ggplot2::geom_rect(data = shade, ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = 0, ymax = Inf),
                       inherit.aes = FALSE, fill = tk$panel) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$cls), size = 1.3, alpha = 0.85) +
    ggplot2::scale_colour_manual(values = wes_sbs_colors(), name = NULL, drop = FALSE) +
    ggplot2::scale_y_log10(labels = function(v) format(v, big.mark = ",", scientific = FALSE)) +
    ggplot2::scale_x_continuous(breaks = off + len / 2, labels = names(len), limits = c(0, sum(as.numeric(len))),
                                expand = c(0.005, 0)) +
    ggplot2::labs(x = sprintf("chromosome (%s)", build), y = "distance to the previous SNV (bp, log)",
                  title = sprintf("Rainfall: %s", tsb),
                  subtitle = sprintf("%s SNVs; runs of close, same-coloured points suggest kataegis", format(nrow(d), big.mark = ","))) +
    omicone_theme(grid = "y") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 7.5), legend.position = "bottom") +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 1, override.aes = list(size = 3)))
  if (isTRUE(changepoints) && has_pkg("changepoint")) {
    kt <- tryCatch(with_null_device(suppressMessages(maftools::rainfallPlot(maf = maf, tsb = tsb,
                                                                            detectChangePoints = TRUE,
                                                                            ref.build = build))),
                   error = function(e) NULL)
    if (is.data.frame(kt) && nrow(kt)) {
      kt <- as.data.frame(kt)
      kc <- sub("^chr", "", as.character(kt$Chromosome), ignore.case = TRUE)
      kt$x0 <- as.numeric(kt$Start_Position) + off[kc]
      kt$x1 <- as.numeric(kt$End_Position) + off[kc]
      p <- p + ggplot2::geom_rect(data = kt, ggplot2::aes(xmin = .data$x0, xmax = .data$x1, ymin = 0, ymax = Inf),
                                  inherit.aes = FALSE, fill = tk$highlight, alpha = 0.25) +
        ggplot2::labs(caption = sprintf("%d kataegis region(s) detected (maftools change points)", nrow(kt)))
    }
  }
  p
}

# ---- TMB versus TCGA ------------------------------------------------------------------

#' Cohort TMB beside the 33 TCGA MC3 cohorts (replaces maftools::tcgaCompare)
#' @param maf A MAF object. @param capture Capture size (Mb).
#' @param cohort Label of this cohort. @param log_scale Log y-axis.
#' @keywords internal
wes_tcga_gg <- function(maf, capture, cohort = "This cohort", log_scale = TRUE) {
  tk <- style_tokens()
  args <- list(maf = maf, cohortName = cohort, capture_size = capture, logscale = log_scale)
  if (wes_has_arg("tcgaCompare", "rm_zero")) args$rm_zero <- TRUE
  tc <- with_null_device(suppressMessages(suppressWarnings(do.call(maftools::tcgaCompare, args))))
  per <- as.data.frame(tc$mutation_burden_perSample)
  med <- as.data.frame(tc$median_mutation_burden)
  med <- med[order(med$Median_Mutations), , drop = FALSE]
  per$cohort <- factor(as.character(per$cohort), levels = med$Cohort)
  per <- per[order(per$cohort, per$total_perMB), , drop = FALSE]
  per$rank <- stats::ave(seq_len(nrow(per)), per$cohort, FUN = function(i) (seq_along(i) - 0.5) / length(i))
  per$x <- as.numeric(per$cohort) - 0.5 + per$rank * 0.8 + 0.1
  per$mine <- as.character(per$cohort) == cohort
  med$x <- seq_len(nrow(med))
  med$mine <- med$Cohort == cohort
  p <- ggplot2::ggplot(per, ggplot2::aes(.data$x, pmax(.data$total_perMB, 1e-3))) +
    ggplot2::geom_rect(data = med[seq(2, nrow(med), 2), , drop = FALSE],
                       ggplot2::aes(xmin = .data$x - 0.5, xmax = .data$x + 0.5, ymin = 0, ymax = Inf),
                       inherit.aes = FALSE, fill = tk$panel) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$mine), size = 0.55, alpha = 0.8, show.legend = FALSE) +
    ggplot2::geom_segment(data = med, ggplot2::aes(x = .data$x - 0.4, xend = .data$x + 0.4,
                                                   y = pmax(.data$Median_Mutations, 1e-3),
                                                   yend = pmax(.data$Median_Mutations, 1e-3)),
                          inherit.aes = FALSE, colour = tk$removed, linewidth = 0.7) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = tk$accent, `FALSE` = "#7f8a96")) +
    ggplot2::scale_x_continuous(breaks = med$x, labels = med$Cohort, expand = c(0.01, 0),
                                sec.axis = ggplot2::dup_axis(labels = med$Cohort_Size, name = NULL)) +
    ggplot2::labs(x = NULL, y = sprintf("TMB (mutations / Mb%s)", if (isTRUE(log_scale)) ", log" else ""),
                  title = "TMB beside the TCGA MC3 cohorts",
                  subtitle = sprintf("%s highlighted; red bars = cohort medians; top numbers = samples. TCGA: MC3 calls over 35.8 Mb.",
                                     cohort)) +
    omicone_theme(grid = "y") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
                   axis.text.x.top = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 0, size = 7,
                                                           colour = tk$muted))
  if (isTRUE(log_scale)) p <- p + ggplot2::scale_y_log10(labels = function(v) format(v, scientific = FALSE, drop0trailing = TRUE))
  p
}

# ---- lollipop -----------------------------------------------------------------------

#' Lollipop plot of one gene (replaces maftools::lollipopPlot)
#'
#' Mutations along the protein, over its Pfam domains (maftools' domain table).
#' Lollipop height = number of mutations at a residue, stacked by class.
#' @param maf A MAF object. @param gene Gene. @param aa_col Protein-change column.
#' @param tx RefSeq transcript (default: the longest).
#' @param label_at Residues to label. @param show_rate Show the mutated share.
#' @keywords internal
wes_lollipop_gg <- function(maf, gene, aa_col = NULL, tx = NULL, label_at = numeric(0), show_rate = TRUE) {
  tk <- style_tokens()
  aa_col <- aa_col %||% wes_guess_aa_col(maf)
  if (is.null(aa_col)) stop("No protein-change column (e.g. HGVSp_Short, Protein_Change) in this MAF.")
  d <- wes_records(maf)
  d <- d[d$Hugo_Symbol == gene, , drop = FALSE]
  if (!nrow(d)) stop(sprintf("%s has no non-synonymous mutation in this cohort.", gene))
  d$aa <- as.character(d[[aa_col]])
  d$pos <- wes_aa_position(d$aa)
  dom <- as.data.frame(readRDS(system.file("extdata", "protein_domains.RDs", package = "maftools")))
  dom <- dom[dom$HGNC == gene, , drop = FALSE]
  if (nrow(dom)) {
    tx <- tx %||% dom$refseq.ID[which.max(dom$aa.length)]
    dom <- dom[dom$refseq.ID == tx, , drop = FALSE]
  }
  len <- if (nrow(dom)) dom$aa.length[1] else max(d$pos, na.rm = TRUE)
  beyond <- sum(d$pos > len, na.rm = TRUE)
  d <- d[!is.na(d$pos) & d$pos <= len, , drop = FALSE]
  if (!nrow(d)) stop("No mutation could be placed on the protein (protein-change format or transcript).")
  agg <- stats::aggregate(list(n = rep(1, nrow(d))), by = list(pos = d$pos, vc = d$Variant_Classification), FUN = sum)
  vcs <- wes_vc_scale_values(agg$vc)
  agg <- agg[order(agg$pos, match(agg$vc, names(vcs))), , drop = FALSE]
  agg$top <- stats::ave(agg$n, agg$pos, FUN = cumsum)
  stem <- stats::aggregate(top ~ pos, agg, max)
  hmax <- max(stem$top)
  dom_cols <- value_colors(unique(dom$Label))
  bar_h <- max(0.6, hmax * 0.09)
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = stem, ggplot2::aes(x = .data$pos, xend = .data$pos, y = 0, yend = .data$top),
                          colour = tk$faint, linewidth = 0.4) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = len, ymin = -bar_h, ymax = 0), fill = "#dfe4ea") +
    ggplot2::geom_point(data = agg, ggplot2::aes(x = .data$pos, y = .data$top, colour = .data$vc),
                        size = 3, alpha = 0.95) +
    ggplot2::scale_colour_manual(values = vcs, labels = wes_vc_label, name = NULL)
  if (nrow(dom)) {
    p <- p + ggplot2::geom_rect(data = dom, ggplot2::aes(xmin = .data$Start, xmax = .data$End,
                                                         ymin = -bar_h * 1.35, ymax = bar_h * 0.35, fill = .data$Label),
                                colour = NA, alpha = 0.95) +
      ggplot2::scale_fill_manual(values = dom_cols, name = "domain")
  }
  if (length(label_at)) {
    lab <- d[d$pos %in% label_at, c("pos", "aa"), drop = FALSE]
    lab <- lab[!duplicated(lab$pos), , drop = FALSE]
    lab$aa <- sub("^p\\.", "", lab$aa)
    lab$y <- stem$top[match(lab$pos, stem$pos)]
    dense <- nrow(lab) >= 12
    p <- p + ggplot2::geom_text(data = lab, ggplot2::aes(x = .data$pos, y = .data$y, label = .data$aa),
                                angle = if (dense) 90 else 0, hjust = if (dense) -0.25 else 0.5,
                                vjust = if (dense) 0.5 else -1.1, size = 3, colour = tk$ink, check_overlap = TRUE)
  }
  n_samp <- length(unique(d$Tumor_Sample_Barcode))
  n_all <- length(wes_samples(maf))
  p + ggplot2::scale_y_continuous(breaks = scales::breaks_pretty(4), expand = ggplot2::expansion(mult = c(0.02, 0.18)),
                                  limits = c(-bar_h * 1.4, NA)) +
    ggplot2::scale_x_continuous(expand = c(0.01, 0)) +
    ggplot2::labs(x = sprintf("amino acid (protein length %s aa%s)", format(len, big.mark = ","),
                              if (!is.null(tx)) paste0(", ", tx) else ""),
                  y = "mutations", title = gene,
                  subtitle = paste0(if (isTRUE(show_rate)) sprintf("Mutated in %d of %d samples (%.1f%%). ", n_samp, n_all, 100 * n_samp / max(1, n_all)) else "",
                                    if (beyond) sprintf("%d mutation(s) beyond this transcript's length not drawn.", beyond) else "")) +
    omicone_theme(grid = "y") +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold.italic"), legend.position = "bottom",
                   legend.box = "vertical")
}

# ---- drivers and interactions ------------------------------------------------------

#' OncodriveCLUST-style driver plot (replaces maftools::plotOncodrive)
#' @param drv oncodrive() result. @param fdr FDR cut-off for labels.
#' @keywords internal
wes_oncodrive_gg <- function(drv, fdr = 0.1) {
  tk <- style_tokens()
  d <- as.data.frame(drv)
  d$y <- -log10(pmax(d$fdr, 1e-300))
  d$sig <- d$fdr < fdr
  ggplot2::ggplot(d, ggplot2::aes(.data$fract_muts_in_clusters, .data$y)) +
    ggplot2::geom_hline(yintercept = -log10(fdr), linetype = 2, colour = tk$faint) +
    ggplot2::geom_point(ggplot2::aes(size = .data$clusters, colour = .data$sig), alpha = 0.8) +
    ggplot2::geom_text(data = d[d$sig, , drop = FALSE], ggplot2::aes(label = .data$Hugo_Symbol), size = 3.2,
                       fontface = "italic", vjust = -1.1, colour = tk$ink, check_overlap = TRUE) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = tk$removed, `FALSE` = tk$faint),
                                 labels = c(`TRUE` = sprintf("FDR < %g", fdr), `FALSE` = "n.s."), name = NULL) +
    ggplot2::scale_size_area(max_size = 7, name = "clusters") +
    ggplot2::scale_x_continuous(limits = c(0, 1.05), labels = function(v) paste0(100 * v, "%")) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.03, 0.12))) +
    ggplot2::labs(x = "share of the gene's mutations in positional clusters", y = "-log10 FDR",
                  title = "Genes with clustered (hotspot) mutations",
                  subtitle = "Positional clustering finds oncogene-like hotspots; tumour suppressors with scattered truncations are missed") +
    omicone_theme()
}

#' Co-occurrence / mutual exclusivity of gene pairs (replaces the
#' somaticInteractions() drawing)
#' @param int [wes_interactions()] result. @param genes Gene order.
#' @param q Threshold marked with an asterisk (on pAdj when present).
#' @keywords internal
wes_interactions_gg <- function(int, genes = NULL, q = 0.05) {
  tk <- style_tokens()
  d <- as.data.frame(int)
  pcol <- if ("pAdj" %in% names(d)) "pAdj" else "pValue"
  genes <- genes %||% unique(c(d$gene1, d$gene2))
  genes <- genes[genes %in% c(d$gene1, d$gene2)]
  ix <- function(g) match(g, genes)
  a <- ix(d$gene1)
  b <- ix(d$gene2)
  d$row <- genes[pmax(a, b)]
  d$col <- genes[pmin(a, b)]
  d$score <- ifelse(d$Event == "Co_Occurence", 1, -1) * pmin(-log10(pmax(d[[pcol]], 1e-10)), 3)
  d$mark <- ifelse(d[[pcol]] < q, "*", ifelse(d[[pcol]] < 2 * q, "\u00b7", ""))
  d$row <- factor(d$row, levels = rev(genes))
  d$col <- factor(d$col, levels = genes)
  ggplot2::ggplot(d, ggplot2::aes(.data$col, .data$row, fill = .data$score)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = .data$mark), size = 5, colour = tk$ink, vjust = 0.75) +
    ggplot2::scale_fill_gradientn(colours = style_div_colors(), limits = c(-3, 3), breaks = c(-3, -1.3, 0, 1.3, 3),
                                  labels = c(">3  exclusive", "1.3", "0", "1.3", ">3  co-occurring"),
                                  name = sprintf("-log10 %s", if (pcol == "pAdj") "FDR" else "p")) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::coord_equal() +
    ggplot2::labs(x = NULL, y = NULL, title = "Somatic interactions",
                  caption = sprintf("Fisher's exact test per pair; * %s < %g, \u00b7 < %g", if (pcol == "pAdj") "FDR" else "p", q, 2 * q)) +
    wes_tile_theme() +
    ggplot2::theme(axis.text.x.top = ggplot2::element_text(angle = 90, hjust = 0, vjust = 0.5, face = "italic"),
                   axis.text.y = ggplot2::element_text(face = "italic"))
}

# ---- clinical ------------------------------------------------------------------------

#' Genes enriched in clinical groups (replaces maftools::plotEnrichmentResults)
#' @param sig Significant groupwise rows (Hugo_Symbol, Group1, n_mutated_group1,
#'   n_mutated_group2, OR, fdr). @param feature Clinical column name.
#' @keywords internal
wes_enrichment_gg <- function(sig, feature = "") {
  tk <- style_tokens()
  d <- as.data.frame(sig)
  frac <- function(x) {
    p <- strsplit(as.character(x), " of ", fixed = TRUE)
    vapply(p, function(v) as.numeric(v[1]) / as.numeric(v[2]), 0)
  }
  d$in_group <- frac(d$n_mutated_group1)
  d$rest <- frac(d$n_mutated_group2)
  d$label <- sprintf("%s  (%s)", d$Hugo_Symbol, d$Group1)
  d <- d[order(d$Group1, -d$in_group), , drop = FALSE]
  d$label <- factor(d$label, levels = rev(unique(d$label)))
  long <- rbind(data.frame(label = d$label, group = d$Group1, who = "in the group", pct = 100 * d$in_group),
                data.frame(label = d$label, group = d$Group1, who = "all other samples", pct = 100 * d$rest))
  gcols <- value_colors(d$Group1)
  long$fill <- ifelse(long$who == "in the group", long$group, "rest")
  ggplot2::ggplot(long, ggplot2::aes(.data$pct, .data$label, fill = .data$fill)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.7) +
    ggplot2::geom_text(data = d, ggplot2::aes(x = 100 * pmax(.data$in_group, .data$rest), y = .data$label,
                                              label = sprintf("FDR %s", formatC(.data$fdr, format = "g", digits = 2))),
                       inherit.aes = FALSE, hjust = -0.15, size = 2.9, colour = tk$muted) +
    ggplot2::scale_fill_manual(values = c(gcols, rest = tk$ns), breaks = c(names(gcols), "rest"),
                               labels = c(names(gcols), "all other samples"), name = feature) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.25)), labels = function(v) paste0(v, "%")) +
    ggplot2::labs(x = "samples mutated", y = NULL, title = sprintf("Genes enriched in a %s group", feature),
                  subtitle = "Each group against all other samples (Fisher's exact test, BH across genes x groups)") +
    omicone_theme(grid = "x") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(face = "italic"))
}

#' Oncogenic pathways (replaces the OncogenicPathways / pathways drawing)
#' @param pw Pathway table (Pathway, N, n_affected_genes, fraction_affected,
#'   Mutated_samples, Fraction_mutated_samples).
#' @keywords internal
wes_pathways_gg <- function(pw) {
  tk <- style_tokens()
  d <- as.data.frame(pw)
  d <- d[order(d$Fraction_mutated_samples), , drop = FALSE]
  d$Pathway <- factor(d$Pathway, levels = d$Pathway)
  ggplot2::ggplot(d, ggplot2::aes(100 * .data$Fraction_mutated_samples, .data$Pathway)) +
    ggplot2::geom_col(fill = tk$kept, width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%d / %d genes", .data$n_affected_genes, .data$N)),
                       hjust = -0.1, size = 3, colour = tk$muted) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.3)), labels = function(v) paste0(v, "%")) +
    ggplot2::labs(x = "samples with the pathway mutated", y = NULL, title = "Oncogenic signalling pathways",
                  subtitle = "Label: genes of the pathway mutated in at least one sample") +
    omicone_theme(grid = "x")
}

#' Druggable gene categories (replaces the drugInteractions() drawing)
#' @param dg drugInteractions() table (Gene, category).
#' @keywords internal
wes_drug_gg <- function(dg) {
  tk <- style_tokens()
  d <- as.data.frame(dg)
  n <- tapply(d$Gene, d$category, function(g) length(unique(g)))
  top <- lapply(split(d$Gene, d$category), function(g) paste(utils::head(unique(g), 5), collapse = ", "))
  t <- data.frame(category = names(n), n = as.numeric(n), genes = unlist(top[names(n)]))
  t <- t[order(t$n), , drop = FALSE]
  t$category <- factor(tools::toTitleCase(tolower(t$category)), levels = tools::toTitleCase(tolower(t$category)))
  ggplot2::ggplot(t, ggplot2::aes(.data$n, .data$category)) +
    ggplot2::geom_col(fill = tk$kept, width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = .data$genes), hjust = -0.05, size = 2.8, colour = tk$muted, fontface = "italic") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.9))) +
    ggplot2::labs(x = "mutated genes", y = NULL, title = "Druggable gene categories (DGIdb)",
                  subtitle = "A category is a lead to look up, not a treatment recommendation") +
    omicone_theme(grid = "x")
}

#' Two cohorts side by side (replaces maftools::coBarplot)
#' @param m1,m2 MAF objects. @param l1,l2 Names. @param genes Genes.
#' @keywords internal
wes_cobar_gg <- function(m1, m2, l1, l2, genes) {
  side <- function(m, lab, sign) {
    d <- wes_records(m)
    cells <- wes_mut_cells(d, genes)
    n <- length(wes_samples(m))
    t <- as.data.frame(table(gene = factor(cells$gene, levels = genes), vc = cells$vc), stringsAsFactors = FALSE)
    names(t)[3] <- "k"
    t$pct <- sign * 100 * t$k / max(1, n)
    t$cohort <- sprintf("%s (n = %d)", lab, n)
    t
  }
  a <- side(m1, l1, -1)
  b <- side(m2, l2, 1)
  d <- rbind(a, b)
  vcs <- wes_vc_scale_values(d$vc)
  d$gene <- factor(d$gene, levels = rev(genes))
  d$vc <- factor(d$vc, levels = rev(names(vcs)))
  lim <- max(abs(tapply(d$pct, paste(d$gene, sign(d$pct)), sum)), 1)
  ggplot2::ggplot(d, ggplot2::aes(.data$pct, .data$gene, fill = .data$vc)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_vline(xintercept = 0, colour = style_tokens()$muted, linewidth = 0.4) +
    ggplot2::annotate("text", x = c(-lim, lim), y = length(genes) + 0.8, label = c(unique(a$cohort), unique(b$cohort)),
                      hjust = c(0, 1), size = 3.4, fontface = "bold", colour = style_tokens()$ink) +
    ggplot2::scale_fill_manual(values = vcs, labels = wes_vc_label, name = NULL) +
    ggplot2::scale_x_continuous(limits = c(-lim, lim) * 1.05, labels = function(v) paste0(abs(v), "%")) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = c(0.6, 1.3))) +
    ggplot2::labs(x = "samples mutated", y = NULL, title = "Mutation frequencies, cohort by cohort") +
    omicone_theme(grid = "x") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(face = "italic"), legend.position = "bottom")
}

# ---- heterogeneity -------------------------------------------------------------------

#' VAF clusters of one sample (replaces maftools::plotClusters)
#' @param het inferHeterogeneity() result. @param tsb Sample.
#' @keywords internal
wes_hetero_gg <- function(het, tsb) {
  tk <- style_tokens()
  d <- as.data.frame(het$clusterData)
  d <- d[as.character(d$Tumor_Sample_Barcode) == tsb, , drop = FALSE]
  if (!nrow(d)) stop("No clustered variants for this sample.")
  d$cluster <- as.character(d$cluster)
  lv <- level_order(d$cluster[d$cluster != "outlier"])
  cols <- c(value_colors(lv, levels = lv), outlier = tk$faint, CN_altered = tk$warn)
  d$cluster <- factor(d$cluster, levels = intersect(names(cols), unique(d$cluster)))
  mn <- as.data.frame(het$clusterMeans)
  mn <- mn[as.character(mn$Tumor_Sample_Barcode) == tsb & mn$cluster %in% lv, , drop = FALSE]
  math <- if ("MATH" %in% names(d)) unique(d$MATH)[1] else NA
  ggplot2::ggplot(d, ggplot2::aes(.data$t_vaf)) +
    ggplot2::geom_density(data = d[d$cluster %in% lv, , drop = FALSE], ggplot2::aes(fill = .data$cluster),
                          colour = NA, alpha = 0.45, adjust = 0.9) +
    ggplot2::geom_rug(ggplot2::aes(colour = .data$cluster), linewidth = 0.8, length = grid::unit(0.06, "npc")) +
    ggplot2::geom_vline(data = mn, ggplot2::aes(xintercept = .data$meanVaf, colour = .data$cluster), linetype = 2) +
    ggplot2::scale_fill_manual(values = cols, name = "cluster") +
    ggplot2::scale_colour_manual(values = cols, guide = "none") +
    ggplot2::scale_x_continuous(limits = c(0, 1), labels = function(v) paste0(100 * v, "%")) +
    ggplot2::labs(x = "variant allele fraction", y = "density", title = sprintf("VAF clusters: %s", tsb),
                  subtitle = sprintf("%d variants%s; dashed = cluster mean VAF; ticks = variants",
                                     nrow(d), if (is.finite(math)) sprintf(", MATH = %.1f", math) else "")) +
    omicone_theme(grid = "y")
}

# ---- signatures ----------------------------------------------------------------------

# the 96 trinucleotide contexts split into class and context ("A[C>A]A")
wes_sbs96_parts <- function(ctx) {
  cls <- sub("^.\\[(.>.)\\].$", "\\1", ctx)
  tri <- paste0(substr(ctx, 1, 1), substr(ctx, 3, 3), substr(ctx, 7, 7))
  data.frame(context = ctx, cls = cls, tri = tri, stringsAsFactors = FALSE)
}

#' Signature profiles on the 96 contexts (replaces maftools::plotSignatures)
#' @param sig extractSignatures() result. @param matches [wes_sig_matches()]
#'   output (best COSMIC match per signature) or NULL.
#' @keywords internal
wes_signatures_gg <- function(sig, matches = NULL) {
  w <- as.matrix(sig$signatures)
  w <- sweep(w, 2, pmax(colSums(w), 1e-12), "/")
  parts <- wes_sbs96_parts(rownames(w))
  ord <- order(match(parts$cls, names(wes_sbs_colors())), parts$tri)
  long <- do.call(rbind, lapply(colnames(w), function(s) {
    data.frame(sig = s, parts[ord, ], share = w[ord, s], stringsAsFactors = FALSE)
  }))
  long$context <- factor(long$context, levels = rownames(w)[ord])
  lab <- colnames(w)
  if (!is.null(matches) && nrow(matches)) {
    m <- matches[match(lab, matches$signature), , drop = FALSE]
    lab <- ifelse(is.na(m$best_match), lab,
                  sprintf("%s  -  best match %s (cosine %.2f)%s", lab, m$best_match, m$cosine,
                          ifelse(is.na(m$aetiology) | !nzchar(m$aetiology), "", paste0(": ", m$aetiology))))
  }
  long$sig <- factor(long$sig, levels = colnames(w), labels = lab)
  ggplot2::ggplot(long, ggplot2::aes(.data$context, 100 * .data$share, fill = .data$cls)) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::facet_wrap(~ sig, ncol = 1, scales = "free_y") +
    ggplot2::scale_fill_manual(values = wes_sbs_colors(), name = NULL) +
    ggplot2::scale_x_discrete(labels = function(x) wes_sbs96_parts(x)$tri) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.08))) +
    ggplot2::labs(x = NULL, y = "% of the signature", title = "Mutational signatures (SBS96)") +
    omicone_theme(base_size = 11, grid = "y") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, size = 5.5, family = "mono"),
                   legend.position = "top") +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1))
}

#' Signature exposures per sample (replaces plotSignatures(contributions = TRUE))
#' @param sig extractSignatures() result.
#' @keywords internal
wes_exposures_gg <- function(sig) {
  h <- as.matrix(sig$contributions)
  h <- sweep(h, 2, pmax(colSums(h), 1e-12), "/")
  dom <- apply(h, 2, which.max)
  ord <- order(dom, -apply(h, 2, max))
  long <- data.frame(sample = factor(rep(colnames(h), each = nrow(h)), levels = colnames(h)[ord]),
                     sig = factor(rep(rownames(h), ncol(h)), levels = rownames(h)),
                     share = as.vector(h))
  ggplot2::ggplot(long, ggplot2::aes(.data$sample, .data$share, fill = .data$sig)) +
    ggplot2::geom_col(width = 1) +
    scale_group(value_colors(rownames(h), levels = rownames(h)), "fill", name = NULL) +
    ggplot2::scale_y_continuous(expand = c(0, 0), labels = function(v) paste0(100 * v, "%")) +
    ggplot2::labs(x = sprintf("%d samples, grouped by their main signature", ncol(h)), y = "share of the sample's SNVs",
                  title = "Signature exposures") +
    omicone_theme(grid = "none") +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
                   legend.position = "bottom")
}

#' Cophenetic correlation over ranks (replaces maftools::plotCophenetic)
#' @param est estimateSignatures() result. @param best Chosen rank.
#' @keywords internal
wes_cophenetic_gg <- function(est, best = NULL) {
  tk <- style_tokens()
  d <- as.data.frame(est$nmfSummary)
  p <- ggplot2::ggplot(d, ggplot2::aes(.data$rank, .data$cophenetic)) +
    ggplot2::geom_line(colour = tk$faint) +
    ggplot2::geom_point(colour = tk$kept, size = 2.6) +
    ggplot2::scale_x_continuous(breaks = d$rank) +
    ggplot2::labs(x = "number of signatures (rank)", y = "cophenetic correlation",
                  title = "How many signatures?",
                  subtitle = "Pick the rank just before the correlation drops sharply") +
    omicone_theme(grid = "y")
  if (!is.null(best) && best %in% d$rank) {
    p <- p + ggplot2::geom_point(data = d[d$rank == best, , drop = FALSE], colour = tk$removed, size = 4.2, shape = 21,
                                 stroke = 1.2, fill = NA)
  }
  p
}

#' APOBEC-enriched versus other samples (replaces maftools::plotApobecDiff)
#'
#' Same analysis as maftools: samples are APOBEC-enriched when the one-sided
#' Fisher test calls them so (trinucleotideMatrix()), and gene frequencies are
#' compared between the two groups with mafCompare().
#' @param tnm trinucleotideMatrix() result. @param maf A MAF object.
#' @param p p-value threshold for the genes shown (as maftools: 0.05).
#' @keywords internal
wes_apobec_gg <- function(tnm, maf, p = 0.05) {
  tk <- style_tokens()
  s <- as.data.frame(tnm$APOBEC_scores)
  s <- s[!is.na(s$APOBEC_Enriched), , drop = FALSE]
  yes <- as.character(s$Tumor_Sample_Barcode[s$APOBEC_Enriched == "yes"])
  no <- as.character(s$Tumor_Sample_Barcode[s$APOBEC_Enriched == "no"])
  if (!length(yes)) stop("No sample is APOBEC-enriched (one-sided Fisher test): that mutational process is essentially absent here.")
  s$group <- factor(ifelse(s$APOBEC_Enriched == "yes", "enriched", "not enriched"), levels = c("enriched", "not enriched"))
  s <- s[order(-s$APOBEC_Enrichment), , drop = FALSE]
  s$rank <- seq_len(nrow(s))
  gcol <- c(enriched = tk$removed, `not enriched` = tk$ns)
  p1 <- ggplot2::ggplot(s, ggplot2::aes(.data$rank, .data$APOBEC_Enrichment, colour = .data$group)) +
    ggplot2::geom_hline(yintercept = 2, linetype = 2, colour = tk$faint) +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::scale_colour_manual(values = gcol, name = NULL) +
    ggplot2::labs(x = "samples, sorted", y = "APOBEC enrichment score", title = "APOBEC enrichment per sample",
                  subtitle = sprintf("%d of %d samples enriched", length(yes), nrow(s))) +
    omicone_theme(grid = "y") +
    ggplot2::theme(legend.position = "bottom")
  frac <- stats::aggregate(fraction_APOBEC_mutations ~ group, s, mean)
  p2 <- ggplot2::ggplot(frac, ggplot2::aes(.data$group, 100 * .data$fraction_APOBEC_mutations, fill = .data$group)) +
    ggplot2::geom_col(width = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.0f%%", 100 * .data$fraction_APOBEC_mutations)), vjust = -0.4,
                       size = 3.2, colour = tk$text) +
    ggplot2::scale_fill_manual(values = gcol, guide = "none") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
    ggplot2::labs(x = NULL, y = "mean % tCw mutations", title = "tCw load") +
    omicone_theme(grid = "y")
  m1 <- maftools::subsetMaf(maf = maf, tsb = yes, mafObj = TRUE)
  out <- list(p1, p2)
  if (length(no)) {
    m2 <- maftools::subsetMaf(maf = maf, tsb = no, mafObj = TRUE)
    mc <- tryCatch(maftools::mafCompare(m1 = m1, m2 = m2, m1Name = "enriched", m2Name = "not enriched", minMut = 2),
                   error = function(e) NULL)
    r <- if (!is.null(mc)) as.data.frame(mc$results) else NULL
    if (!is.null(r) && nrow(r) && any(r$pval < p)) {
      genes <- utils::head(r$Hugo_Symbol[order(r$pval)][r$pval[order(r$pval)] < p], 10)
      share <- function(m, n) {
        gs <- as.data.frame(maftools::getGeneSummary(m))
        v <- gs$MutatedSamples[match(genes, gs$Hugo_Symbol)]
        v[is.na(v)] <- 0
        100 * v / n
      }
      g <- rbind(data.frame(gene = genes, group = "enriched", pct = share(m1, length(yes))),
                 data.frame(gene = genes, group = "not enriched", pct = share(m2, length(no))))
      g$gene <- factor(g$gene, levels = genes)
      g$group <- factor(g$group, levels = names(gcol))
      out[[3]] <- ggplot2::ggplot(g, ggplot2::aes(.data$gene, .data$pct, fill = .data$group)) +
        ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.7) +
        ggplot2::scale_fill_manual(values = gcol, name = NULL) +
        ggplot2::labs(x = NULL, y = "samples mutated (%)", title = sprintf("Genes differing (Fisher p < %g)", p)) +
        omicone_theme(grid = "y") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(face = "italic", angle = 45, hjust = 1), legend.position = "none")
    }
  }
  if (length(out) == 3) {
    compose_grid(list(compose_grid(list(p1, p2), 1, 2, widths = c(2, 1)), out[[3]]), nrow = 2, ncol = 1,
                 title = "APOBEC-enriched samples")
  } else {
    compose_grid(list(p1, p2), 1, 2, widths = c(2, 1), title = "APOBEC-enriched samples",
                 subtitle = "No gene's mutation frequency differs between the two groups")
  }
}
