#' Module: Differential abundance of cell types
#'
#' Test whether the share of each cell type differs between conditions, with
#' the sample as the unit: proportions per sample, logit-transformed, compared
#' with a limma linear model (the propeller method, Phipson et al. 2022).
#' Cell counts pooled across samples (a chi-square on the total table) are not
#' used: they treat cells as replicates and ignore sample-to-sample variation.
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_abundance
NULL

#' @rdname mod_abundance
#' @keywords internal
mod_abundance_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = list(en = "Differential abundance", zh = "细胞组成差异"),
    what = list(
      en = "Ask whether a cell type makes up a larger or smaller share of the cells in one condition than in another.",
      zh = "检验某种细胞类型在一种条件下占细胞总数的比例，是否比另一种条件下更高或更低。"),
    why = list(
      en = "Composition is often the first biological signal: more CD8 T cells in responders, fewer B cells
            after treatment. But samples vary a lot on their own, and pooling all cells into one table
            hides that. Testing per-sample proportions keeps the patient as the replicate.",
      zh = "细胞组成往往是最先出现的生物学信号：应答者 CD8 T 细胞更多，治疗后 B 细胞更少。但样本之间本身差异就很大，把所有细胞合并成一张表会掩盖这一点。按样本计算比例再检验，才能以患者为重复单位。"),
    how = list(
      en = "Pick the <b>sample</b>, <b>condition</b> and <b>cell-type</b> columns. Each sample's proportions get a
            0.5 pseudo-count and a logit transform, then limma compares conditions with robust empirical
            Bayes (a moderated t test for two conditions, an F test for more); BH across cell types. This is
            the method of <code>speckle::propeller</code>, so speckle is not needed.",
      zh = "选择<b>样本</b>、<b>条件</b>和<b>细胞类型</b>列。每个样本的比例加 0.5 伪计数后做 logit 变换，再用 limma 的稳健经验贝叶斯比较条件（两种条件用 moderated t 检验，多于两种用 F 检验），并在细胞类型之间做 BH 校正。这正是 <code>speckle::propeller</code> 的方法，因此无需安装 speckle。"),
    read = list(
      en = "One panel per cell type: each dot is a sample's share of that cell type, by condition, with the
            FDR above. Look for a shift that most samples follow, not one extreme sample. Proportions are
            linked: if one cell type expands, the others shrink, so read several panels together.",
      zh = "每种细胞类型一个面板：每个点是一个样本中该类型细胞的占比，按条件分组，上方标注 FDR。要看大多数样本是否一致地偏移，而不是某一个极端样本。比例之间是相互关联的：一种细胞扩增，其他细胞的比例就会下降，所以要把几个面板放在一起看。"),
    example = list(
      en = "With 5 responders and 5 non-responders, CD8 T cells at 18% vs 9% of cells with FDR 0.03 is a
            composition difference between patients, not just between cells.",
      zh = "5 名应答者和 5 名非应答者中，CD8 T 细胞占比 18% 对 9%，FDR 0.03：这是患者之间的组成差异，而不只是细胞之间的差异。")
  )
  tbl_out <- if (has_pkg("DT")) DT::dataTableOutput(ns("tbl")) else shiny::verbatimTextOutput(ns("tbl"))
  controls <- shiny::tagList(
    sample_design_controls(ns, levels = FALSE),
    run_button(ns("run"), "Test abundance", "检验细胞组成")
  )
  step_container(id = id, 
    title     = list(en = "Differential abundance", zh = "细胞组成差异"),
    subtitle  = list(en = "Do cell-type proportions differ between conditions?",
                     zh = "细胞类型的比例在条件之间是否不同？"),
    explainer = explainer,
    controls  = controls,
    summary   = shiny::uiOutput(ns("summary")),
    preview   = shiny::tagList(
      shiny::uiOutput(ns("insight")),
      bslib::navset_card_tab(
        bslib::nav_panel(i18n("Proportions", "比例"),
          preview_plot_ui(ns("preview"), download = TRUE,
            guide = list(en = "Each sample's share of every cell type, by condition, will be drawn here.",
                         zh = "运行后，这里将按条件绘制每个样本中各细胞类型的占比。"),
            caption = list(en = "One dot = one sample; boxes summarise the samples of a condition; FDR from the propeller test (BH across cell types).",
                           zh = "每个点 = 一个样本；箱线总结同一条件的样本；FDR 来自 propeller 检验（细胞类型间 BH 校正）。"))),
        bslib::nav_panel(i18n("Composition", "组成"),
          preview_plot_ui(ns("stack"), download = TRUE,
            caption = list(en = "Cell-type composition of every sample, grouped by condition.",
                           zh = "每个样本的细胞类型组成，按条件分组。"))),
        bslib::nav_panel(i18n("Table", "表格"),
          shiny::div(class = "omicone-table",
                     shiny::downloadButton(ns("csv"), i18n("Download table (CSV)", "下载表格（CSV）"),
                                           class = "btn-sm"),
                     tbl_out))
      )
    )
  )
}

#' @rdname mod_abundance
#' @keywords internal
mod_abundance_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    res <- step_results(rv, "sc", df = NULL, design = NULL)
    design <- sample_design_server(input, output, session, rv, levels = FALSE)

    shiny::observeEvent(input$run, {
      shiny::req(rv$obj)
      d <- design()
      md <- obj_meta(rv$obj)
      if (!all(c(d$sample, d$condition, d$group) %in% names(md))) {
        shiny::showNotification(i18n("Choose the sample, condition and cell-type columns first.",
                                     "请先选择样本、条件和细胞类型列。"), type = "error")
        return(NULL)
      }
      if (!require_pkgs("limma", "Differential abundance")) return(NULL)
      df <- with_progress_notify({
        samples <- sample_table(md, d$sample, d$condition)
        abundance_test(md, d$sample, d$group, d$condition, samples = samples)
      }, message = "Testing cell-type proportions...")
      if (is.null(df)) return(NULL)
      res$df <- df
      res$design <- d
      mark_done(rv, "abundance")
      n <- attr(df, "n")
      log_step(log_rv, "Differential abundance",
               params = list(sample = d$sample, condition = d$condition, cell_type = d$group,
                             samples = paste(sprintf("%s = %d", names(n), as.integer(n)), collapse = ", "),
                             method = "propeller (logit, limma robust eBayes)"),
               code = abundance_log_code(d$sample, d$group, d$condition, n_levels = length(n)))
      shiny::showNotification(
        i18n(sprintf("%d of %d cell types differ at FDR < 0.05.", sum(df$fdr < 0.05), nrow(df)),
             sprintf("%d / %d 种细胞类型在 FDR < 0.05 下有差异。", sum(df$fdr < 0.05), nrow(df))),
        type = "message")
    })

    output$summary <- shiny::renderUI({
      df <- res$df
      if (is.null(df)) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("Choose the sample and condition columns, then click <b>Test abundance</b>.",
                               "选择样本列和条件列后，点击<b>检验细胞组成</b>。")))
      }
      n <- attr(df, "n")
      shiny::tagList(
        stat_tile(i18n("Samples", "样本"), paste(sprintf("%s %d", names(n), as.integer(n)), collapse = " · ")),
        stat_tile(i18n("Cell types", "细胞类型"), nrow(df)),
        stat_tile(i18n("FDR < 0.05", "FDR < 0.05"), sum(df$fdr < 0.05)),
        stat_tile(i18n("Test", "检验"), if (length(n) == 2) "moderated t" else "moderated F")
      )
    })

    output$insight <- shiny::renderUI({
      df <- res$df
      if (is.null(df)) return(NULL)
      n <- attr(df, "n")
      top <- df[1, ]
      small <- min(n) < 3
      dir <- if (length(n) == 2) {
        sprintf(" (%.1f%% vs %.1f%%, %s vs %s)", 100 * top[[3]], 100 * top[[2]], names(n)[2], names(n)[1])
      } else ""
      insight_bar(
        sprintf("%d of %d cell types differ at FDR < 0.05; the strongest is %s%s, FDR = %s.%s Proportions are linked: one population expanding makes the others shrink.",
                sum(df$fdr < 0.05), nrow(df), top$group, dir, formatC(top$fdr, format = "g", digits = 2),
                if (small) " With fewer than 3 samples in a condition the test has little power." else ""),
        sprintf("%d / %d 种细胞类型在 FDR < 0.05 下有差异；最明显的是 %s%s，FDR = %s。%s比例相互关联：一个细胞群扩增，其他细胞群的比例就会下降。",
                sum(df$fdr < 0.05), nrow(df), top$group, dir, formatC(top$fdr, format = "g", digits = 2),
                if (small) "某条件少于 3 个样本时，检验力很低。" else ""))
    })

    render_step_plot(output, input, "preview", function() {
      df <- res$df
      shiny::req(df)
      d <- res$design
      abundance_plot(df, cond_colors = if (d$condition %in% obj_meta_cols(rv$obj)) group_colors(rv$obj, d$condition))
    }, name = "abundance_proportions", width = 11, height = 7)

    render_step_plot(output, input, "stack", function() {
      df <- res$df
      shiny::req(df)
      d <- res$design
      abundance_stack_plot(df, group_colors = if (d$group %in% obj_meta_cols(rv$obj)) group_colors(rv$obj, d$group))
    }, name = "abundance_composition", width = 11, height = 6)

    output$tbl <- render_tbl_wrap(function() {
      df <- res$df
      shiny::req(df)
      df
    })

    output$csv <- shiny::downloadHandler(
      filename = function() sprintf("abundance_%s.csv", format(Sys.time(), "%Y%m%d_%H%M%S")),
      content = function(file) {
        df <- res$df
        utils::write.csv(if (is.null(df)) data.frame() else df, file, row.names = FALSE)
      }
    )
  })
}
