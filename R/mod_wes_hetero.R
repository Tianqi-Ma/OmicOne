#' WES module 12: Tumour heterogeneity
#'
#' Cluster one sample's variants by their allele frequency. Distinct VAF modes
#' correspond to distinct clones, and the spread of the highest mode gives the
#' MATH score — a single-number summary of how heterogeneous the tumour is.
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_wes_hetero
NULL

#' @rdname mod_wes_hetero
#' @keywords internal
mod_wes_hetero_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = list(en = "Heterogeneity", zh = "肿瘤异质性"),
    what = list(
      en = "Group one sample's mutations by the fraction of reads that carried
            them. Mutations present in every tumour cell sit at a high VAF;
            mutations acquired later by a subset sit lower.",
      zh = "把某个样本的突变按其支持读段占比分组。存在于所有肿瘤细胞中的突变 VAF 较高；后来才被一部分细胞获得的突变则较低。"),
    why  = list(
      en = "A tumour is rarely one clone. Knowing which mutations are clonal
            (present everywhere, and therefore a real therapeutic target) versus
            subclonal changes what you would treat. High heterogeneity is
            associated with worse outcome and resistance.",
      zh = "肿瘤很少只有一个克隆。区分克隆性突变（无处不在，因而是真正的治疗靶点）与亚克隆突变，会改变治疗决策。高异质性与更差的预后和耐药相关。"),
    how  = list(
      en = "This step needs a <b>VAF column</b> in the MAF; if the dropdown shows
            none, your caller did not report allele frequencies and this analysis
            cannot run. Copy-number-altered regions distort VAF, so treat the
            clusters as descriptive.",
      zh = "该步骤需要 MAF 中有 <b>VAF 列</b>；如果下拉框里没有可选项，说明你的变异检出工具未报告等位基因频率，本分析无法进行。拷贝数改变区域会扭曲 VAF，因此这些聚类结果只作描述性参考。"),
    example = list(
      en = "A MATH score above roughly 30 is commonly treated as a heterogeneous
               tumour; below ~20 the tumour is close to clonal.",
      zh = "MATH 分数大约高于 30 通常视为高异质性肿瘤；低于 20 左右则接近单克隆。")
  )
  controls <- shiny::tagList(
    shiny::uiOutput(ns("sample_ui")),
    shiny::uiOutput(ns("vaf_ui")),
    run_button(ns("run"), "Infer clones", "推断克隆结构")
  )
  step_container(
    title     = list(en = "Heterogeneity", zh = "肿瘤异质性"),
    explainer = explainer,
    controls  = controls,
    summary   = shiny::uiOutput(ns("summary")),
    preview   = bslib::navset_card_tab(
      bslib::nav_panel(i18n("Clusters", "克隆聚类"), preview_plot_ui(ns("plot"), download = TRUE)),
      bslib::nav_panel(i18n("Variants", "变异明细"), shiny::uiOutput(ns("tbl_slot")))
    )
  )
}

#' @rdname mod_wes_hetero
#' @keywords internal
mod_wes_hetero_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    ns  <- session$ns
    res <- shiny::reactiveValues(het = NULL, sample = NULL)

    output$sample_ui <- shiny::renderUI({
      shiny::req(rv$maf)
      s <- wes_samples(rv$maf)
      shiny::tagList(
        label_with_help("Sample", "Heterogeneity is inferred one sample at a time.",
                        label_zh = "样本", tip_zh = "异质性是逐个样本推断的。"),
        shiny::selectizeInput(ns("tsb"), NULL, choices = s,
                              selected = if (length(s)) s[1] else NULL)
      )
    })

    output$vaf_ui <- shiny::renderUI({
      shiny::req(rv$maf)
      f <- wes_fields(rv$maf)
      guess <- wes_guess_vaf_col(rv$maf)
      shiny::tagList(
        label_with_help("VAF column", "The MAF column with the variant allele frequency.",
                        label_zh = "VAF 列", tip_zh = "MAF 中存放变异等位基因频率的列。"),
        shiny::selectInput(ns("vaf_col"), NULL,
                           choices = c(stats::setNames("", "(none)"), f),
                           selected = guess %||% "")
      )
    })

    shiny::observeEvent(input$run, {
      shiny::req(rv$maf, input$tsb)
      if (!require_pkgs("maftools", "Heterogeneity")) return(NULL)
      if (!nzchar(input$vaf_col %||% "")) {
        shiny::showNotification(
          "Pick a VAF column first — this analysis is entirely based on allele frequencies.",
          type = "error", duration = 12)
        return(NULL)
      }
      het <- with_progress_notify(
        wes_heterogeneity(rv$maf, sample = input$tsb, vaf_col = input$vaf_col),
        message = "Clustering variants by allele frequency...")
      if (is.null(het)) return(NULL)
      res$het    <- het
      res$sample <- input$tsb
      mark_done(rv, "wes_hetero")
      log_step(log_rv, "WES heterogeneity",
               params = list(sample = input$tsb, vafCol = input$vaf_col),
               code = sprintf('het <- maftools::inferHeterogeneity(maf, tsb = "%s", vafCol = "%s")',
                              input$tsb, input$vaf_col))
    })

    het_df <- shiny::reactive({
      shiny::req(res$het)
      as.data.frame(res$het$clusterData)
    })

    output$summary <- shiny::renderUI({
      if (is.null(rv$maf)) return(wes_no_maf())
      if (is.null(res$het)) {
        return(wes_prompt("Pick a sample and a VAF column, then click <b>Infer clones</b>.",
                          "选择样本与 VAF 列，然后点击<b>推断克隆结构</b>。"))
      }
      d <- tryCatch(het_df(), error = function(e) NULL)
      math <- tryCatch(round(unique(d$MATH)[1], 1), error = function(e) NA_real_)
      n_cl <- tryCatch(length(unique(d$cluster)), error = function(e) NA_integer_)
      bslib::layout_columns(
        col_widths = c(4, 4, 4),
        stat_tile(i18n("Sample", "样本"), res$sample),
        stat_tile(i18n("Clusters", "克隆数"), if (is.na(n_cl)) "-" else n_cl),
        stat_tile("MATH", if (is.na(math)) "-" else math)
      )
    })

    draw_het <- with_text_boost(function() {
      shiny::req(res$het)
      maftools::plotClusters(clusters = res$het, tsb = res$sample)
    })
    output$plot <- render_base_plot(draw_het)
    register_figure_download(output, input, "plot", draw_het,
                             "wes_heterogeneity", width = 10, height = 7)

    output$tbl_slot <- shiny::renderUI({
      if (is.null(res$het)) return(wes_no_maf())
      if (has_pkg("DT")) DT::dataTableOutput(ns("tbl"))
      else shiny::verbatimTextOutput(ns("tbl_txt"))
    })
    view <- shiny::reactive({
      d <- het_df()
      num <- vapply(d, is.numeric, logical(1))
      d[num] <- lapply(d[num], function(x) round(x, 3))
      d
    })
    if (has_pkg("DT")) {
      output$tbl <- DT::renderDataTable(
        DT::datatable(view(), rownames = FALSE, filter = "top",
                      options = list(pageLength = 15, scrollX = TRUE)))
    } else {
      output$tbl_txt <- shiny::renderPrint(utils::head(view(), 20))
    }
  })
}
