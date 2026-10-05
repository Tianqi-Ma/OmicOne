#' WES module 12: Tumour heterogeneity
#'
#' Cluster one sample's variants by their allele frequency. Distinct VAF modes
#' correspond to distinct clones; the MATH score summarises the spread of the
#' VAF distribution in one number, reported here as a percentile of the
#' cohort's own MATH scores rather than against fixed cut-offs.
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
            them (a Gaussian mixture, mclust). Mutations present in every tumour
            cell sit at a high VAF; mutations acquired later by a subset sit
            lower.",
      zh = "用高斯混合模型（mclust）把某个样本的突变按其支持读段占比分组。存在于所有肿瘤细胞中的突变 VAF 较高；后来才被一部分细胞获得的突变则较低。"),
    why  = list(
      en = "A tumour is rarely one clone. Knowing which mutations are clonal
            (present everywhere, and therefore a better therapeutic target)
            versus subclonal changes what you would treat. High heterogeneity
            has been associated with worse outcome in several cancers.",
      zh = "肿瘤很少只有一个克隆。区分克隆性突变（无处不在，因而是更好的治疗靶点）与亚克隆突变，会改变治疗决策。在多种癌症中，高异质性与更差的预后相关。"),
    how  = list(
      en = "This step needs allele fractions: a <b>VAF column</b>, or
            <code>t_ref_count</code>/<code>t_alt_count</code> from which maftools
            computes them (choose <i>auto</i>). Purity and copy number move VAF:
            without a copy-number segment file, read the clusters descriptively.
            With fewer than about 10 variants the mixture fit is unstable.",
      zh = "该步骤需要等位基因频率：一个 <b>VAF 列</b>，或 <code>t_ref_count</code>/<code>t_alt_count</code>（选择<i>自动</i>，maftools 会据此计算）。肿瘤纯度与拷贝数会改变 VAF：没有拷贝数分段文件时，这些聚类只作描述性参考。变异少于约 10 个时，混合模型的拟合不稳定。"),
    read = list(
      en = "Bottom panel: the density of the sample's VAFs, with each variant
            coloured by cluster; top panels: one boxplot per cluster. A clonal
            heterozygous mutation in a diploid region sits near half the tumour
            purity (0.5 only in a pure tumour); clusters well below that are
            subclones. Outliers are drawn but not counted as clusters. MATH =
            100 × 1.4826 × MAD / median VAF; the app reports where this sample
            falls in the cohort's MATH distribution.",
      zh = "下方面板：该样本 VAF 的密度曲线，每个变异按聚类着色；上方面板：每个聚类一个箱线图。二倍体区域中的克隆性杂合突变位于肿瘤纯度的一半附近（只有纯肿瘤才是 0.5）；明显低于此值的聚类是亚克隆。离群点会画出，但不计为聚类。MATH = 100 × 1.4826 × MAD / VAF 中位数；本应用报告该样本在本队列 MATH 分布中的百分位。"),
    example = list(
      en = "TCGA-AB-2972 in the LAML demo: 19 variants, two clusters plus one
               outlier, the main cluster near VAF 0.45.",
      zh = "LAML 演示数据中的 TCGA-AB-2972：19 个变异，两个聚类外加一个离群点，主聚类位于 VAF 0.45 附近。")
  )
  controls <- shiny::tagList(
    shiny::uiOutput(ns("sample_ui")),
    shiny::uiOutput(ns("vaf_ui")),
    run_button(ns("run"), "Infer clones", "推断克隆结构")
  )
  step_container(id = id, 
    title     = list(en = "Heterogeneity", zh = "肿瘤异质性"),
    subtitle  = list(en = "Clonal structure from allele fractions (MATH score).",
                     zh = "从等位基因频率推断克隆结构（MATH 分数）。"),
    explainer = explainer,
    controls  = controls,
    summary   = shiny::uiOutput(ns("summary")),
    preview   = shiny::tagList(
      shiny::uiOutput(ns("insight")),
      bslib::navset_card_tab(
        bslib::nav_panel(i18n("Clusters", "克隆聚类"),
                         preview_plot_ui(ns("plot"), download = TRUE,
                                         guide = list(en = "The chosen sample's VAF clusters will be drawn here.",
                                                      zh = "运行后，这里将绘制所选样本的 VAF 聚类图。"),
                                         caption = list(en = "Bottom: VAF density, variants coloured by cluster; top: one boxplot per cluster. Clonal heterozygous variants sit near purity / 2.",
                                                        zh = "下方：VAF 密度，变异按聚类着色；上方：每个聚类一个箱线图。克隆性杂合变异位于纯度的一半附近。"))),
        bslib::nav_panel(i18n("Variants", "变异明细"), shiny::uiOutput(ns("tbl_slot")))
      ))
  )
}

#' @rdname mod_wes_hetero
#' @keywords internal
mod_wes_hetero_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    ns  <- session$ns
    res <- step_results(rv, "wes", het = NULL, sample = NULL, math = NULL, vaf = NULL)

    output$sample_ui <- shiny::renderUI({
      shiny::req(rv$maf)
      s <- wes_samples(rv$maf)
      wes_col_select(ns, "tsb",
                     label = list(en = "Sample", zh = "样本"),
                     tip = list(en = "Heterogeneity is inferred one sample at a time.",
                                zh = "异质性是逐个样本推断的。"),
                     choices = s, selected = if (length(s)) s[1] else NULL, selectize = TRUE)
    })

    output$vaf_ui <- shiny::renderUI({
      shiny::req(rv$maf)
      f <- wes_fields(rv$maf)
      vc <- wes_vaf_choice(f)
      wes_col_select(ns, "vaf_col",
                     label = list(en = "VAF column", zh = "VAF 列"),
                     tip = list(en = "The tumour variant allele fraction. Auto = t_vaf, or t_alt_count / (t_ref_count + t_alt_count).",
                                zh = "肿瘤变异等位基因频率。自动＝t_vaf，或 t_alt_count / (t_ref_count + t_alt_count)。"),
                     choices = f, selected = vc$selected, none = vc$none)
    })

    shiny::observeEvent(input$run, {
      shiny::req(rv$maf, input$tsb)
      if (!require_pkgs(c("maftools", "mclust"), "Heterogeneity")) return(NULL)
      sample <- input$tsb
      vaf <- if (nzchar(input$vaf_col %||% "")) input$vaf_col else NULL
      if (is.null(vaf) && !wes_vaf_auto(wes_fields(rv$maf))) {
        wes_notify("Pick a VAF column first — this MAF has neither t_vaf nor t_ref_count / t_alt_count, and the analysis is entirely based on allele fractions.",
                   "请先选择 VAF 列——该 MAF 既没有 t_vaf，也没有 t_ref_count / t_alt_count，而本分析完全基于等位基因频率。",
                   duration = 12)
        return(NULL)
      }
      maf <- rv$maf
      het <- with_progress_notify(wes_heterogeneity(maf, sample = sample, vaf_col = vaf),
                                  message = "Clustering variants by allele frequency...")
      if (is.null(het)) return(NULL)
      math <- tryCatch(wes_math_scores(maf, vaf), error = function(e) NULL)
      res$het    <- het
      res$sample <- sample
      res$math   <- math
      res$vaf    <- vaf
      mark_done(rv, "wes_hetero")
      log_step(log_rv, "WES heterogeneity",
               params = list(sample = sample, vafCol = vaf %||% "(auto)"),
               code = c(wes_code("maftools::inferHeterogeneity",
                                 list(maf = quote(maf), tsb = sample, vafCol = vaf), assign = "het"),
                        wes_code("maftools::plotClusters", list(clusters = quote(het), tsb = sample)),
                        "table(het$clusterData$cluster)   # 'outlier' / 'CN_altered' are not clones",
                        wes_code("maftools::math.score", list(maf = quote(maf), vafCol = vaf),
                                 assign = "math"),
                        sprintf("100 * mean(math$MATH <= math$MATH[math$Tumor_Sample_Barcode == %s])   # percentile",
                                deparse(sample))))
    })

    het_df <- shiny::reactive({
      shiny::req(res$het)
      as.data.frame(res$het$clusterData)
    })

    stats <- shiny::reactive({
      d <- het_df()
      pc <- if (!is.null(res$math)) wes_math_percentile(res$math, res$sample)
            else list(math = NA_real_, percentile = NA_real_, n = 0L)
      cl <- as.character(d$cluster)
      list(n_var = nrow(d), clusters = wes_clone_count(cl),
           outliers = sum(cl %in% "outlier"), cn = sum(cl %in% "CN_altered"),
           math = pc$math, pct = pc$percentile, n_math = pc$n)
    })

    output$summary <- shiny::renderUI({
      if (is.null(rv$maf)) return(wes_no_maf())
      if (is.null(res$het)) {
        return(wes_prompt("Pick a sample and a VAF column, then click <b>Infer clones</b>.",
                          "选择样本与 VAF 列，然后点击<b>推断克隆结构</b>。"))
      }
      s <- stats()
      shiny::tagList(
        stat_tile(i18n("Sample", "样本"), res$sample),
        stat_tile(i18n("Variants", "变异数"), s$n_var),
        stat_tile(i18n("Clusters", "克隆数"), s$clusters),
        stat_tile(i18n("MATH (cohort pct.)", "MATH（队列百分位）"),
                  if (is.na(s$math)) "-" else sprintf("%.1f (%.0f%%)", s$math, s$pct))
      )
    })

    output$insight <- shiny::renderUI({
      if (is.null(res$het)) return(NULL)
      s <- stats()
      few_en <- if (s$n_var < 10) sprintf(" Only %d variants: the mixture fit is unstable, treat the cluster count as tentative.", s$n_var) else ""
      few_zh <- if (s$n_var < 10) sprintf("只有 %d 个变异：混合模型拟合不稳定，克隆数仅供参考。", s$n_var) else ""
      math_en <- if (is.na(s$math)) " No MATH score (math.score needs at least 5 variants with VAF ≥ 0.075)."
                 else sprintf(" MATH %.1f is at the %.0fth percentile of the %d scored samples in this cohort.",
                              s$math, s$pct, s$n_math)
      math_zh <- if (is.na(s$math)) "没有 MATH 分数（math.score 需要至少 5 个 VAF ≥ 0.075 的变异）。"
                 else sprintf("MATH %.1f 位于本队列 %d 个可评分样本的第 %.0f 百分位。",
                              s$math, s$n_math, s$pct)
      insight_bar(
        sprintf("Sample <b>%s</b>: %d variants in <b>%d</b> VAF cluster(s) (%d outlier(s)%s not counted).%s Without purity and copy-number data, read the clusters descriptively.%s",
                res$sample, s$n_var, s$clusters, s$outliers,
                if (s$cn) sprintf(", %d copy-number-altered", s$cn) else "", math_en, few_en),
        sprintf("样本 <b>%s</b>：%d 个变异分为 <b>%d</b> 个 VAF 聚类（%d 个离群点%s不计入）。%s没有纯度和拷贝数数据时，这些聚类只作描述性解读。%s",
                res$sample, s$n_var, s$clusters, s$outliers,
                if (s$cn) sprintf("、%d 个拷贝数改变的变异", s$cn) else "", math_zh, few_zh))
    })

    draw_het <- function() {
      shiny::req(res$het)
      wes_hetero_gg(res$het, res$sample)
    }
    render_step_plot(output, input, "plot", draw_het, name = "wes_heterogeneity", width = 9, height = 5.5)

    view <- shiny::reactive({
      d <- het_df()
      # inferHeterogeneity's own MATH uses every variant; the app reports
      # math.score's (VAF >= 0.075), so only one MATH number is on screen
      d <- d[, setdiff(names(d), c("MATH", "MedianAbsoluteDeviation")), drop = FALSE]
      num <- vapply(d, is.numeric, logical(1))
      d[num] <- lapply(d[num], function(x) round(x, 3))
      d
    })
    tb <- wes_table(ns, "tbl", view, "wes_heterogeneity",
                    has_maf = function() !is.null(rv$maf),
                    ready = function() !is.null(res$het))
    output$tbl_slot <- tb$slot
    output$tbl <- tb$table
    output$tbl_dl <- tb$download
  })
}
