#' WES module 4: TiTv, VAF and rainfall
#'
#' Three views of the mutation spectrum rather than the gene list: the
#' transition/transversion balance across the cohort, the allele-frequency
#' distribution per gene, and the genomic distribution of one sample's variants.
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_wes_titv
NULL

#' @rdname mod_wes_titv
#' @keywords internal
mod_wes_titv_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = list(en = "TiTv / VAF / rainfall", zh = "TiTv / VAF / rainfall"),
    what = list(
      en = "<b>TiTv</b>: the balance of transitions (A↔G, C↔T) versus
            transversions. <b>VAF</b>: what fraction of reads carried the variant,
            per gene. <b>Rainfall</b>: each variant plotted along the genome
            against the distance to the previous one.",
      zh = "<b>TiTv</b>：转换（A↔G、C↔T）与颠换的比例。<b>VAF</b>：每个基因的变异等位基因频率分布。<b>Rainfall</b>：把每个变异按其在基因组上的位置、以及与前一个变异的距离画出来。"),
    why  = list(
      en = "The TiTv spectrum is a fingerprint of the underlying mutational
            process. VAF separates clonal (~50% in a pure diploid tumour) from
            subclonal variants. Rainfall reveals <i>kataegis</i> — localised
            hypermutation showing up as a tight cluster low on the plot.",
      zh = "TiTv 谱是潜在突变过程的指纹。VAF 能区分克隆性变异（纯二倍体肿瘤中约 50%）与亚克隆变异。Rainfall 可揭示 <i>kataegis</i>——局部超突变，表现为图上贴近底部的密集簇。"),
    how  = list(
      en = "Turn <b>include synonymous</b> off to look only at coding impact.
            VAF and rainfall need a <b>VAF column</b>; many MAFs do not have one,
            in which case those two tabs will say so.",
      zh = "关闭<b>包含同义突变</b>可只看编码影响。VAF 和 rainfall 需要一个 <b>VAF 列</b>；很多 MAF 没有该列，这两个页签会给出提示。"),
    read = list(
      en = "<b>TiTv</b>: the boxplots spread each sample's six base changes; the
            bars summarise the cohort. Transitions normally outnumber
            transversions ~2:1 — a flipped ratio means a mutagen or a QC
            problem. <b>VAF</b>: one cloud per gene; clonal heterozygous
            mutations sit near 0.5, subclones trail lower. <b>Rainfall</b>: most
            variants lie far apart (high on the y-axis); a tight cluster hugging
            the bottom is a kataegis hotspot.",
      zh = "<b>TiTv</b>：箱线图展开每个样本的六类碱基替换，柱条汇总整个队列。正常情况下转换约为颠换的 2 倍——比例倒挂提示诱变剂暴露或质控问题。<b>VAF</b>：每个基因一朵云；克隆性杂合突变位于 0.5 附近，亚克隆拖在低处。<b>Rainfall</b>：多数变异彼此相距很远（纵轴高处）；贴底的密集簇即 kataegis 热点。"),
    example = list(
      en = "A C>T dominated spectrum in a skin tumour points at UV damage; in a
               lung tumour a C>A excess points at tobacco.",
      zh = "皮肤肿瘤中以 C>T 为主的谱指向紫外损伤；肺肿瘤中 C>A 偏多则指向烟草。")
  )
  controls <- shiny::tagList(
    shiny::checkboxInput(ns("use_syn"),
                         i18n("Include synonymous variants", "包含同义突变"), value = TRUE),
    shiny::uiOutput(ns("vaf_ui")),
    shiny::uiOutput(ns("sample_ui")),
    shiny::checkboxInput(ns("changepoints"),
                         i18n("Detect kataegis change points", "检测 kataegis 变化点"),
                         value = TRUE),
    run_button(ns("run"), "Draw spectra", "绘制图谱")
  )
  step_container(
    title     = list(en = "TiTv / VAF / rainfall", zh = "TiTv / VAF / rainfall"),
    subtitle  = list(en = "Three spectra: base-change balance, allele fractions, genome-wide spacing.",
                     zh = "三种图谱：碱基替换平衡、等位基因频率、全基因组间距。"),
    explainer = explainer,
    controls  = controls,
    summary   = shiny::uiOutput(ns("summary")),
    preview   = shiny::tagList(
      shiny::uiOutput(ns("insight")),
      bslib::navset_card_tab(
      bslib::nav_panel("TiTv",     preview_plot_ui(ns("titv"), download = TRUE,
                                       guide = list(en = "Transition/transversion spectra for every sample will be drawn here.",
                                                    zh = "运行后，这里将绘制每个样本的转换/颠换图谱。"),
                                       caption = list(en = "Boxplots: per-sample share of the six base changes; bars: cohort totals.",
                                                      zh = "箱线图：每样本六类碱基替换的占比；柱条：队列汇总。"))),
      bslib::nav_panel("VAF",      preview_plot_ui(ns("vaf"), download = TRUE)),
      bslib::nav_panel(i18n("Rainfall", "Rainfall"), preview_plot_ui(ns("rain"), download = TRUE))
    ))
  )
}

#' @rdname mod_wes_titv
#' @keywords internal
mod_wes_titv_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    ns  <- session$ns
    res <- shiny::reactiveValues(titv = NULL, cfg = NULL)

    output$vaf_ui <- shiny::renderUI({
      shiny::req(rv$maf)
      f <- wes_fields(rv$maf)
      guess <- wes_guess_vaf_col(rv$maf)
      shiny::tagList(
        label_with_help("VAF column",
                        "Which MAF column holds the variant allele frequency. Needed by the VAF and rainfall tabs.",
                        label_zh = "VAF 列",
                        tip_zh = "MAF 中存放变异等位基因频率的列，VAF 与 rainfall 页签需要它。"),
        shiny::selectInput(ns("vaf_col"), NULL,
                           choices = c(stats::setNames("", "(none)"), f),
                           selected = guess %||% "")
      )
    })

    output$sample_ui <- shiny::renderUI({
      shiny::req(rv$maf)
      s <- wes_samples(rv$maf)
      shiny::tagList(
        label_with_help("Sample for the rainfall plot",
                        "Rainfall is per-sample: pick the one you want to inspect.",
                        label_zh = "Rainfall 的样本",
                        tip_zh = "Rainfall 是逐样本的：选择要查看的样本。"),
        shiny::selectInput(ns("tsb"), NULL, choices = s,
                           selected = if (length(s)) s[1] else NULL)
      )
    })

    shiny::observeEvent(input$run, {
      shiny::req(rv$maf)
      if (!require_pkgs("maftools", "Mutation spectra")) return(NULL)
      tv <- with_progress_notify(wes_titv(rv$maf, use_syn = isTRUE(input$use_syn)),
                                 message = "Computing TiTv...")
      if (is.null(tv)) return(NULL)
      res$titv <- tv
      res$cfg  <- list(vaf = if (nzchar(input$vaf_col %||% "")) input$vaf_col else NULL,
                       tsb = input$tsb,
                       cp  = isTRUE(input$changepoints))
      mark_done(rv, "wes_titv")
      log_step(log_rv, "WES mutation spectra",
               params = list(useSyn = input$use_syn, vafCol = input$vaf_col,
                             rainfall_sample = input$tsb),
               code = c(sprintf('titv <- maftools::titv(maf, useSyn = %s, plot = FALSE)',
                                input$use_syn),
                        'maftools::plotTiTv(titv)'))
    })

    output$summary <- shiny::renderUI({
      if (is.null(rv$maf)) return(wes_no_maf())
      if (is.null(res$titv)) {
        return(wes_prompt("Click <b>Draw spectra</b> to compute the TiTv summary.",
                          "点击<b>绘制图谱</b>计算 TiTv 概览。"))
      }
      frac <- tryCatch(as.data.frame(res$titv$fraction.contribution),
                       error = function(e) NULL)
      pick <- function(nm) {
        if (is.null(frac) || !nm %in% colnames(frac)) return("-")
        sprintf("%.1f%%", mean(frac[[nm]], na.rm = TRUE))
      }
      bslib::layout_columns(
        col_widths = c(3, 3, 3, 3),
        stat_tile("C>T", pick("C>T")),
        stat_tile("C>A", pick("C>A")),
        stat_tile("T>C", pick("T>C")),
        stat_tile(i18n("VAF column", "VAF 列"),
                  res$cfg$vaf %||% i18n("none", "无"))
      )
    })

    output$insight <- shiny::renderUI({
      tv <- res$titv
      if (is.null(tv)) return(NULL)
      frac <- tryCatch(as.data.frame(tv$fraction.contribution),
                       error = function(e) NULL)
      if (is.null(frac) || !all(c("C>T", "C>A") %in% colnames(frac))) return(NULL)
      ct <- mean(frac[["C>T"]], na.rm = TRUE)
      ca <- mean(frac[["C>A"]], na.rm = TRUE)
      verdict_en <- if (isTRUE(ca >= 25)) "a C>A share this high points at tobacco exposure"
                    else if (isTRUE(ct >= 40)) "a C>T-led spectrum is the normal ageing background"
                    else "no single base change dominates"
      verdict_zh <- if (isTRUE(ca >= 25)) "C>A 占比如此之高，提示烟草暴露"
                    else if (isTRUE(ct >= 40)) "以 C>T 为主是正常的衰老背景"
                    else "没有单一碱基替换占主导"
      insight_bar(
        sprintf("C>T accounts for <b>%.1f%%</b> of substitutions, C>A for <b>%.1f%%</b> — %s.",
                ct, ca, verdict_en),
        sprintf("C>T 占碱基替换的 <b>%.1f%%</b>，C>A 占 <b>%.1f%%</b>——%s。",
                ct, ca, verdict_zh))
    })

    draw_titv <- with_text_boost(function() {
      shiny::req(res$titv)
      maftools::plotTiTv(res = res$titv)
    })
    output$titv <- render_base_plot(draw_titv)
    register_figure_download(output, input, "titv", draw_titv, "wes_titv",
                             width = 9, height = 6)

    draw_vaf <- with_text_boost(function() {
      shiny::req(rv$maf, res$cfg)
      if (is.null(res$cfg$vaf)) {
        stop("No VAF column selected. Pick one in the control panel, or this MAF does not carry allele frequencies.")
      }
      maftools::plotVaf(maf = rv$maf, vafCol = res$cfg$vaf)
    })
    output$vaf <- render_base_plot(draw_vaf)
    register_figure_download(output, input, "vaf", draw_vaf, "wes_vaf",
                             width = 10, height = 7)

    draw_rain <- with_text_boost(function() {
      shiny::req(rv$maf, res$cfg, res$cfg$tsb)
      maftools::rainfallPlot(maf = rv$maf, tsb = res$cfg$tsb,
                             detectChangePoints = res$cfg$cp, pointSize = 0.6)
    })
    output$rain <- render_base_plot(draw_rain)
    register_figure_download(output, input, "rain", draw_rain, "wes_rainfall",
                             width = 12, height = 6)
  })
}
