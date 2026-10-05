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
      en = "The base-change spectrum reflects the mutational processes and the
            sample handling behind the calls. VAF separates clonal (~50% in a
            pure diploid tumour) from subclonal variants. Rainfall was designed
            to show <i>kataegis</i> — localised hypermutation, a tight cluster
            low on the plot — which is defined on whole genomes: an exome
            samples only ~1–2% of the genome, so a WES rainfall rarely has the
            density to show it.",
      zh = "碱基替换谱反映了变异背后的突变过程与样本处理过程。VAF 能区分克隆性变异（纯二倍体肿瘤中约 50%）与亚克隆变异。Rainfall 原本用于展示 <i>kataegis</i>——局部超突变，表现为图上贴近底部的密集簇——它是在全基因组上定义的：外显子组只覆盖约 1–2% 的基因组，WES 的 rainfall 很少有足够的密度显示它。"),
    how  = list(
      en = "Turn <b>include synonymous</b> off to look only at coding impact.
            VAF and heterogeneity need allele fractions: a <b>VAF column</b>, or
            <code>t_ref_count</code>/<code>t_alt_count</code>, from which
            maftools computes it. Population frequencies (ExAC, gnomAD) are not
            VAFs and are never offered. Set the <b>reference build</b> to the
            MAF's coordinates so the rainfall chromosomes have the right lengths.",
      zh = "关闭<b>包含同义突变</b>可只看编码影响。VAF 与异质性分析需要等位基因频率：一个 <b>VAF 列</b>，或 <code>t_ref_count</code>/<code>t_alt_count</code>（maftools 会据此计算）。人群频率（ExAC、gnomAD）不是 VAF，不会出现在候选中。请把<b>参考基因组</b>设为 MAF 坐标所用的版本，rainfall 的染色体长度才会正确。"),
    read = list(
      en = "<b>TiTv</b>: the boxplots spread each sample's six base changes; the
            bars summarise the cohort. Exome and genome data have different
            expected ratios, and FFPE deamination inflates C>T while oxidative
            damage (8-oxoG) inflates C>A, so judge the spectrum against a
            comparable cohort rather than a fixed ratio. <b>VAF</b>: one cloud
            per gene; clonal heterozygous mutations sit near half the purity,
            subclones trail lower. <b>Rainfall</b>: most variants lie far apart
            (high on the y-axis).",
      zh = "<b>TiTv</b>：箱线图展开每个样本的六类碱基替换，柱条汇总整个队列。外显子组与全基因组的预期比例不同，FFPE 脱氨基会抬高 C>T，氧化损伤（8-oxoG）会抬高 C>A，因此请与可比队列对照，而不是套用固定比例。<b>VAF</b>：每个基因一朵云；克隆性杂合突变位于纯度一半附近，亚克隆拖在低处。<b>Rainfall</b>：多数变异彼此相距很远（纵轴高处）。"),
    example = list(
      en = "A C>T-dominated spectrum in a skin tumour fits UV damage; a C>A excess
               in a lung tumour fits tobacco, but the same excess in an FFPE or
               long-stored library can be 8-oxoG artefact — the signature step
               (SBS4 vs SBS45) tells them apart better.",
      zh = "皮肤肿瘤中以 C>T 为主的谱符合紫外损伤；肺肿瘤中 C>A 偏多符合烟草，但 FFPE 或长期保存文库中同样的偏多可能是 8-oxoG 伪影——突变特征步骤（SBS4 与 SBS45）更能区分两者。")
  )
  controls <- shiny::tagList(
    shiny::checkboxInput(ns("use_syn"),
                         i18n("Include synonymous variants", "包含同义突变"), value = TRUE),
    shiny::uiOutput(ns("vaf_ui")),
    shiny::uiOutput(ns("sample_ui")),
    shiny::uiOutput(ns("build_ui")),
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
        bslib::nav_panel("TiTv", preview_plot_ui(ns("titv"), download = TRUE,
                                                 guide = list(en = "Transition/transversion spectra for every sample will be drawn here.",
                                                              zh = "运行后，这里将绘制每个样本的转换/颠换图谱。"),
                                                 caption = list(en = "Boxplots: per-sample share of the six base changes; bars: cohort totals.",
                                                                zh = "箱线图：每样本六类碱基替换的占比；柱条：队列汇总。"))),
        bslib::nav_panel("VAF", preview_plot_ui(ns("vaf"), download = TRUE,
                                                guide = list(en = "Allele-fraction distributions of the most mutated genes will be drawn here.",
                                                             zh = "运行后，这里将绘制高频突变基因的等位基因频率分布。"),
                                                caption = list(en = "One box per gene: VAF of its mutations across samples.",
                                                               zh = "每个基因一个箱：其突变在各样本中的 VAF。"))),
        bslib::nav_panel(i18n("Rainfall", "Rainfall"),
                         preview_plot_ui(ns("rain"), download = TRUE,
                                         guide = list(en = "The chosen sample's variants along the genome will be drawn here.",
                                                      zh = "运行后，这里将绘制所选样本的变异在基因组上的分布。"),
                                         caption = list(en = "x = genomic position, y = log10 distance to the previous variant; colour = base change. Kataegis is defined on whole genomes.",
                                                        zh = "横轴＝基因组位置，纵轴＝与前一个变异距离的 log10；颜色＝碱基替换。kataegis 是在全基因组上定义的。")))
      ))
  )
}

#' @rdname mod_wes_titv
#' @keywords internal
mod_wes_titv_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    ns  <- session$ns
    res <- step_results(rv, "wes", titv = NULL, cfg = NULL)

    output$vaf_ui <- shiny::renderUI({
      shiny::req(rv$maf)
      f <- wes_fields(rv$maf)
      vc <- wes_vaf_choice(f)
      wes_col_select(ns, "vaf_col",
                     label = list(en = "VAF column", zh = "VAF 列"),
                     tip = list(en = "Tumour variant allele fraction. Auto = maftools uses t_vaf, or computes t_alt_count / (t_ref_count + t_alt_count).",
                                zh = "肿瘤变异等位基因频率。自动＝maftools 使用 t_vaf，或按 t_alt_count / (t_ref_count + t_alt_count) 计算。"),
                     choices = f, selected = vc$selected, none = vc$none)
    })

    output$sample_ui <- shiny::renderUI({
      shiny::req(rv$maf)
      s <- wes_samples(rv$maf)
      wes_col_select(ns, "tsb",
                     label = list(en = "Sample for the rainfall plot", zh = "Rainfall 的样本"),
                     tip = list(en = "Rainfall is per-sample: pick the one you want to inspect.",
                                zh = "Rainfall 是逐样本的：选择要查看的样本。"),
                     choices = s, selected = if (length(s)) s[1] else NULL, selectize = TRUE)
    })

    output$build_ui <- shiny::renderUI({
      shiny::req(rv$maf)
      guess <- wes_guess_build(rv$maf)
      wes_col_select(ns, "build",
                     label = list(en = "Reference build", zh = "参考基因组"),
                     tip = list(en = "Chromosome lengths for the rainfall x-axis. Pre-selected from the MAF's NCBI_Build column.",
                                zh = "Rainfall 横轴的染色体长度。已按 MAF 的 NCBI_Build 列预选。"),
                     choices = c("hg19" = "hg19", "hg38" = "hg38"),
                     selected = if (is.na(guess)) "hg19" else guess)
    })

    shiny::observeEvent(input$run, {
      shiny::req(rv$maf)
      if (!require_pkgs("maftools", "Mutation spectra")) return(NULL)
      use_syn <- isTRUE(input$use_syn)
      vaf <- if (nzchar(input$vaf_col %||% "")) input$vaf_col else NULL
      fields <- wes_fields(rv$maf)
      cfg <- list(vaf = vaf, vaf_ok = !is.null(vaf) || wes_vaf_auto(fields),
                  tsb = input$tsb, cp = isTRUE(input$changepoints),
                  build = input$build %||% "hg19", use_syn = use_syn)
      tv <- with_progress_notify(wes_titv(rv$maf, use_syn = use_syn),
                                 message = "Computing TiTv...")
      if (is.null(tv)) return(NULL)
      res$titv <- tv
      res$cfg  <- cfg
      mark_done(rv, "wes_titv")
      code <- c(wes_code("maftools::titv", list(maf = quote(maf), useSyn = use_syn,
                                                 plot = FALSE), assign = "titv_res"),
                wes_code("maftools::plotTiTv", list(res = quote(titv_res))),
                if (cfg$vaf_ok) wes_code("maftools::plotVaf", list(maf = quote(maf), vafCol = vaf)),
                if (!is.null(cfg$tsb))
                  wes_code("maftools::rainfallPlot",
                           list(maf = quote(maf), tsb = cfg$tsb, detectChangePoints = cfg$cp,
                                ref.build = cfg$build, pointSize = 0.6)))
      log_step(log_rv, "WES mutation spectra",
               params = list(useSyn = use_syn, vafCol = vaf %||% "(auto)",
                             rainfall_sample = cfg$tsb, ref.build = cfg$build),
               code = code)
    })

    stats <- shiny::reactive({
      shiny::req(res$titv)
      wes_titv_pooled(res$titv)
    })

    output$summary <- shiny::renderUI({
      if (is.null(rv$maf)) return(wes_no_maf())
      if (is.null(res$titv)) {
        return(wes_prompt("Click <b>Draw spectra</b> to compute the TiTv summary.",
                          "点击<b>绘制图谱</b>计算 TiTv 概览。"))
      }
      p <- stats()
      pick <- function(nm) if (nm %in% names(p)) sprintf("%.1f%%", p[[nm]]) else "-"
      shiny::tagList(
        stat_tile(i18n("C>T (pooled)", "C>T（合并）"), pick("C>T")),
        stat_tile(i18n("C>A (pooled)", "C>A（合并）"), pick("C>A")),
        stat_tile(i18n("Transitions (pooled)", "转换（合并）"), pick("Ti")),
        stat_tile(i18n("VAF column", "VAF 列"),
                  res$cfg$vaf %||% (if (isTRUE(res$cfg$vaf_ok)) i18n("auto", "自动")
                                    else i18n("none", "无")))
      )
    })

    output$insight <- shiny::renderUI({
      if (is.null(res$titv)) return(NULL)
      p <- stats()
      if (!all(c("C>T", "C>A") %in% names(p))) return(NULL)
      ct <- p[["C>T"]]
      ca <- p[["C>A"]]
      n_snv <- attr(p, "n_snv")
      lead_en <- if (isTRUE(ct >= 40)) "C>T leads, as in most clock-like / ageing spectra (FFPE deamination also raises C>T)"
                 else if (isTRUE(ca >= 25)) "C>A is high; tobacco, 8-oxoG oxidative damage and FFPE / library artefacts all raise it, so check the signature step before naming a cause"
                 else "no single base change dominates"
      lead_zh <- if (isTRUE(ct >= 40)) "C>T 居首，与多数时钟样/衰老谱一致（FFPE 脱氨基同样会抬高 C>T）"
                 else if (isTRUE(ca >= 25)) "C>A 偏高；烟草、8-oxoG 氧化损伤以及 FFPE/建库伪影都会抬高它，命名原因前请先看突变特征步骤"
                 else "没有单一碱基替换占主导"
      insight_bar(
        sprintf("Of %s SNVs pooled across the cohort, C>T is <b>%.1f%%</b> and C>A <b>%.1f%%</b> — %s.",
                wes_fmt(n_snv), ct, ca, lead_en),
        sprintf("全队列合并的 %s 个 SNV 中，C>T 占 <b>%.1f%%</b>，C>A 占 <b>%.1f%%</b>——%s。",
                wes_fmt(n_snv), ct, ca, lead_zh))
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
      if (!isTRUE(res$cfg$vaf_ok)) {
        stop("No allele fractions: pick a VAF column, or provide t_ref_count / t_alt_count in the MAF.")
      }
      maftools::plotVaf(maf = rv$maf, vafCol = res$cfg$vaf)
    })
    output$vaf <- render_base_plot(draw_vaf)
    register_figure_download(output, input, "vaf", draw_vaf, "wes_vaf",
                             width = 10, height = 7)

    draw_rain <- with_text_boost(function() {
      shiny::req(rv$maf, res$cfg, res$cfg$tsb)
      maftools::rainfallPlot(maf = rv$maf, tsb = res$cfg$tsb,
                             detectChangePoints = res$cfg$cp,
                             ref.build = res$cfg$build, pointSize = 0.6)
    })
    output$rain <- render_base_plot(draw_rain)
    register_figure_download(output, input, "rain", draw_rain, "wes_rainfall",
                             width = 12, height = 6)
  })
}
