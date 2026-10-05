#' WES module 6: Lollipop plot / protein domains
#'
#' One gene at a time: where along the protein its mutations land, drawn against
#' the domain structure. Recurrent positions stack into tall lollipops.
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_wes_lolli
NULL

#' @rdname mod_wes_lolli
#' @keywords internal
mod_wes_lolli_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = list(en = "Lollipop / domains", zh = "Lollipop / 结构域"),
    what = list(
      en = "Every non-synonymous mutation in one gene, positioned along the
            protein sequence and drawn over its annotated domains. Stack height
            = how many times that exact change was seen.",
      zh = "某一个基因的全部非同义突变，按其在蛋白序列上的位置绘制，并叠加已注释的结构域。棒棒糖的高度＝该改变出现的次数。"),
    why  = list(
      en = "Where a mutation lands tells you what it does. Oncogenes cluster at
            one or two hotspots (a gain of function); tumour suppressors scatter
            truncating mutations across the whole length (a loss of function).",
      zh = "突变落在哪里决定了它的作用。癌基因会聚集在一两个热点上（功能获得）；抑癌基因则在全长范围内散布截短突变（功能丧失）。"),
    how  = list(
      en = "Pick a gene from the dropdown — it is sorted by mutation frequency.
            The <b>protein change column</b> is auto-detected. The protein is
            drawn on one RefSeq <b>transcript</b> (maftools picks the longest);
            your MAF was annotated on its own transcript, and if the two differ,
            positions shift or fall past the protein's end — pick the matching
            transcript here.",
      zh = "从下拉框选择基因——已按突变频率排序。<b>蛋白改变列</b>会自动识别。蛋白按某一条 RefSeq <b>转录本</b>绘制（maftools 默认取最长的）；你的 MAF 是按它自己的转录本注释的，两者不一致时位点会偏移甚至超出蛋白末端——请在此选择对应的转录本。"),
    read = list(
      en = "The x-axis is the protein sequence, boxes are annotated domains, and
            each lollipop is one observed change — height counts how often it
            was seen, colour its consequence. With <b>label recurrent
            positions</b> on, only positions hit at least twice are labelled. One
            tall lollipop inside a domain is a candidate hotspot (oncogene-style);
            truncating mutations along the whole length suggest loss of function
            (tumour-suppressor-style).",
      zh = "横轴是蛋白序列，方框为注释结构域，每个棒棒糖代表一种观察到的改变——高度是出现次数，颜色是后果。勾选<b>标注复发位点</b>时，只标注至少出现两次的位点。结构域内一根孤立的高棒棒糖是候选热点（癌基因风格）；沿全长散布的截短突变提示功能缺失（抑癌基因风格）。"),
    example = list(
      en = "<code>DNMT3A</code> in TCGA LAML: 27 of its 51 mutations with a
               protein position sit at R882, in the methyltransferase domain;
               only three other positions are hit twice.",
      zh = "TCGA LAML 中的 <code>DNMT3A</code>：51 个有蛋白位置的突变中有 27 个位于甲基转移酶结构域的 R882；另外只有三个位点各出现两次。")
  )
  controls <- shiny::tagList(
    shiny::uiOutput(ns("gene_ui")),
    shiny::uiOutput(ns("aa_ui")),
    shiny::uiOutput(ns("tx_ui")),
    shiny::checkboxInput(ns("show_rate"),
                         i18n("Show mutation rate", "显示突变率"), value = TRUE),
    shiny::checkboxInput(ns("label_pos"),
                         i18n("Label recurrent positions (seen ≥ 2 times)",
                              "标注复发位点（出现 ≥ 2 次）"), value = TRUE),
    run_button(ns("run"), "Draw lollipop", "绘制 Lollipop")
  )
  step_container(id = id, 
    title     = list(en = "Lollipop / domains", zh = "Lollipop / 结构域"),
    subtitle  = list(en = "Where on the protein the mutations land.",
                     zh = "突变落在蛋白的哪个位置。"),
    explainer = explainer,
    controls  = controls,
    summary   = shiny::uiOutput(ns("summary")),
    preview   = shiny::tagList(
      shiny::uiOutput(ns("insight")),
      preview_plot_ui(ns("plot"), download = TRUE,
                      guide = list(en = "The chosen gene's mutation map over its protein domains will be drawn here.",
                                   zh = "运行后，这里将绘制所选基因在蛋白结构域上的突变分布图。"),
                      caption = list(en = "Lollipop height = times that amino-acid change was seen; colour = consequence; boxes = protein domains of the transcript drawn.",
                                     zh = "棒棒糖高度＝该氨基酸改变出现的次数；颜色＝后果；方框＝所绘转录本的蛋白结构域。")))
  )
}

#' @rdname mod_wes_lolli
#' @keywords internal
mod_wes_lolli_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    ns  <- session$ns
    cfg <- step_result(rv, "wes")

    output$gene_ui <- shiny::renderUI({
      shiny::req(rv$maf)
      g <- wes_genes(rv$maf, n = 300)
      wes_col_select(ns, "gene",
                     label = list(en = "Gene", zh = "基因"),
                     tip = list(en = "Sorted by how many samples carry a mutation in it.",
                                zh = "按携带该基因突变的样本数排序。"),
                     choices = g, selected = if (length(g)) g[1] else NULL, selectize = TRUE)
    })

    output$aa_ui <- shiny::renderUI({
      shiny::req(rv$maf)
      f <- wes_fields(rv$maf)
      wes_col_select(ns, "aa_col",
                     label = list(en = "Protein change column", zh = "蛋白改变列"),
                     tip = list(en = "The MAF column with the amino-acid change (e.g. HGVSp_Short, AAChange). Auto-detected.",
                                zh = "MAF 中存放氨基酸改变的列（如 HGVSp_Short、AAChange），会自动识别。"),
                     choices = f, selected = wes_guess_aa_col(rv$maf) %||% "", none = "(auto)")
    })

    output$tx_ui <- shiny::renderUI({
      shiny::req(rv$maf, input$gene)
      tx <- wes_gene_transcripts(input$gene)
      if (!nrow(tx)) return(NULL)
      ch <- stats::setNames(tx$refseq.ID,
                            sprintf("%s / %s (%d aa)", tx$refseq.ID, tx$protein.ID,
                                    as.integer(tx$aa.length)))
      wes_col_select(ns, "tx",
                     label = list(en = "Transcript", zh = "转录本"),
                     tip = list(en = "Protein model to draw (refSeqID). Default: maftools' choice, the longest. Match it to the transcript your MAF was annotated on.",
                                zh = "绘制所用的蛋白模型（refSeqID）。默认：maftools 的选择，即最长的转录本。请与 MAF 注释所用的转录本保持一致。"),
                     choices = ch, none = "(longest — maftools default)")
    })

    aa_for <- function(maf, aa) {
      if (!is.null(aa)) return(aa)
      f <- wes_fields(maf)
      hit <- intersect(c("HGVSp_Short", "Protein_Change", "AAChange"), f)
      if (length(hit)) hit[1] else NULL
    }

    draw_with <- function(maf, c0) {
      wes_lollipop_gg(maf, c0$gene, aa_col = aa_for(maf, c0$aa), tx = c0$tx, label_at = c0$label_at,
                      show_rate = c0$rate)
    }

    shiny::observeEvent(input$run, {
      shiny::req(rv$maf, input$gene)
      if (!require_pkgs("maftools", "Lollipop plot")) return(NULL)
      maf <- rv$maf
      gene <- input$gene
      aa <- if (nzchar(input$aa_col %||% "")) input$aa_col else NULL
      tx <- if (nzchar(input$tx %||% "")) input$tx else NULL
      rp <- tryCatch(wes_recurrent_positions(maf, gene, aa_for(maf, aa)),
                     error = function(e) list(positions = numeric(0), all = numeric(0)))
      label <- isTRUE(input$label_pos)
      txs <- wes_gene_transcripts(gene)
      len <- if (nrow(txs)) {
        if (is.null(tx)) max(txs$aa.length) else txs$aa.length[txs$refseq.ID == tx][1]
      } else NA_real_
      c0 <- list(gene = gene, aa = aa, tx = tx, rate = isTRUE(input$show_rate),
                 label = label, label_at = if (label) rp$positions else numeric(0),
                 n_pos = length(rp$all), n_recurrent = length(rp$positions),
                 beyond = if (is.finite(len)) sum(rp$all > len) else NA_integer_,
                 aa_length = len)
      ok <- with_progress_notify(wes_dry_run(function() draw_plot_object(draw_with(maf, c0))),
                                 message = "Drawing the lollipop plot...")
      if (is.null(ok)) return(NULL)
      cfg(c0)
      mark_done(rv, "wes_lolli")
      log_step(log_rv, "WES lollipop",
               params = list(gene = gene, AACol = aa %||% "(auto)", refSeqID = tx %||% "(longest)",
                             labelPos = if (length(c0$label_at))
                               paste(c0$label_at, collapse = ", ") else "(none)"),
               code = wes_code("maftools::lollipopPlot",
                               list(maf = quote(maf), gene = gene, AACol = aa, refSeqID = tx,
                                    labelPos = if (length(c0$label_at)) c0$label_at else NULL,
                                    showMutationRate = c0$rate, repel = TRUE)))
    })

    stats <- shiny::reactive({
      c0 <- cfg()
      shiny::req(rv$maf, c0)
      fr <- wes_gene_freq(rv$maf, c0$gene)
      tot <- tryCatch({
        gs <- as.data.frame(maftools::getGeneSummary(rv$maf))
        gs$total[gs$Hugo_Symbol == c0$gene][1]
      }, error = function(e) NA_real_)
      list(mut = if (nrow(fr)) fr$mutated[1] else NA_integer_,
           pct = if (nrow(fr)) fr$pct[1] else NA_real_,
           n = wes_overview(rv$maf)$samples, total = tot)
    })

    output$summary <- shiny::renderUI({
      if (is.null(rv$maf)) return(wes_no_maf())
      c0 <- cfg()
      if (is.null(c0)) {
        return(wes_prompt("Pick a gene and click <b>Draw lollipop</b>.",
                          "选择基因后点击<b>绘制 Lollipop</b>。"))
      }
      s <- stats()
      shiny::tagList(
        stat_tile(i18n("Gene", "基因"), c0$gene),
        stat_tile(i18n("Mutated samples", "突变样本数"), wes_fmt(s$mut)),
        stat_tile(i18n("Cohort frequency", "队列频率"),
                  if (is.na(s$pct)) "-" else sprintf("%.1f%%", s$pct)),
        stat_tile(i18n("Recurrent positions", "复发位点数"), c0$n_recurrent)
      )
    })

    output$insight <- shiny::renderUI({
      c0 <- cfg()
      if (is.null(rv$maf) || is.null(c0)) return(NULL)
      s <- stats()
      if (is.na(s$mut)) return(NULL)
      tx_en <- ""
      tx_zh <- ""
      if (isTRUE(c0$beyond > 0)) {
        tx_en <- sprintf(" %d mutation(s) sit past the drawn protein's end (%d aa): the MAF was probably annotated on another transcript — pick it in <b>Transcript</b>.",
                         c0$beyond, as.integer(c0$aa_length))
        tx_zh <- sprintf("有 %d 个突变位于所绘蛋白末端（%d aa）之后：MAF 可能是按另一条转录本注释的——请在<b>转录本</b>中选择它。",
                         c0$beyond, as.integer(c0$aa_length))
      }
      insight_bar(
        sprintf("<b>%s</b> is mutated in <b>%s</b> of %s samples (%.1f%%; %s non-synonymous variants); %d position(s) are hit at least twice.%s",
                c0$gene, wes_fmt(s$mut), wes_fmt(s$n), s$pct, wes_fmt(s$total),
                c0$n_recurrent, tx_en),
        sprintf("在 %s 个样本中，有 <b>%s</b> 个样本携带 <b>%s</b> 突变（%.1f%%；共 %s 个非同义变异）；%d 个位点至少出现两次。%s",
                wes_fmt(s$n), wes_fmt(s$mut), c0$gene, s$pct, wes_fmt(s$total),
                c0$n_recurrent, tx_zh))
    })

    draw_lolli <- function() {
      c0 <- cfg()
      shiny::req(rv$maf, c0)
      draw_with(rv$maf, c0)
    }
    render_step_plot(output, input, "plot", draw_lolli, name = "wes_lollipop", width = 12, height = 6)
  })
}
