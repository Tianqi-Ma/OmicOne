#' WES module 8: Mutational signatures
#'
#' Decompose the cohort's trinucleotide mutation spectrum into de-novo
#' signatures and match them against COSMIC, which names the mutational process
#' behind each one.
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_wes_sig
NULL

#' @rdname mod_wes_sig
#' @keywords internal
mod_wes_sig_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = list(en = "Mutational signatures", zh = "突变特征"),
    what = list(
      en = "Classify every single-base substitution by its base change and the
            two bases flanking it (96 categories), then factorise that matrix
            into a small number of recurring patterns.",
      zh = "把每一个单碱基替换按照碱基改变类型及其两侧碱基分类（96 类），再把这个矩阵分解成少数几种反复出现的模式。"),
    why  = list(
      en = "Each mutational process leaves its own 96-category fingerprint. The
            decomposition recovers those processes — UV, tobacco, APOBEC,
            defective mismatch repair, platinum chemotherapy — from the mutations
            alone.",
      zh = "每一种突变过程都会留下自己独特的 96 类指纹。这个分解仅凭突变本身就能还原出这些过程——紫外、烟草、APOBEC、错配修复缺陷、铂类化疗等。"),
    how  = list(
      en = "This needs the reference genome as a <b>BSgenome</b> package so the
            flanking bases can be looked up — install
            <code>BSgenome.Hsapiens.UCSC.hg19</code> (or hg38) first; it is a
            large download. Then pick how many signatures to extract: start at 3
            and raise it only if the cohort is large.",
      zh = "该步骤需要参考基因组的 <b>BSgenome</b> 包以查询侧翼碱基——请先安装 <code>BSgenome.Hsapiens.UCSC.hg19</code>（或 hg38），这是个较大的下载。然后选择要提取的特征数：从 3 开始，只有队列很大时才增加。"),
    read = list(
      en = "Left: each extracted signature's 96-category fingerprint. Right: how
            much of each signature every sample carries — the stacked bars
            should differ between samples, otherwise one process dominates the
            cohort. The <b>COSMIC match</b> tab names each signature: SBS1 is
            ageing, SBS4 tobacco, SBS2/SBS13 APOBEC, SBS6/SBS15 mismatch-repair
            deficiency. Cosine similarity above ~0.85 is a confident match;
            below ~0.6 the label is a guess.",
      zh = "左图：每个提取特征的 96 类指纹。右图：每个样本携带各特征的比例——堆叠条应在样本间有差异，否则说明该队列由单一突变过程主导。<b>COSMIC 匹配</b>页签给出每个特征的名称：SBS1 为衰老，SBS4 为烟草，SBS2/SBS13 为 APOBEC，SBS6/SBS15 为错配修复缺陷。余弦相似度高于约 0.85 才算可信匹配；低于 0.6 时标签仅供参考。"),
    example = list(
      en = "A cohort dominated by <b>SBS4</b> is a smoking cohort; <b>SBS2/SBS13</b>
               together mean APOBEC activity; <b>SBS6/SBS15</b> point at mismatch
               repair deficiency, which matters for immunotherapy.",
      zh = "以 <b>SBS4</b> 为主的队列是吸烟队列；<b>SBS2/SBS13</b> 同时出现意味着 APOBEC 活性；<b>SBS6/SBS15</b> 指向错配修复缺陷——这与免疫治疗相关。")
  )
  controls <- shiny::tagList(
    label_with_help("Reference genome",
                    "Must match the coordinates in your MAF. Needs the matching BSgenome package installed.",
                    label_zh = "参考基因组",
                    tip_zh = "必须与 MAF 中的坐标一致，并已安装对应的 BSgenome 包。"),
    shiny::selectInput(ns("build"), NULL, c("hg19" = "hg19", "hg38" = "hg38"),
                       selected = "hg19"),
    label_with_help("Number of signatures",
                    "How many de-novo signatures to extract. Too many on a small cohort produces noise.",
                    label_zh = "特征数量",
                    tip_zh = "提取多少个 de-novo 特征。小队列上取太多只会得到噪声。"),
    shiny::numericInput(ns("n_sig"), NULL, value = 3, min = 2, max = 8, step = 1),
    run_button(ns("run"), "Extract signatures", "提取突变特征")
  )
  step_container(
    title     = list(en = "Mutational signatures", zh = "突变特征"),
    subtitle  = list(en = "Decompose the mutation spectrum into known mutational processes.",
                     zh = "把突变谱分解为已知的突变过程。"),
    explainer = explainer,
    controls  = controls,
    summary   = shiny::uiOutput(ns("summary")),
    preview   = shiny::tagList(
      shiny::uiOutput(ns("insight")),
      bslib::navset_card_tab(
      bslib::nav_panel(i18n("Signatures", "特征谱"),
                       preview_plot_ui(ns("sig"), download = TRUE,
                                       guide = list(en = "The extracted signature fingerprints and per-sample exposures will be drawn here.",
                                                    zh = "运行后，这里将绘制提取出的特征指纹与各样本的暴露量。"),
                                       caption = list(en = "Left: signature fingerprints over 96 trinucleotide categories; right: per-sample exposure.",
                                                      zh = "左：96 类三核苷酸特征指纹；右：每样本暴露量。"))),
      bslib::nav_panel(i18n("COSMIC match", "COSMIC 匹配"), shiny::uiOutput(ns("tbl_slot"))),
      bslib::nav_panel(i18n("APOBEC enrichment", "APOBEC 富集"), preview_plot_ui(ns("apo"), download = TRUE))
    ))
  )
}

#' @rdname mod_wes_sig
#' @keywords internal
mod_wes_sig_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    ns  <- session$ns
    res <- shiny::reactiveValues(tnm = NULL, sig = NULL, cmp = NULL, n = NA_integer_)

    shiny::observeEvent(input$run, {
      shiny::req(rv$maf)
      if (!require_pkgs("maftools", "Mutational signatures")) return(NULL)
      bs <- wes_bsgenome_pkg(input$build)
      if (!require_pkgs(c(bs, "NMF"), "Mutational signatures")) return(NULL)

      tnm <- with_progress_notify(wes_trinuc(rv$maf, build = input$build),
                                  message = "Building the trinucleotide matrix...")
      if (is.null(tnm)) return(NULL)
      out <- with_progress_notify(
        wes_signatures(tnm, n = max(2, as.integer(input$n_sig))),
        message = "Extracting signatures...")
      if (is.null(out)) return(NULL)

      res$tnm <- tnm; res$sig <- out$sig; res$cmp <- out$cmp
      res$n <- max(2, as.integer(input$n_sig))
      mark_done(rv, "wes_sig")
      log_step(log_rv, "WES mutational signatures",
               params = list(genome = input$build, n = input$n_sig),
               code = c(sprintf('tnm <- maftools::trinucleotideMatrix(maf, ref_genome = "%s")', bs),
                        sprintf('sig <- maftools::extractSignatures(tnm, n = %s)', input$n_sig),
                        'maftools::compareSignatures(sig, sig_db = "SBS")'))
    })

    output$summary <- shiny::renderUI({
      if (is.null(rv$maf)) return(wes_no_maf())
      if (is.null(res$sig)) {
        return(wes_prompt("Pick a reference genome and click <b>Extract signatures</b>. This step needs a BSgenome package and takes a minute.",
                          "选择参考基因组后点击<b>提取突变特征</b>。该步骤需要 BSgenome 包，耗时约一分钟。"))
      }
      best <- tryCatch({
        txt <- paste(unlist(lapply(res$cmp$best_match, unlist)), collapse = " ")
        hits <- unique(regmatches(txt, gregexpr("SBS[0-9]+[a-z]?", txt))[[1]])
        if (length(hits)) paste(utils::head(hits, 4), collapse = ", ") else "-"
      }, error = function(e) "-")
      bslib::layout_columns(
        col_widths = c(4, 4, 4),
        stat_tile(i18n("Signatures", "特征数"), res$n),
        stat_tile(i18n("Genome", "基因组"), input$build),
        stat_tile(i18n("Best COSMIC match", "最佳 COSMIC 匹配"), best)
      )
    })

    output$insight <- shiny::renderUI({
      if (is.null(res$sig)) return(NULL)
      hits <- tryCatch({
        txt <- paste(unlist(lapply(res$cmp$best_match, unlist)), collapse = " ")
        h <- unique(regmatches(txt, gregexpr("SBS[0-9]+[a-z]?", txt))[[1]])
        if (length(h)) paste(utils::head(h, 4), collapse = ", ") else "none"
      }, error = function(e) "-")
      insight_bar(
        sprintf("Extracted <b>%d</b> signatures; closest COSMIC matches: <b>%s</b>. If a name is unfamiliar, compare the 96-category fingerprints themselves, not just the labels.",
                res$n, hits),
        sprintf("提取出 <b>%d</b> 个特征；最接近的 COSMIC 匹配：<b>%s</b>。对陌生的名称，请直接比对 96 类指纹本身，而不只是看标签。",
                res$n, hits))
    })

    draw_sig <- with_text_boost(function() {
      shiny::req(res$sig)
      maftools::plotSignatures(nmfRes = res$sig, title_size = 1.0,
                               sig_db = "SBS")
    })
    output$sig <- render_base_plot(draw_sig)
    register_figure_download(output, input, "sig", draw_sig, "wes_signatures",
                             width = 11, height = function() max(4, 2 + 1.6 * res$n))

    draw_apo <- with_text_boost(function() {
      shiny::req(res$tnm)
      msg <- tryCatch({ maftools::plotApobecDiff(tnm = res$tnm, maf = rv$maf); NULL },
                      error = function(e) conditionMessage(e))
      if (!is.null(msg)) {
        if (grepl("differentially mutated genes", msg, ignore.case = TRUE)) {
          stop("No APOBEC enrichment found in this cohort — that mutational process is essentially absent here (expected for many tumour types, e.g. leukaemia).")
        }
        stop(msg)
      }
    })
    output$apo <- render_base_plot(draw_apo)
    register_figure_download(output, input, "apo", draw_apo, "wes_apobec",
                             width = 10, height = 7)

    output$tbl_slot <- shiny::renderUI({
      if (is.null(res$cmp)) return(wes_no_maf())
      if (has_pkg("DT")) DT::dataTableOutput(ns("tbl"))
      else shiny::verbatimTextOutput(ns("tbl_txt"))
    })
    # compareSignatures() returns a per-signature list whose exact shape has
    # shifted between maftools versions; flatten it defensively rather than
    # assuming field names.
    view <- shiny::reactive({
      cmp <- res$cmp; shiny::req(cmp)
      bm <- cmp$best_match
      flat <- function(x) paste(unlist(lapply(x, as.character)), collapse = " | ")
      data.frame(
        signature = names(bm) %||% paste0("Signature_", seq_along(bm)),
        match = vapply(bm, flat, character(1)),
        stringsAsFactors = FALSE, row.names = NULL
      )
    })
    if (has_pkg("DT")) {
      output$tbl <- DT::renderDataTable(
        DT::datatable(view(), rownames = FALSE,
                      options = list(pageLength = 10, scrollX = TRUE)))
    } else {
      output$tbl_txt <- shiny::renderPrint(view())
    }
  })
}
