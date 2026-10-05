#' WES module 8: Mutational signatures
#'
#' Decompose the cohort's trinucleotide mutation spectrum into de-novo
#' signatures and match them against COSMIC, which names the mutational process
#' behind each one. Before any NMF, the reference build is checked against the
#' MAF (NCBI_Build, and the reference bases of a sample of SNVs when BSgenome
#' is installed), because a wrong build still yields a plausible-looking matrix.
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
            into a small number of recurring patterns (NMF).",
      zh = "把每一个单碱基替换按照碱基改变类型及其两侧碱基分类（96 类），再用 NMF 把这个矩阵分解成少数几种反复出现的模式。"),
    why  = list(
      en = "Each mutational process leaves its own 96-category fingerprint. The
            decomposition recovers those processes — UV, tobacco, APOBEC,
            defective mismatch repair, platinum chemotherapy — from the mutations
            alone. Exomes carry far fewer SNVs than genomes, so signatures from a
            small WES cohort are noisy.",
      zh = "每一种突变过程都会留下自己独特的 96 类指纹。这个分解仅凭突变本身就能还原出这些过程——紫外、烟草、APOBEC、错配修复缺陷、铂类化疗等。外显子组的 SNV 远少于全基因组，小型 WES 队列得到的特征噪声较大。"),
    how  = list(
      en = "This needs the reference genome as a <b>BSgenome</b> package so the
            flanking bases can be looked up — install
            <code>BSgenome.Hsapiens.UCSC.hg19</code> (or hg38) first. The build is
            pre-selected from the MAF's <code>NCBI_Build</code>; a mismatch is
            flagged, and a sample of reference bases is checked before running.
            Tick <b>estimate the number of signatures</b> to see the cophenetic
            curve (slow); otherwise start at 3 and raise it only on large
            cohorts.",
      zh = "该步骤需要参考基因组的 <b>BSgenome</b> 包以查询侧翼碱基——请先安装 <code>BSgenome.Hsapiens.UCSC.hg19</code>（或 hg38）。参考基因组按 MAF 的 <code>NCBI_Build</code> 预选；不一致时会提示，运行前还会抽查一批参考碱基。勾选<b>估计特征数</b>可查看 cophenetic 曲线（较慢）；否则从 3 开始，只有队列很大时才增加。"),
    read = list(
      en = "<b>Signatures</b>: each extracted signature's 96-category fingerprint,
            titled with its closest COSMIC match. <b>Exposures</b>: how much of
            each signature every sample carries — the stacked bars should differ
            between samples, otherwise one process dominates the cohort. Cosine
            similarity ≥ 0.85 is a confident match; below that the label is
            shown in grey and is a guess. <b>Rank estimate</b>: pick the rank
            just before the cophenetic coefficient drops.",
      zh = "<b>特征谱</b>：每个提取特征的 96 类指纹，标题为最接近的 COSMIC 匹配。<b>暴露量</b>：每个样本携带各特征的比例——堆叠条应在样本间有差异，否则说明该队列由单一突变过程主导。余弦相似度 ≥ 0.85 才算可信匹配；低于此值的标签以灰色显示，仅供参考。<b>特征数估计</b>：选择 cophenetic 系数开始下降之前的那个秩。"),
    example = list(
      en = "A cohort dominated by <b>SBS4</b> is a smoking cohort; <b>SBS2/SBS13</b>
               together mean APOBEC activity; <b>SBS6/SBS15</b> point at mismatch
               repair deficiency, which matters for immunotherapy.",
      zh = "以 <b>SBS4</b> 为主的队列是吸烟队列；<b>SBS2/SBS13</b> 同时出现意味着 APOBEC 活性；<b>SBS6/SBS15</b> 指向错配修复缺陷——这与免疫治疗相关。")
  )
  controls <- shiny::tagList(
    shiny::uiOutput(ns("build_ui")),
    label_with_help("Number of signatures",
                    "How many de-novo signatures to extract. Too many on a small cohort produces noise.",
                    label_zh = "特征数量",
                    tip_zh = "提取多少个 de-novo 特征。小队列上取太多只会得到噪声。"),
    shiny::numericInput(ns("n_sig"), NULL, value = 3, min = 2, max = 8, step = 1),
    shiny::checkboxInput(ns("estimate"),
                         i18n("Also estimate the number of signatures (ranks 2–6, slow)",
                              "同时估计特征数（秩 2–6，较慢）"),
                         value = FALSE),
    run_button(ns("run"), "Extract signatures", "提取突变特征")
  )
  step_container(id = id, 
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
                                         guide = list(en = "The extracted signature fingerprints will be drawn here.",
                                                      zh = "运行后，这里将绘制提取出的特征指纹。"),
                                         caption = list(en = "One panel per signature: share of its mutations in each of the 96 trinucleotide categories; title = closest COSMIC match.",
                                                        zh = "每个面板一个特征：其突变在 96 类三核苷酸类别中的占比；标题＝最接近的 COSMIC 匹配。"))),
        bslib::nav_panel(i18n("Exposures", "暴露量"),
                         preview_plot_ui(ns("contrib"), download = TRUE,
                                         guide = list(en = "Each sample's share of every signature will be drawn here.",
                                                      zh = "运行后，这里将绘制每个样本各特征所占的比例。"),
                                         caption = list(en = "One bar per sample: fraction of its SNVs attributed to each signature.",
                                                        zh = "每个样本一根柱：其 SNV 归属于各特征的比例。"))),
        bslib::nav_panel(i18n("Rank estimate", "特征数估计"),
                         preview_plot_ui(ns("rank"), download = TRUE,
                                         guide = list(en = "Tick 'estimate the number of signatures' and run to draw the cophenetic curve here.",
                                                      zh = "勾选“估计特征数”并运行后，这里将绘制 cophenetic 曲线。"),
                                         caption = list(en = "Cophenetic correlation per NMF rank; the line marks the rank you extracted.",
                                                        zh = "每个 NMF 秩的 cophenetic 相关系数；竖线标记你提取时使用的秩。"))),
        bslib::nav_panel(i18n("COSMIC match", "COSMIC 匹配"), shiny::uiOutput(ns("tbl_slot"))),
        bslib::nav_panel(i18n("APOBEC enrichment", "APOBEC 富集"),
                         preview_plot_ui(ns("apo"), download = TRUE,
                                         guide = list(en = "APOBEC-enriched versus non-enriched samples will be compared here.",
                                                      zh = "运行后，这里将比较 APOBEC 富集与非富集样本。"),
                                         caption = list(en = "APOBEC enrichment (Roberts et al. 2013, one-sided Fisher test per sample) and the genes mutated differently between enriched and non-enriched samples.",
                                                        zh = "APOBEC 富集（Roberts 等 2013，每样本单侧 Fisher 检验），以及富集与非富集样本之间突变频率不同的基因。")))
      ))
  )
}

#' @rdname mod_wes_sig
#' @keywords internal
mod_wes_sig_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    ns  <- session$ns
    res <- step_results(rv, "wes", tnm = NULL, sig = NULL, cmp = NULL, n = NA_integer_,
                        build = NULL, db = NULL, pconstant = NULL, est = NULL,
                        check = NULL, maf_build = NA_character_)

    output$build_ui <- shiny::renderUI({
      shiny::req(rv$maf)
      guess <- wes_guess_build(rv$maf)
      wes_col_select(ns, "build",
                     label = list(en = "Reference genome", zh = "参考基因组"),
                     tip = list(en = "Must match the coordinates in your MAF (pre-selected from NCBI_Build). Needs the matching BSgenome package installed.",
                                zh = "必须与 MAF 中的坐标一致（已按 NCBI_Build 预选），并已安装对应的 BSgenome 包。"),
                     choices = c("hg19" = "hg19", "hg38" = "hg38"),
                     selected = keep_selected(shiny::isolate(input$build), c("hg19", "hg38"),
                                              if (is.na(guess)) "hg19" else guess))
    })

    shiny::observeEvent(input$run, {
      shiny::req(rv$maf)
      if (!require_pkgs("maftools", "Mutational signatures")) return(NULL)
      build <- input$build %||% "hg19"
      n <- wes_int(input$n_sig, 3, lo = 2, hi = 8)
      estimate <- isTRUE(input$estimate)
      bs <- wes_bsgenome_pkg(build)
      if (!require_pkgs(c(bs, "NMF"), "Mutational signatures")) return(NULL)
      maf <- rv$maf

      maf_build <- wes_guess_build(maf)
      if (!is.na(maf_build) && !identical(maf_build, build)) {
        wes_notify(sprintf("The MAF's NCBI_Build says %s, but %s is selected. Check the coordinates before trusting the signatures.",
                           maf_build, build),
                   sprintf("MAF 的 NCBI_Build 为 %s，但当前选择的是 %s。请先核对坐标，再采信突变特征结果。",
                           maf_build, build),
                   type = "warning", duration = 15)
      }
      chk <- tryCatch(wes_ref_check(maf, build), error = function(e) NULL)
      if (!is.null(chk) && chk$n > 0 && chk$concordance < 0.95) {
        wes_notify(sprintf("Only %.0f%% of %d checked SNVs match the %s reference base: the build is probably wrong (or the MAF's Reference_Allele is not the + strand). Stopped before extraction.",
                           100 * chk$concordance, chk$n, build),
                   sprintf("抽查的 %d 个 SNV 中只有 %.0f%% 与 %s 参考碱基一致：参考基因组很可能选错了（或 MAF 的 Reference_Allele 不是正链）。已在提取前停止。",
                           chk$n, 100 * chk$concordance, build),
                   duration = 20)
        return(NULL)
      }

      tnm <- with_progress_notify(wes_trinuc(maf, build = build),
                                  message = "Building the trinucleotide matrix...")
      if (is.null(tnm)) return(NULL)
      out <- with_progress_notify(wes_signatures(tnm, n = n),
                                  message = "Extracting signatures...")
      if (is.null(out)) return(NULL)
      est <- NULL
      if (estimate) {
        est <- with_progress_notify(wes_estimate_rank(tnm, n_try = 6, nrun = 10),
                                    message = "Estimating the number of signatures...")
      }

      res$tnm <- tnm
      res$sig <- out$sig
      res$cmp <- out$cmp
      res$db <- out$db
      res$pconstant <- out$pconstant
      res$n <- n
      res$build <- build
      res$est <- est
      res$check <- chk
      res$maf_build <- maf_build
      mark_done(rv, "wes_sig")
      tn_args <- list(maf = quote(maf), ref_genome = bs)
      if (isTRUE(attr(tnm, "prefix_added"))) {
        tn_args$prefix <- "chr"
        tn_args$add <- TRUE
      }
      log_step(log_rv, "WES mutational signatures",
               params = list(genome = build, n = n, sig_db = out$db,
                             pConstant = out$pconstant %||% "(none)",
                             ref_check = if (is.null(chk)) "skipped (no BSgenome)"
                                         else sprintf("%.1f%% of %d", 100 * chk$concordance, chk$n)),
               code = c(wes_code("maftools::trinucleotideMatrix", tn_args, assign = "tnm"),
                        paste0("if (!\"package:NMF\" %in% search()) attachNamespace(\"NMF\")",
                               "   # NMF's seed generic, see wes_signatures()"),
                        if (estimate)
                          wes_code("maftools::estimateSignatures",
                                   list(mat = quote(tnm), nMin = 2, nTry = 6, nrun = 10,
                                        parallel = wes_raw("NULL"), pConstant = est$pconstant),
                                   assign = "est"),
                        wes_code("maftools::extractSignatures",
                                 list(mat = quote(tnm), n = n, parallel = wes_raw("NULL"),
                                      pConstant = out$pconstant),
                                 assign = "sig"),
                        wes_code("maftools::compareSignatures",
                                 list(nmfRes = quote(sig), sig_db = out$db), assign = "cmp"),
                        wes_code("maftools::plotSignatures",
                                 list(nmfRes = quote(sig), sig_db = out$db)),
                        wes_code("maftools::plotSignatures",
                                 list(nmfRes = quote(sig), contributions = TRUE)),
                        wes_code("maftools::plotApobecDiff", list(tnm = quote(tnm), maf = quote(maf)))))
    })

    stats <- shiny::reactive({
      shiny::req(res$sig)
      m <- tryCatch(wes_sig_matches(res$cmp), error = function(e) NULL)
      list(matches = m, n = res$n, build = res$build,
           n_confident = if (!is.null(m)) sum(m$cosine >= 0.85) else NA_integer_)
    })

    output$summary <- shiny::renderUI({
      if (is.null(rv$maf)) return(wes_no_maf())
      if (is.null(res$sig)) {
        return(wes_prompt("Pick a reference genome and click <b>Extract signatures</b>. This step needs a BSgenome package and takes a minute.",
                          "选择参考基因组后点击<b>提取突变特征</b>。该步骤需要 BSgenome 包，耗时约一分钟。"))
      }
      s <- stats()
      m <- s$matches
      cos_tiles <- if (!is.null(m)) {
        lapply(seq_len(nrow(m)), function(i) {
          val <- sprintf("%s (%.2f)", m$best_match[i], m$cosine[i])
          if (m$cosine[i] < 0.85) val <- shiny::span(class = "omicone-muted", val)
          stat_tile(m$signature[i], val)
        })
      }
      shiny::tagList(
        stat_tile(i18n("Signatures", "特征数"), s$n),
        stat_tile(i18n("Genome (run)", "基因组（运行时）"), s$build),
        utils::head(cos_tiles, 3)
      )
    })

    output$insight <- shiny::renderUI({
      if (is.null(res$sig)) return(NULL)
      s <- stats()
      m <- s$matches
      txt <- if (!is.null(m)) paste(sprintf("%s → %s (cosine %.2f)", m$signature, m$best_match,
                                            m$cosine), collapse = "; ") else "-"
      weak <- if (!is.null(m)) sum(m$cosine < 0.85) else 0
      chk <- res$check
      chk_en <- if (is.null(chk)) " Reference bases were not checked (BSgenome not installed)."
                else sprintf(" %.1f%% of %d sampled SNVs matched the %s reference.",
                             100 * chk$concordance, chk$n, res$build)
      chk_zh <- if (is.null(chk)) "未检查参考碱基（未安装 BSgenome）。"
                else sprintf("抽查的 %d 个 SNV 中 %.1f%% 与 %s 参考一致。",
                             chk$n, 100 * chk$concordance, res$build)
      pc_en <- if (!is.null(res$pconstant)) " NMF needed pConstant = 1e-4 (empty categories)." else ""
      pc_zh <- if (!is.null(res$pconstant)) "NMF 需要 pConstant = 1e-4（存在空类别）。" else ""
      if (!is.na(res$maf_build) && !identical(res$maf_build, res$build)) {
        chk_en <- paste0(chk_en, sprintf(" <b>The MAF declares %s (NCBI_Build) but %s was used.</b>",
                                         res$maf_build, res$build))
        chk_zh <- paste0(chk_zh, sprintf("<b>MAF 声明的是 %s（NCBI_Build），但运行时使用了 %s。</b>",
                                         res$maf_build, res$build))
      }
      weak_en <- if (weak > 0) sprintf(" %d match(es) fall below cosine 0.85 and are tentative.", weak)
                 else " All matches reach cosine 0.85."
      weak_zh <- if (weak > 0) sprintf("其中 %d 个匹配的余弦相似度低于 0.85，仅供参考。", weak)
                 else "所有匹配的余弦相似度均达到 0.85。"
      insight_bar(
        sprintf("Extracted <b>%d</b> signatures (%s): %s.%s%s%s",
                s$n, res$db, txt, weak_en, chk_en, pc_en),
        sprintf("提取出 <b>%d</b> 个特征（%s）：%s。%s%s%s",
                s$n, res$db, txt, weak_zh, chk_zh, pc_zh))
    })

    draw_sig <- function() {
      shiny::req(res$sig)
      wes_signatures_gg(res$sig, tryCatch(wes_sig_matches(res$cmp), error = function(e) NULL))
    }
    render_step_plot(output, input, "sig", draw_sig, name = "wes_signatures", width = 12,
                     height = function() max(4, 1.8 + 1.9 * res$n))

    draw_contrib <- function() {
      shiny::req(res$sig)
      wes_exposures_gg(res$sig)
    }
    render_step_plot(output, input, "contrib", draw_contrib, name = "wes_signature_exposures",
                     width = 12, height = 5.5)

    draw_rank <- function() {
      shiny::req(res$sig)
      if (is.null(res$est)) {
        stop("The rank was not estimated in this run: tick 'Also estimate the number of signatures' and run again.")
      }
      wes_cophenetic_gg(res$est$res, best = res$n)
    }
    render_step_plot(output, input, "rank", draw_rank, name = "wes_signature_rank", width = 7, height = 5)

    draw_apo <- function() {
      shiny::req(res$tnm, rv$maf)
      wes_apobec_gg(res$tnm, rv$maf)
    }
    render_step_plot(output, input, "apo", draw_apo, name = "wes_apobec", width = 11, height = 8)

    view <- shiny::reactive({
      shiny::req(res$cmp)
      m <- stats()$matches
      shiny::req(m)
      m$confident <- m$cosine >= 0.85
      m
    })
    tb <- wes_table(ns, "tbl", view, "wes_cosmic_matches",
                    has_maf = function() !is.null(rv$maf),
                    ready = function() !is.null(res$cmp))
    output$tbl_slot <- tb$slot
    output$tbl <- tb$table
    output$tbl_dl <- tb$download
  })
}
