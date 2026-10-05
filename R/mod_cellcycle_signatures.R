#' Module: Cell cycle & signature scoring
#'
#' Two related sub-actions on the working object:
#'  (a) Cell-cycle scoring -- adds S.Score, G2M.Score and a Phase call
#'      (Seurat::CellCycleScoring with Seurat's 2019 gene lists).
#'  (b) Signature scoring -- scores built-in or custom gene sets per cell
#'      (UCell or Seurat::AddModuleScore), after a gene-coverage check.
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_cellcycle_signatures
NULL

#' Built-in example signature gene sets (small, illustrative, human symbols)
#'
#' Short hand-picked lists for a first look. They are not curated signatures:
#' for publication use a published set, e.g. MSigDB Hallmark
#' (HALLMARK_E2F_TARGETS, HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,
#' HALLMARK_HYPOXIA, HALLMARK_INFLAMMATORY_RESPONSE), pasted as a custom set.
#' @keywords internal
cellcycle_example_sets <- function() {
  list(
    proliferation = c("MKI67", "TOP2A", "PCNA", "CCNB1", "CCNB2", "CDK1",
                      "BIRC5", "AURKB", "CENPF"),
    EMT           = c("VIM", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1",
                      "FN1", "CDH2", "MMP2", "MMP9"),
    hypoxia       = c("HIF1A", "VEGFA", "CA9", "SLC2A1", "LDHA", "PGK1",
                      "ALDOA", "ENO1", "BNIP3"),
    inflammation  = c("IL6", "IL1B", "TNF", "CXCL8", "CCL2", "NFKB1",
                      "STAT1", "IRF1", "PTGS2")
  )
}

#' Parse "SetName: GENE1, GENE2" lines into a named list of gene sets
#' @param txt Text-area value.
#' @return Named list (names made safe for metadata columns).
#' @keywords internal
parse_custom_sets <- function(txt) {
  if (is.null(txt) || !nzchar(trimws(txt))) return(list())
  lines <- strsplit(txt, "\n", fixed = TRUE)[[1]]
  lines <- trimws(lines)
  lines <- lines[nzchar(lines) & grepl(":", lines, fixed = TRUE)]
  if (!length(lines)) return(list())
  nm <- trimws(sub(":.*$", "", lines))
  genes <- lapply(sub("^[^:]*:", "", lines), parse_genes)
  keep <- nzchar(nm) & lengths(genes) > 0
  stats::setNames(genes[keep], safe_set_names(nm[keep]))
}

#' @rdname mod_cellcycle_signatures
#' @keywords internal
mod_cellcycle_signatures_ui <- function(id) {
  ns <- shiny::NS(id)
  sets <- cellcycle_example_sets()
  labels <- sprintf("%s (illustrative, %d genes)", names(sets), lengths(sets))
  explainer <- explainer_card(
    title = list(en = "Cell cycle & signatures", zh = "细胞周期与信号评分"),
    what = list(
      en = "Score each cell's cell-cycle phase and the activity of gene sets such as
            proliferation, EMT or hypoxia.",
      zh = "为每个细胞评定细胞周期时相，以及增殖、EMT、缺氧等基因集的活性。"),
    why  = list(
      en = "Cell-cycle differences can dominate clustering and be mistaken for
            biology; a signature score turns a gene list into one comparable
            per-cell value.",
      zh = "细胞周期差异可能主导聚类并被误认为生物学差异；信号评分把一份基因列表变成一个可比较的单细胞数值。"),
    how  = list(
      en = "<b>Score cell cycle</b> uses Seurat's 2019 S / G2M lists (case-converted
            for mouse). For signatures tick built-in sets and/or type
            <code>SetName: GENE1, GENE2, ...</code> per line. Every set needs at
            least 3 genes present in the data. <b>UCell</b> is rank-based and
            robust to depth; <b>AddModuleScore</b> is the Seurat default.",
      zh = "<b>细胞周期评分</b>使用 Seurat 2019 版 S / G2M 基因列表（小鼠数据做大小写转换）。信号评分可勾选内置基因集，和/或按每行 <code>集合名: GENE1, GENE2, ...</code> 输入。每个基因集在数据中至少要有 3 个基因。<b>UCell</b> 基于排名，对测序深度稳健；<b>AddModuleScore</b> 为 Seurat 默认。"),
    read = list(
      en = "Each cell gets one score per set, painted on the embedding; Phase is a
            categorical call. A cluster dominated by G2M/S scores is a cycling
            state, not a separate cell type.",
      zh = "每个细胞对每个基因集得到一个分数，绘制在嵌入图上；Phase 为分类判定。被 G2M/S 分数主导的簇是增殖状态，而不是独立的细胞类型。"),
    example = list(
      en = "A tumour cluster high on proliferation and hypoxia and mostly in G2M
            points to a growing, oxygen-starved niche.",
      zh = "某肿瘤簇在增殖与缺氧上得分高且大多处于 G2M，提示一个正在生长且缺氧的微环境。")
  )
  controls <- shiny::tagList(
    shiny::tags$h6(i18n("1. Cell cycle", "1. 细胞周期")),
    label_with_help("Cell-cycle scoring",
                    "Seurat::CellCycleScoring with Seurat::cc.genes.updated.2019.",
                    label_zh = "细胞周期评分",
                    tip_zh = "使用 Seurat::cc.genes.updated.2019 运行 Seurat::CellCycleScoring。"),
    run_button(ns("run_cc"), "Score cell cycle", "细胞周期评分"),
    shiny::tags$hr(),
    shiny::tags$h6(i18n("2. Signature scoring", "2. 信号评分")),
    label_with_help("Built-in example sets",
                    "Short illustrative lists (human symbols, case-converted for mouse). For publication paste an MSigDB Hallmark set as a custom set.",
                    label_zh = "内置示例基因集",
                    tip_zh = "简短的示例列表（人类基因符号，小鼠数据做大小写转换）。用于发表时，请把 MSigDB Hallmark 基因集作为自定义集合粘贴进来。"),
    shiny::checkboxGroupInput(ns("builtin"), NULL,
                              choices = stats::setNames(names(sets), labels)),
    label_with_help("Custom gene sets",
                    "One set per line, format: SetName: GENE1, GENE2, GENE3",
                    label_zh = "自定义基因集",
                    tip_zh = "每行一个，格式：集合名: GENE1, GENE2, GENE3"),
    shiny::textAreaInput(ns("custom"), NULL, value = "", rows = 4,
                         placeholder = "MySet: CD3D, CD3E, TRAC"),
    label_with_help("Method",
                    "UCell = rank-based, robust to depth (needs the UCell package). AddModuleScore = Seurat default.",
                    label_zh = "方法",
                    tip_zh = "UCell = 基于排名、对测序深度稳健（需要 UCell 包）。AddModuleScore = Seurat 默认。"),
    shiny::selectInput(ns("method"), NULL,
                       c("UCell" = "UCell", "AddModuleScore" = "AddModuleScore")),
    run_button(ns("run_sig"), "Score signatures", "信号评分"),
    shiny::tags$hr(),
    label_with_help("Preview feature",
                    "A computed score, or Phase, to display on the embedding.",
                    label_zh = "预览特征",
                    tip_zh = "选择要在嵌入图上显示的评分或 Phase。"),
    shiny::uiOutput(ns("preview_ui"))
  )
  step_container(title = list(en = "Cell cycle & signatures", zh = "细胞周期与信号评分"),
                 subtitle = list(en = "Score cell-cycle phase and gene-set activity per cell.",
                                 zh = "为每个细胞评定细胞周期时相与基因集活性。"),
                 explainer = explainer, controls = controls,
                 summary = shiny::uiOutput(ns("summary")),
                 preview = shiny::tagList(
                   shiny::uiOutput(ns("insight")),
                   preview_plot_ui(ns("preview"), download = TRUE,
                     guide = list(en = "Per-cell scores will be drawn here.",
                                  zh = "运行后，这里将绘制每个细胞的评分。"),
                     caption = list(en = "The chosen score (colour) or Phase (category) per cell on the embedding.",
                                    zh = "所选评分（颜色）或 Phase（类别）在嵌入图上的每细胞分布。"))))
}

#' @rdname mod_cellcycle_signatures
#' @keywords internal
mod_cellcycle_signatures_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    res <- step_results(rv, "sc", cc_done = FALSE, sig_done = FALSE, cc_cov = NULL,
                        species = NULL, sig_cov = NULL, method = NULL, last = NULL,
                        run_id = 0L, converted = FALSE)
    seen_run <- 0L

    # (a) Cell-cycle scoring ---------------------------------------------------
    shiny::observeEvent(input$run_cc, {
      shiny::req(rv$obj)
      if (!require_pkgs("Seurat", "Cell-cycle scoring")) return(NULL)
      species <- guess_species(rv$obj)
      out <- with_progress_notify({
        sc_cellcycle(rv$obj, species = species)
      }, message = "Scoring cell cycle (S / G2M / Phase)...")
      if (is.null(out)) return(NULL)
      rv$obj <- out$obj
      res$cc_done <- TRUE
      res$cc_cov  <- out$coverage
      res$species <- species
      res$last    <- "Phase"
      res$run_id  <- res$run_id + 1L
      mark_done(rv, "cellcycle")
      log_step(log_rv, "Cell-cycle scoring",
               params = list(gene_lists = "Seurat::cc.genes.updated.2019", species = species,
                             S_genes_found = sprintf("%d of %d", out$coverage$n_found[1],
                                                     out$coverage$n_input[1]),
                             G2M_genes_found = sprintf("%d of %d", out$coverage$n_found[2],
                                                       out$coverage$n_input[2])),
               code = cellcycle_log_code(species))
      shiny::showNotification(
        i18n("Cell-cycle scoring done: S.Score, G2M.Score and Phase added.",
             "细胞周期评分完成：已添加 S.Score、G2M.Score 和 Phase。"),
        type = "message")
    })

    # (b) Signature scoring ----------------------------------------------------
    shiny::observeEvent(input$run_sig, {
      shiny::req(rv$obj)
      method <- input$method
      builtin <- input$builtin
      custom <- parse_custom_sets(input$custom)
      species <- guess_species(rv$obj)
      chosen <- cellcycle_example_sets()[builtin %||% character(0)]
      converted <- identical(species, "mouse") && length(chosen) > 0
      if (converted) chosen <- lapply(chosen, mouse_case)
      sets <- utils::modifyList(chosen, custom)
      if (length(sets) == 0) {
        shiny::showNotification(
          i18n("Select a built-in set or enter a custom gene set first.",
               "请先选择内置基因集或输入自定义基因集。"),
          type = "warning")
        return(NULL)
      }
      cov <- geneset_coverage(sets, rownames(rv$obj))
      low <- cov[cov$n_found < 3, , drop = FALSE]
      if (nrow(low)) {
        txt <- paste(sprintf("%s (%d of %d)", low$set, low$n_found, low$n_input), collapse = ", ")
        shiny::showNotification(
          i18n(sprintf("Each gene set needs at least 3 genes in the data. Too few found: %s.", txt),
               sprintf("每个基因集在数据中至少需要 3 个基因。以下基因集找到的基因太少：%s。", txt)),
          type = "error", duration = 12)
        return(NULL)
      }
      pkgs <- if (identical(method, "UCell")) c("Seurat", "UCell") else "Seurat"
      if (!require_pkgs(pkgs, "Signature scoring")) return(NULL)
      features <- stats::setNames(cov$found, cov$set)
      obj <- with_progress_notify({
        sc_modulescore(rv$obj, features = features, method = method)
      }, message = "Scoring gene signatures...")
      if (is.null(obj)) return(NULL)
      rv$obj <- obj
      res$sig_done  <- TRUE
      res$sig_cov   <- cov
      res$method    <- method
      res$converted <- converted
      res$last      <- modulescore_cols(cov$set, method)[1]
      res$run_id    <- res$run_id + 1L
      mark_done(rv, "cellcycle")
      log_step(log_rv, "Signature scoring",
               params = list(method = method,
                             coverage = paste(sprintf("%s %d/%d", cov$set, cov$n_found, cov$n_input),
                                              collapse = "; "),
                             mouse_case_converted = converted),
               code = modulescore_log_code(features, method))
      if (converted) {
        shiny::showNotification(
          i18n("Mouse data: built-in human gene sets were converted by case (MKI67 -> Mki67). Check the coverage above the plot.",
               "小鼠数据：内置的人类基因集已做大小写转换（MKI67 -> Mki67）。请查看图上方的覆盖情况。"),
          type = "warning", duration = 10)
      }
      shiny::showNotification(
        i18n(sprintf("Scored %d signature set(s).", length(features)),
             sprintf("已评分 %d 个信号基因集。", length(features))),
        type = "message")
    })

    # Preview feature choices: Phase / S / G2M and the score columns. After a
    # run the new column is selected; otherwise the user's choice is kept.
    output$preview_ui <- shiny::renderUI({
      cols <- obj_meta_cols(rv$obj)
      score_cols <- cols[cols %in% c("Phase", "S.Score", "G2M.Score") |
                           grepl("_UCell$|_AMS$", cols)]
      if (length(score_cols) == 0) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("Run scoring to enable the preview.",
                               "运行评分以启用预览。")))
      }
      sel <- keep_selected(shiny::isolate(input$feature), score_cols)
      if (res$run_id != seen_run && !is.null(res$last) && res$last %in% score_cols) {
        sel <- res$last
      }
      seen_run <<- res$run_id
      shiny::selectInput(session$ns("feature"), NULL, choices = score_cols, selected = sel)
    })

    output$summary <- shiny::renderUI({
      if (!isTRUE(res$cc_done) && !isTRUE(res$sig_done)) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("Score the cell cycle and/or gene signatures.",
                               "对细胞周期和/或基因信号进行评分。")))
      }
      shiny::tagList(
        stat_tile(i18n("Cell cycle", "细胞周期"),
                  if (isTRUE(res$cc_done)) i18n("Scored", "已评分") else i18n("No", "否")),
        stat_tile(i18n("Signatures", "信号集"),
                  if (isTRUE(res$sig_done)) nrow(res$sig_cov) else 0L),
        stat_tile(i18n("Method", "方法"), res$method %||% "-")
      )
    })

    output$insight <- shiny::renderUI({
      parts_en <- character(0)
      parts_zh <- character(0)
      if (isTRUE(res$cc_done)) {
        cc <- res$cc_cov
        conv_en <- if (identical(res$species, "mouse")) ", case-converted for mouse" else ""
        conv_zh <- if (identical(res$species, "mouse")) "，已为小鼠做大小写转换" else ""
        parts_en <- c(parts_en, sprintf(
          "Cell cycle: S genes %d of %d and G2M genes %d of %d found (Seurat 2019 lists%s).",
          cc$n_found[1], cc$n_input[1], cc$n_found[2], cc$n_input[2], conv_en))
        parts_zh <- c(parts_zh, sprintf(
          "细胞周期：S 基因找到 %d / %d 个，G2M 基因找到 %d / %d 个（Seurat 2019 列表%s）。",
          cc$n_found[1], cc$n_input[1], cc$n_found[2], cc$n_input[2], conv_zh))
      }
      if (isTRUE(res$sig_done)) {
        sc <- res$sig_cov
        cov <- paste(sprintf("%s %d/%d", sc$set, sc$n_found, sc$n_input), collapse = ", ")
        parts_en <- c(parts_en, sprintf(
          "%s scores, genes found per set: %s.%s Built-in sets are illustrative, not curated signatures.",
          res$method, cov,
          if (isTRUE(res$converted)) " Built-in human sets were case-converted for mouse." else ""))
        parts_zh <- c(parts_zh, sprintf(
          "%s 评分，各基因集找到的基因数：%s。%s内置基因集仅为示例，并非经整理的标准签名。",
          res$method, cov,
          if (isTRUE(res$converted)) "内置人类基因集已为小鼠做大小写转换。" else ""))
      }
      if (!length(parts_en)) return(NULL)
      insight_bar(paste(parts_en, collapse = " "), paste(parts_zh, collapse = ""))
    })

    render_step_plot(output, input, "preview", function() {
      obj <- rv$obj
      shiny::req(obj)
      feature <- input$feature
      shiny::req(feature, feature %in% obj_meta_cols(obj))
      if (identical(feature, "Phase")) {
        sc_dimplot(obj, group_by = "Phase")
      } else {
        sc_featureplot(obj, features = feature)
      }
    }, name = "cellcycle_signatures")
  })
}
