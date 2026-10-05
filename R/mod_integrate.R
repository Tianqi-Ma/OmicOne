#' Module: Batch integration (optional)
#'
#' Align cells from different samples/batches so that shared cell types overlap
#' instead of forming separate, technically-driven clumps. Harmony is the
#' one-click default; "none" leaves the data uncorrected (and removes any
#' earlier integrated space).
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_integrate
NULL

#' @rdname mod_integrate
#' @keywords internal
mod_integrate_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = list(en = "Batch integration (optional)", zh = "批次整合（可选）"),
    what = list(
      en = "Correct for technical differences between samples/batches so that the
            same cell type from different samples lines up.",
      zh = "校正样本/批次之间的技术差异，使来自不同样本的相同细胞类型对齐。"),
    why  = list(
      en = "Without correction, cells often cluster by their sample of origin
            rather than by biology, which confounds downstream clustering.",
      zh = "若不校正，细胞常常按来源样本而非生物学聚在一起，从而干扰下游聚类。"),
    how  = list(
      en = "Pick the categorical column that identifies your batch/sample, then a
            method. <b>Harmony</b> is fast and a good default; <b>CCA/RPCA</b>
            need Seurat 5 and a LogNormalize-normalised object. Choose
            <b>none</b> for a single sample; it also removes an earlier
            integration so clustering returns to PCA.",
      zh = "选择标识批次/样本的分类列，再选择方法。<b>Harmony</b> 快速且是不错的默认；<b>CCA/RPCA</b> 需要 Seurat 5 和 LogNormalize 归一化的对象。单样本请选 <b>none</b>；它也会删除先前的整合结果，使聚类回到 PCA。"),
    read = list(
      en = "Left: the first two uncorrected PCs; right: the first two dims of the
            integrated space; colour = batch. The panel titles give the
            batch-mixing score (kNN batch entropy relative to the batch
            composition; 1 = fully mixed). Higher mixing is wanted, but
            distinct cell types must stay apart: mixing alone cannot show that
            biology survived.",
      zh = "左：未校正的前两个主成分；右：整合空间的前两维；颜色 = 批次。面板标题给出批次混合度（kNN 批次熵相对于批次构成的熵；1 = 完全混合）。混合度越高越好，但不同细胞类型必须仍然分开：仅凭混合度无法说明生物学信号得以保留。"),
    example = list(
      en = "Two samples whose cells sit apart on PCA (mixing 0.3) overlap after
               Harmony (mixing 0.9) while T and B cells stay separate.",
      zh = "两个样本的细胞在 PCA 上彼此分开（混合度 0.3），经 Harmony 后重叠（混合度 0.9），而 T 细胞与 B 细胞仍保持分离。")
  )
  controls <- shiny::tagList(
    label_with_help("Batch / sample column",
                    "A categorical metadata column (2-200 levels, no missing values) that identifies each sample or batch.",
                    "批次/样本列",
                    "标识每个样本或批次的分类元数据列（2-200 个水平，无缺失值）。"),
    shiny::uiOutput(ns("batch_ui")),
    label_with_help("Method",
                    "Harmony is fast and robust. CCA/RPCA are Seurat v5 anchor-based (IntegrateLayers). none = no correction.",
                    "方法",
                    "Harmony 快速且稳健。CCA/RPCA 是 Seurat v5 基于锚点的方法（IntegrateLayers）。none = 不做校正。"),
    shiny::selectInput(ns("method"), NULL,
                       c("none (no correction)" = "none", "Harmony" = "harmony",
                         "CCA" = "CCA", "RPCA" = "RPCA")),
    label_with_help("PCs used", "Leading PCs passed to the integration (capped at the PCs computed).",
                    "使用的主成分数", "传入整合的前若干个主成分（不超过已计算的主成分数）。"),
    shiny::numericInput(ns("dims"), NULL, value = 30, min = 2, max = 200, step = 1),
    run_button(ns("run"), "Run integration", "运行整合")
  )
  step_container(title = list(en = "Batch integration", zh = "批次整合"),
                 subtitle = list(en = "Align batches so the same cell type sits together.",
                                 zh = "对齐批次，使相同细胞类型聚在一起。"),
                 explainer = explainer, controls = controls,
                 summary = shiny::uiOutput(ns("summary")),
                 preview = shiny::tagList(
                   shiny::uiOutput(ns("insight")),
                   preview_plot_ui(ns("preview"), download = TRUE,
                     guide = list(en = "Batches on the uncorrected PCA and on the integrated space will be drawn side by side here.",
                                  zh = "运行后，这里将并排绘制批次在未校正 PCA 与整合空间上的分布。"),
                     caption = list(en = "One point = one cell, coloured by batch; panel titles give the kNN batch-mixing score.",
                                    zh = "每个点为一个细胞，按批次着色；面板标题给出 kNN 批次混合度。"))))
}

#' @rdname mod_integrate
#' @keywords internal
mod_integrate_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    res <- step_results(rv, "sc", done = FALSE, method = NULL, batch = NULL,
                        reduction = NULL, pd = NULL)

    # Only categorical columns can be a batch: never a continuous QC metric,
    # never a single-level or partly missing column.
    output$batch_ui <- shiny::renderUI({
      md <- obj_meta(rv$obj)
      if (!ncol(md)) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("Load a dataset to choose a batch column.",
                               "加载数据集以选择批次列。")))
      }
      cols <- sc_group_cols(md)
      if (!length(cols)) {
        return(shiny::div(class = "omicone-note",
                          i18n("No categorical column with 2 or more levels: this looks like a single sample, so integration can be skipped (or run with method none).",
                               "没有 2 个及以上水平的分类列：看起来是单样本，可以跳过整合（或以 none 方法运行）。")))
      }
      default <- guess_batch_col(md[, cols, drop = FALSE]) %||% cols[1]
      shiny::selectInput(session$ns("batch"), NULL, choices = cols,
                         selected = keep_selected(shiny::isolate(input$batch), cols, default))
    })

    shiny::observeEvent(input$run, {
      shiny::req(rv$obj)
      method <- input$method
      batch <- input$batch
      dims <- int_input(input$dims, 2, 500)
      if (is.na(dims)) {
        shiny::showNotification(i18n("Enter a whole number of PCs.", "请输入整数的主成分数。"),
                                type = "error")
        return(NULL)
      }
      if (method != "none" && (is.null(batch) || !batch %in% sc_group_cols(obj_meta(rv$obj)))) {
        shiny::showNotification(i18n("No usable batch column (single sample): integration skipped.",
                                     "没有可用的批次列（单样本）：已跳过整合。"),
                                type = "warning")
        return(NULL)
      }
      pkgs <- c("Seurat", if (method == "harmony") "harmony")
      if (!require_pkgs(pkgs, "Integration")) return(NULL)
      seed <- 42
      out <- with_progress_notify({
        o <- integrate_obj(rv$obj, batch = batch, method = method, dims = dims, seed = seed)
        red <- obj_misc(o, "omicone_reduction") %||% "pca"
        pd <- if (!is.null(batch) && has_reduction(o, "pca")) {
          tryCatch(integrate_plot_data(o, batch, red, dims), error = function(e) NULL)
        }
        list(obj = o, red = red, pd = pd)
      }, message = "Integrating batches...")
      if (is.null(out)) return(NULL)
      rv$obj <- out$obj
      res$done <- TRUE
      res$method <- method
      res$batch <- batch
      res$reduction <- out$red
      res$pd <- out$pd
      mark_done(rv, "integrate")
      kw <- if (method %in% c("CCA", "RPCA")) {
        integrate_k_weight(as.character(obj_meta(rv$obj)[[batch]]))
      } else 100
      log_step(log_rv, "Integration",
               params = list(method = method, batch = if (method == "none") NULL else batch,
                             dims = dims, seed = seed),
               code = integrate_log_code(method, batch, dims, seed, kw))
      shiny::showNotification(
        if (method == "none") {
          i18n("No correction: downstream steps use PCA.", "未做校正：下游步骤使用 PCA。")
        } else {
          i18n(sprintf("Integration done using %s on '%s'.", method, batch),
               sprintf("已用 %s 按“%s”完成整合。", method, batch))
        },
        type = "message")
    })

    output$summary <- shiny::renderUI({
      if (!isTRUE(res$done)) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("Choose a batch column and method, then click Run integration.",
                               "选择批次列和方法，然后点击运行整合。")))
      }
      shiny::tagList(
        stat_tile(i18n("Method", "方法"), res$method),
        stat_tile(i18n("Batch column", "批次列"), res$batch %||% "-"),
        stat_tile(i18n("Reduction for clustering", "聚类所用降维"), res$reduction)
      )
    })

    output$insight <- shiny::renderUI({
      pd <- res$pd
      if (!isTRUE(res$done)) return(NULL)
      if (is.null(pd)) {
        return(insight_bar("No batch correction was applied; Cluster and Embed default to PCA.",
                           "未做批次校正；聚类与降维可视化默认使用 PCA。"))
      }
      ent <- pd$entropy
      if (all(is.na(ent))) return(NULL)
      n_b <- length(unique(pd$df$batch))
      if (length(ent) == 1) {
        return(insight_bar(
          sprintf("Batch mixing on uncorrected PCA: %.2f across %d batches (1 = fully mixed).",
                  ent[[1]], n_b),
          sprintf("未校正 PCA 上的批次混合度：%.2f（%d 个批次；1 = 完全混合）。", ent[[1]], n_b)))
      }
      insight_bar(
        sprintf("Batch mixing %.2f on PCA vs %.2f on %s across %d batches (kNN entropy, 1 = fully mixed). Check that distinct cell types stay apart; mixing also rises when integration over-corrects.",
                ent[["pca"]], ent[[2]], names(ent)[2], n_b),
        sprintf("%d 个批次的混合度：PCA 上 %.2f，%s 上 %.2f（kNN 熵，1 = 完全混合）。请确认不同细胞类型仍然分开；整合过度时混合度同样会升高。",
                n_b, ent[["pca"]], names(ent)[2], ent[[2]]))
    })

    render_step_plot(output, input, "preview", function() {
      shiny::req(res$done, res$pd)
      integrate_plot(res$pd, res$batch %||% "batch")
    }, name = "integrate")
  })
}
