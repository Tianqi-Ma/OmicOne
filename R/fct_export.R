#' Shared "export current data" download handlers
#'
#' Wires the top-bar export menu (available on every step) to the working data
#' of the active pipeline. Single-cell: the Seurat object, the cell metadata,
#' the counts matrix and the embeddings; plus the clinical cohort whenever one
#' is loaded. The buttons are rendered into `output$export_items`, so the menu
#' only offers what the current pipeline can export. Another omics plugs in its
#' own buttons through `omics_items` (OmicOne passes the WES set from
#' fct_wes.R) and registers its handlers itself.
#'
#' @param input,output,session Shiny server context.
#' @param rv Shared reactive hub (uses `rv$obj`, `rv$clinical`, `rv$omics`).
#' @param omics_items Named list of zero-argument functions returning the menu
#'   buttons for a non-single-cell omics, keyed by omics.
#' @keywords internal
register_exports <- function(input, output, session, rv, omics_items = list()) {

  stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

  output$export_items <- shiny::renderUI({
    clin <- if (!is.null(rv$clinical))
      export_button("dl_clin", "Clinical cohort (.csv)", "临床队列 (.csv)")
    om <- rv$omics %||% "sc"
    if (!identical(om, "sc") && is.function(omics_items[[om]])) {
      return(shiny::tagList(omics_items[[om]](), clin))
    }
    shiny::tagList(
      export_button("dl_rds",    "Object (.rds)",        "对象 (.rds)"),
      export_button("dl_meta",   "Cell metadata (.csv)", "细胞元数据 (.csv)"),
      export_button("dl_matrix", "Counts matrix (.rds)", "表达矩阵 (.rds)"),
      export_button("dl_embed",  "Embeddings (.csv)",    "降维坐标 (.csv)"),
      clin)
  })

  need_obj <- function() {
    if (is.null(rv$obj)) {
      shiny::showNotification("Load data first, then export.", type = "warning")
      return(FALSE)
    }
    TRUE
  }

  # 1. Full Seurat object (.rds)
  output$dl_rds <- shiny::downloadHandler(
    filename = function() paste0("omicone_object_", stamp(), ".rds"),
    content = function(file) {
      if (!need_obj()) { saveRDS(NULL, file); return() }
      saveRDS(rv$obj, file)
    }
  )

  # 2. Cell metadata (.csv)
  output$dl_meta <- shiny::downloadHandler(
    filename = function() paste0("omicone_metadata_", stamp(), ".csv"),
    content = function(file) {
      if (!need_obj()) { utils::write.csv(data.frame(), file); return() }
      utils::write.csv(obj_meta(rv$obj), file, row.names = TRUE)
    }
  )

  # 3. Expression matrix — sparse counts as .rds (keeps gene/cell names, compact)
  output$dl_matrix <- shiny::downloadHandler(
    filename = function() paste0("omicone_counts_", stamp(), ".rds"),
    content = function(file) {
      if (!need_obj()) { saveRDS(NULL, file); return() }
      m <- tryCatch(SeuratObject::LayerData(rv$obj, layer = "counts"),
                    error = function(e) tryCatch(SeuratObject::GetAssayData(rv$obj, slot = "counts"),
                                                 error = function(e2) NULL))
      if (is.null(m)) { shiny::showNotification("No counts matrix found.", type = "error"); saveRDS(NULL, file); return() }
      saveRDS(m, file)
    }
  )

  # 4. Dimensional-reduction embeddings (.csv), all reductions side by side
  output$dl_embed <- shiny::downloadHandler(
    filename = function() paste0("omicone_embeddings_", stamp(), ".csv"),
    content = function(file) {
      if (!need_obj()) { utils::write.csv(data.frame(), file); return() }
      reds <- obj_reductions(rv$obj)
      if (!length(reds)) {
        shiny::showNotification("No reductions yet (run PCA/UMAP first).", type = "warning")
        utils::write.csv(data.frame(), file); return()
      }
      parts <- lapply(reds, function(r) {
        e <- tryCatch(SeuratObject::Embeddings(rv$obj, reduction = r), error = function(x) NULL)
        if (is.null(e)) return(NULL)
        colnames(e) <- paste0(r, "_", seq_len(ncol(e)))
        as.data.frame(e)
      })
      parts <- Filter(Negate(is.null), parts)
      out <- do.call(cbind, parts)
      utils::write.csv(out, file, row.names = TRUE)
    }
  )

  # --- Clinical cohort (shared across omics) ---------------------------------
  output$dl_clin <- shiny::downloadHandler(
    filename = function() paste0("omicone_clinical_", stamp(), ".csv"),
    content = function(file) {
      utils::write.csv(rv$clinical %||% data.frame(), file, row.names = FALSE)
    }
  )
}

#' One full-width button in the top-bar export menu
#' @param id Download output id. @param en,zh Label.
#' @keywords internal
export_button <- function(id, en, zh) {
  shiny::downloadButton(id, i18n(en, zh), class = "btn-sm w-100 mb-1")
}
