#' Module: Export & reproducibility
#'
#' Save the processed object for later use and download a runnable R script that
#' reproduces every step you performed, built from the reproducibility log.
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_export
NULL

#' @rdname mod_export
#' @keywords internal
mod_export_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = list(en = "Export & reproducibility", zh = "导出与可复现"),
    what = list(
      en = "Download your processed object and a script that reproduces the whole
            analysis.",
      zh = "下载你处理后的对象，以及一个可复现整个分析流程的脚本。"),
    why  = list(
      en = "Reproducibility is the point: anyone (including future you) should be
            able to regenerate these results from the raw data and the script.",
      zh = "可复现是关键：任何人（包括未来的你）都应能凭原始数据和脚本重现这些结果。"),
    how  = list(
      en = "Choose <b>RDS</b> to reload the object in R/Seurat, or <b>.h5ad</b> for
            Python/Scanpy (best-effort, needs SeuratDisk). The R script lists every
            step in order -- you only need to set the input path where it says so.",
      zh = "选择 <b>RDS</b> 以便在 R/Seurat 中重新载入对象，或选择 <b>.h5ad</b> 用于 Python/Scanpy（尽力而为，需要 SeuratDisk）。该 R 脚本会按顺序列出每一步，你只需在提示处设置输入路径。"),
    read = list(
      en = "Every export replays the current object state. The log lists each
            step with its parameters — enough for anyone to reproduce the
            analysis from scratch.",
      zh = "每次导出都保存当前对象的状态。日志列出每一步及其参数——足以让任何人从头复现该分析。"),
    example = list(
      en = "Re-run with <code>source(\"omicone_analysis.R\")</code> after
               editing the <code>input_path</code> line at the top.",
      zh = "在编辑顶部的 <code>input_path</code> 行之后，使用 <code>source(\"omicone_analysis.R\")</code> 重新运行。")
  )
  controls <- shiny::tagList(
    label_with_help("Object format",
                    "RDS = native R/Seurat. .h5ad = AnnData for Python/Scanpy (needs SeuratDisk). Figures = guidance only.",
                    "对象格式",
                    "RDS = 原生 R/Seurat。.h5ad = 用于 Python/Scanpy 的 AnnData（需要 SeuratDisk）。Figures = 仅为说明。"),
    shiny::selectInput(ns("fmt"), NULL,
                       choices = c("RDS (.rds)"   = "rds",
                                   "AnnData (.h5ad)" = "h5ad",
                                   "Figures (note)"  = "figures"),
                       selected = "rds"),
    shiny::downloadButton(ns("download_obj"),
                          i18n("Download object", "下载对象"), class = "w-100"),
    shiny::tags$hr(),
    label_with_help("Reproducibility script",
                    "A commented R script rebuilding every step you ran, in order.",
                    "可复现脚本",
                    "一个带注释的 R 脚本，按顺序重建你运行过的每一步。"),
    shiny::downloadButton(ns("download_script"),
                          i18n("Download R script", "下载 R 脚本"), class = "w-100")
  )
  step_container(
    title     = list(en = "Export & reproducibility", zh = "导出与复现"),
    subtitle  = list(en = "Take the object, the tables, and a reproducible log with you.",
                     zh = "带走对象、表格与可复现日志。"),
    explainer = explainer,
    controls  = controls,
    summary   = shiny::uiOutput(ns("summary")),
    preview   = shiny::uiOutput(ns("preview"))
  )
}

#' @rdname mod_export
#' @keywords internal
mod_export_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {

    # Build the reproducibility R script text from the log entries.
    build_script <- function(entries) {
      header <- c(
        "# OmicOne reproducibility script",
        paste0("# Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        "#",
        "# IMPORTANT: set the input path below to your own data before running.",
        'input_path <- "PATH/TO/YOUR/DATA"  # <-- edit this',
        "library(Seurat)",
        ""
      )
      if (is.null(entries) || !length(entries)) {
        return(paste(c(header, "# No steps were recorded."), collapse = "\n"))
      }
      body <- unlist(lapply(seq_along(entries), function(i) {
        e <- entries[[i]]
        params <- if (length(e$params)) {
          paste(vapply(names(e$params), function(k)
            sprintf("%s=%s", k, paste(deparse(e$params[[k]]), collapse = "")),
            character(1)), collapse = ", ")
        } else ""
        c(sprintf("# Step %d: %s", i, e$step),
          if (nzchar(params)) sprintf("#   params: %s", params) else NULL,
          if (length(e$code)) e$code else "# (no code recorded)",
          "")
      }))
      paste(c(header, body), collapse = "\n")
    }

    output$summary <- shiny::renderUI({
      entries <- log_rv()
      if (is.null(entries) || !length(entries)) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("No steps recorded yet. Run some analysis steps first.",
                               "尚未记录任何步骤。请先运行一些分析步骤。")))
      }
      shiny::tagList(
        stat_tile(i18n("Steps performed", "已执行步骤数"), length(entries)),
        shiny::tags$ol(class = "omicone-steps",
          lapply(entries, function(e) {
            shiny::tags$li(shiny::tags$b(e$step),
                           shiny::tags$span(class = "omicone-muted",
                                            paste0("  (", e$time, ")")))
          }))
      )
    })

    output$preview <- shiny::renderUI({
      shiny::div(class = "omicone-note",
                 i18n(paste0("No plot for this step. Use the buttons on the left to download ",
                             "your processed object and the reproducibility script. ",
                             "Remember to set input/output paths when you re-run the script."),
                      paste0("这一步没有图表。请使用左侧的按钮下载你处理后的对象和可复现脚本。",
                             "重新运行脚本时，记得设置输入/输出路径。")))
    })

    output$download_obj <- shiny::downloadHandler(
      filename = function() {
        fmt <- input$fmt %||% "rds"
        ext <- switch(fmt, rds = "rds", h5ad = "h5ad", figures = "txt")
        paste0("omicone_object_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".", ext)
      },
      content = function(file) {
        fmt <- input$fmt %||% "rds"
        obj <- rv$obj
        if (is.null(obj)) {
          writeLines("No object available to export.", file)
          return(invisible(NULL))
        }
        if (fmt == "rds") {
          saveRDS(obj, file)
        } else if (fmt == "h5ad") {
          if (!require_pkgs(c("SeuratDisk", "Seurat"), "AnnData (.h5ad) export")) {
            writeLines("SeuratDisk not installed; could not write .h5ad.", file)
            return(invisible(NULL))
          }
          tmp <- tempfile(fileext = ".h5Seurat")
          SeuratDisk::SaveH5Seurat(obj, filename = tmp, overwrite = TRUE)
          SeuratDisk::Convert(tmp, dest = "h5ad", overwrite = TRUE)
          file.copy(sub("\\.h5Seurat$", ".h5ad", tmp), file, overwrite = TRUE)
        } else {
          writeLines(c(
            "Figures are downloaded individually from the Visualize step",
            "using its 'Download plot' button (PNG/PDF)."), file)
        }
        mark_done(rv, "export")
      }
    )

    output$download_script <- shiny::downloadHandler(
      filename = function()
        paste0("omicone_analysis_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".R"),
      content = function(file) {
        writeLines(build_script(log_rv()), file)
        mark_done(rv, "export")
      }
    )
  })
}
