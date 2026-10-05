#' Module: Export & reproducibility
#'
#' Save the processed object for later use and download a runnable R script that
#' reproduces every step you performed, built from the reproducibility log.
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_export
NULL

#' Packages the logged code calls (`pkg::fun`)
#'
#' Base packages attached by default are left out; they need no library().
#' @param entries Log entries (each with a `code` character vector).
#' @return Sorted unique package names, Seurat / SeuratObject first.
#' @keywords internal
log_code_packages <- function(entries) {
  code <- unlist(lapply(entries %||% list(), function(e) e$code), use.names = FALSE)
  if (!length(code)) return(character(0))
  hits <- unlist(regmatches(code, gregexpr("[A-Za-z][A-Za-z0-9.]*(?=:::?)", code, perl = TRUE)))
  base <- c("base", "stats", "utils", "methods", "graphics", "grDevices", "datasets", "tools")
  pk <- sort(unique(setdiff(hits, base)))
  first <- intersect(c("Seurat", "SeuratObject"), pk)
  c(first, setdiff(pk, first))
}

#' The R session the analysis ran in, as lines of text
#'
#' Versions of the packages the log calls (scop included when it is installed)
#' and the full sessionInfo().
#' @param pkgs Package names whose versions are listed.
#' @return Character vector.
#' @keywords internal
session_info_lines <- function(pkgs = character(0)) {
  pkgs <- unique(c("OmicOne", "Seurat", "SeuratObject", "scop", pkgs))
  ver <- vapply(pkgs, function(p) {
    tryCatch(as.character(utils::packageVersion(p)), error = function(e) "not installed")
  }, character(1))
  c(sprintf("%s %s", R.version.string, R.version$platform),
    sprintf("%s: %s", pkgs, ver),
    "",
    utils::capture.output(utils::sessionInfo()))
}

#' The reproducibility script text for a set of log entries
#'
#' The header defines `input_path`, loads every package the logged code calls
#' (each call is also namespaced) and records the app's sessionInfo() as
#' comments, so a reader can see which versions produced the results.
#' @param entries Log entries. @param app App name for the title line.
#' @param session Lines from [session_info_lines()] (comment block).
#' @return One string.
#' @keywords internal
export_script_text <- function(entries, app = "OmicOne",
                               session = session_info_lines(log_code_packages(entries))) {
  pkgs <- log_code_packages(entries)
  header <- c(
    sprintf("# %s reproducibility script", app),
    paste0("# Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "#",
    "# IMPORTANT: set the input path below to your own data before running.",
    "input_path <- \"PATH/TO/YOUR/DATA\"  # <-- edit this",
    "",
    "# Packages the steps below call (every call is also written as pkg::fun).",
    if (length(pkgs)) sprintf("library(%s)", pkgs) else "# (no package calls recorded)",
    "",
    "# Session the analysis ran in (sessionInfo() when this script was exported):",
    paste0("#   ", session),
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

#' A copy of the object whose assays are all classic (v3) Seurat assays
#'
#' SeuratDisk predates Seurat 5's `Assay5`; its support for the layered
#' assay is not documented, so assays are joined and converted first.
#' @param obj A Seurat object.
#' @keywords internal
export_as_v3 <- function(obj) {
  for (a in obj_assays(obj)) {
    if (methods::is(obj[[a]], "Assay5")) {
      obj[[a]] <- SeuratObject::JoinLayers(obj[[a]])
      obj[[a]] <- methods::as(obj[[a]], "Assay")
    }
  }
  obj
}

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
            Python/Scanpy (needs SeuratDisk; Seurat 5 assays are converted to the
            classic format first). The R script loads the packages it uses, lists
            every step in order and records the session's package versions.",
      zh = "选择 <b>RDS</b> 以便在 R/Seurat 中重新载入对象，或选择 <b>.h5ad</b> 用于 Python/Scanpy（需要 SeuratDisk；Seurat 5 的 assay 会先转换为经典格式）。R 脚本会加载所用的包、按顺序列出每一步，并记录当前会话的包版本。"),
    read = list(
      en = "The step list shows what the script will replay: each step with its
            parameters and code, from the log of the steps you ran.",
      zh = "步骤列表即脚本将重放的内容：依据你运行过的步骤的日志，列出每一步及其参数和代码。"),
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

    output$summary <- shiny::renderUI({
      entries <- log_entries_for(log_rv(), "sc")
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
      explain_scene("export",
                 paste0("No plot for this step. Use the buttons on the left to download ",
                             "your processed object and the reproducibility script. ",
                             "Remember to set input/output paths when you re-run the script."),
                 paste0("这一步没有图表。请使用左侧的按钮下载你处理后的对象和可复现脚本。",
                        "重新运行脚本时，记得设置输入/输出路径。"))
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
          msg <- tryCatch({
            tmp <- tempfile(fileext = ".h5Seurat")
            SeuratDisk::SaveH5Seurat(export_as_v3(obj), filename = tmp, overwrite = TRUE)
            SeuratDisk::Convert(tmp, dest = "h5ad", overwrite = TRUE)
            file.copy(sub("\\.h5Seurat$", ".h5ad", tmp), file, overwrite = TRUE)
            NULL
          }, error = function(e) conditionMessage(e))
          if (!is.null(msg)) {
            writeLines(c("The .h5ad export failed:", msg,
                         "Download the RDS instead and convert it in R, e.g. with",
                         "zellkonverter::writeH5AD() on Seurat::as.SingleCellExperiment(obj)."),
                       file)
            shiny::showNotification(
              i18n(paste("The .h5ad export failed:", msg,
                         "The downloaded file explains the error; use RDS instead."),
                   paste("导出 .h5ad 失败：", msg, "下载的文件中给出了错误说明；请改用 RDS。")),
              type = "error", duration = 15)
            return(invisible(NULL))
          }
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
        writeLines(export_script_text(log_entries_for(log_rv(), "sc")), file)
        mark_done(rv, "export")
      }
    )
  })
}
