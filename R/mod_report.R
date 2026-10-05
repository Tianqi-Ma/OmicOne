#' Module: Export report
#'
#' Generate a narrated, scCancer-style HTML report of the analysis: the steps
#' run, their parameters and equivalent code, and the R session, assembled from
#' the reproducibility log. No heavy compute happens here -- it only reads what
#' previous steps logged.
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_report
NULL

#' Report sections and the step keys each one covers
#'
#' Log entries carry the key of the step that wrote them (`$key`), so a
#' section matches steps by key, never by a word in the step's name.
#' @return Named list: en, zh, keys.
#' @keywords internal
report_sections <- function() {
  wes_keys <- tryCatch(vapply(steps_wes(), function(s) s$v, character(1)),
                       error = function(e) character(0))
  list(
    qc         = list(en = "Quality control",     zh = "质量控制",     keys = "qc"),
    doublet    = list(en = "Doublets",            zh = "双细胞",       keys = "doublet"),
    normalize  = list(en = "Normalization",       zh = "标准化",       keys = "normalize"),
    reduce     = list(en = "Features & PCA",      zh = "特征与 PCA",   keys = "reduce"),
    integrate  = list(en = "Integration",         zh = "整合",         keys = "integrate"),
    cluster    = list(en = "Clustering",          zh = "聚类",         keys = "cluster"),
    embed      = list(en = "Embedding",           zh = "降维可视化",   keys = "embed"),
    markers    = list(en = "Marker genes",        zh = "标志基因",     keys = "markers"),
    annotation = list(en = "Annotation",          zh = "注释",         keys = "annotate"),
    enrichment = list(en = "Enrichment / GSEA",   zh = "富集 / GSEA",  keys = "enrichment"),
    trajectory = list(en = "Trajectory, velocity & dynamics", zh = "轨迹、速率与动态",
                      keys = c("trajectory", "velocity", "dynamic")),
    signatures = list(en = "Cell cycle & signatures", zh = "细胞周期与信号", keys = "cellcycle"),
    cellcomm   = list(en = "Cell communication",  zh = "细胞通讯",     keys = "cellcomm"),
    malignancy = list(en = "Malignant / CNV",     zh = "恶性 / CNV",   keys = "malignancy"),
    survival   = list(en = "Clinical & survival", zh = "临床与生存",   keys = "clinical"),
    wes        = list(en = "Somatic mutations",   zh = "体细胞突变",   keys = wes_keys)
  )
}

#' Keep the log entries of the ticked report sections
#'
#' An entry whose step key belongs to no section (or that has no key) is
#' always kept, so nothing that ran is silently left out.
#' @param entries Log entries. @param sections Ticked section names.
#' @param defs Output of [report_sections()].
#' @keywords internal
report_filter_entries <- function(entries, sections, defs = report_sections()) {
  if (is.null(entries) || !length(entries)) return(entries)
  all_keys <- unlist(lapply(defs, `[[`, "keys"), use.names = FALSE)
  sel_keys <- unlist(lapply(defs[intersect(sections, names(defs))], `[[`, "keys"),
                     use.names = FALSE)
  Filter(function(e) {
    k <- e$key %||% ""
    k %in% sel_keys || !k %in% all_keys
  }, entries)
}

#' @rdname mod_report
#' @keywords internal
mod_report_ui <- function(id) {
  ns <- shiny::NS(id)
  secs <- report_sections()
  explainer <- explainer_card(
    title = list(en = "Export report", zh = "导出报告"),
    what = list(
      en = "Generate an HTML report of your analysis: for each step you ran, when it
            ran, its parameters and its R code, followed by the R session and
            package versions.",
      zh = "生成一份 HTML 分析报告：对你运行过的每一步，列出运行时间、参数与 R 代码，最后附上 R 会话与包版本。"),
    why  = list(
      en = "A shareable, human-readable record (scCancer-style) of what was done and
            how, for collaborators, supervisors or a methods section.",
      zh = "一份可分享、易读的记录（scCancer 风格），说明做了什么以及如何做，供合作者、导师或方法学部分使用。"),
    how  = list(
      en = "Tick the sections to include and give the report a title, then click
            <b>Download report</b>. The report is built from the log of the steps
            you ran; a step re-run replaces its earlier entry.",
      zh = "勾选要包含的章节并为报告命名，然后点击<b>下载报告</b>。报告依据你运行过的步骤的日志生成；重新运行某一步会替换它之前的记录。"),
    read = list(
      en = "Each section is one step: time, parameters and code. The report holds no
            figures; download those from each step's figure button.",
      zh = "每个章节对应一步：时间、参数与代码。报告不含图表；图表请在各步骤中通过下载按钮获取。"),
    example = list(
      en = "A report titled <i>PBMC 3k analysis</i> with QC, Clustering and
               Annotation sections, ready to attach to an email.",
      zh = "一份题为 <i>PBMC 3k analysis</i> 的报告，包含质量控制、聚类和注释章节，可直接作为邮件附件。")
  )
  controls <- shiny::tagList(
    label_with_help("Report title", "Shown as the report heading.",
                    label_zh = "报告标题", tip_zh = "作为报告的标题显示。"),
    shiny::textInput(ns("title"), NULL, value = "OmicOne analysis report"),
    label_with_help("Sections to include",
                    "Only the ticked sections are written (matched by the key of the steps you ran).",
                    label_zh = "包含的章节",
                    tip_zh = "只写入勾选的章节（按你运行过的步骤的标识匹配）。"),
    shiny::checkboxGroupInput(
      ns("sections"), NULL,
      choiceNames = unname(lapply(secs, function(x) i18n(x$en, x$zh))),
      choiceValues = names(secs),
      selected = names(secs)),
    shiny::downloadButton(ns("download_report"),
                          i18n("Download report", "下载报告"), class = "w-100")
  )
  step_container(
    title     = list(en = "Export report", zh = "导出报告"),
    explainer = explainer,
    controls  = controls,
    summary   = shiny::uiOutput(ns("summary")),
    preview   = shiny::uiOutput(ns("preview"))
  )
}

#' @rdname mod_report
#' @keywords internal
mod_report_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {

    # Single-cell steps, plus the WES steps when that section is ticked.
    report_entries <- shiny::reactive({
      all <- log_rv()
      entries <- c(log_entries_for(all, "sc"), log_entries_for(all, "wes"))
      report_filter_entries(entries, input$sections %||% character(0))
    })

    # Format one step's parameters as "k = v, k = v".
    params_text <- function(params) {
      if (!length(params)) return("")
      paste(vapply(names(params), function(k)
        sprintf("%s = %s", k, paste(deparse(params[[k]]), collapse = "")),
        character(1)), collapse = ", ")
    }

    # Assemble the report body as an R Markdown string (used by rmarkdown, or
    # rendered to a self-contained HTML fallback).
    build_rmd <- function(title, entries) {
      lines <- c(
        "---",
        sprintf('title: "%s"', gsub('"', "'", title)),
        sprintf('date: "%s"', format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        "output:",
        "  html_document:",
        "    self_contained: true",
        "    toc: true",
        "---",
        "",
        "> Generated by OmicOne from the log of the steps you ran.",
        "")
      if (is.null(entries) || !length(entries)) {
        return(paste(c(lines, "No analysis steps were recorded yet."), collapse = "\n"))
      }
      body <- unlist(lapply(seq_along(entries), function(i) {
        e <- entries[[i]]
        p <- params_text(e$params)
        c(sprintf("## %d. %s", i, e$step),
          sprintf("*Run at %s.*", e$time %||% "unknown time"),
          "",
          if (nzchar(p)) c("**Parameters:**", "", sprintf("`%s`", p), "") else NULL,
          if (length(e$code)) c("**Code:**", "", "```r", e$code, "```", "") else NULL)
      }))
      sess <- c("## R session and package versions", "", "```",
                session_info_lines(log_code_packages(entries)), "```", "")
      paste(c(lines, body, sess), collapse = "\n")
    }

    # Self-contained HTML fallback when rmarkdown/pandoc is unavailable.
    build_html <- function(title, entries) {
      esc <- function(x) {
        x <- gsub("&", "&amp;", x, fixed = TRUE)
        x <- gsub("<", "&lt;",  x, fixed = TRUE)
        gsub(">", "&gt;", x, fixed = TRUE)
      }
      css <- paste(
        "body{font-family:system-ui,Arial,sans-serif;max-width:860px;margin:2rem auto;",
        "padding:0 1rem;color:#1f2933;line-height:1.6}",
        "h1{border-bottom:2px solid #3b6ea5;padding-bottom:.3rem}",
        "h2{color:#3b6ea5;margin-top:2rem}",
        ".muted{color:#66788a;font-size:.9em}",
        "pre{background:#f4f6f8;padding:.8rem;border-radius:6px;overflow:auto}",
        "code{background:#f4f6f8;padding:.1rem .3rem;border-radius:4px}",
        collapse = "")
      head <- c(
        "<!DOCTYPE html><html><head><meta charset='utf-8'>",
        sprintf("<title>%s</title>", esc(title)),
        sprintf("<style>%s</style></head><body>", css),
        sprintf("<h1>%s</h1>", esc(title)),
        sprintf("<p class='muted'>Generated by OmicOne from the log of the steps you ran &middot; %s</p>",
                format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
      if (is.null(entries) || !length(entries)) {
        body <- "<p>No analysis steps were recorded yet.</p>"
      } else {
        body <- unlist(lapply(seq_along(entries), function(i) {
          e <- entries[[i]]
          p <- params_text(e$params)
          c(sprintf("<h2>%d. %s</h2>", i, esc(e$step)),
            sprintf("<p class='muted'>Run at %s.</p>", esc(e$time %||% "unknown time")),
            if (nzchar(p)) sprintf("<p><b>Parameters:</b> <code>%s</code></p>", esc(p)) else "",
            if (length(e$code))
              sprintf("<p><b>Code:</b></p><pre><code>%s</code></pre>",
                      esc(paste(e$code, collapse = "\n"))) else "")
        }))
      }
      sess <- c("<h2>R session and package versions</h2>",
                sprintf("<pre>%s</pre>",
                        esc(paste(session_info_lines(log_code_packages(entries)),
                                  collapse = "\n"))))
      paste(c(head, body, sess, "</body></html>"), collapse = "\n")
    }

    output$summary <- shiny::renderUI({
      entries <- report_entries()
      if (is.null(entries) || !length(entries)) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("No steps recorded yet (or none in the ticked sections). Run some analysis steps first.",
                               "尚未记录任何步骤（或勾选的章节中没有步骤）。请先运行一些分析步骤。")))
      }
      shiny::tagList(
        stat_tile(i18n("Steps in report", "报告中的步骤数"), length(entries)),
        shiny::tags$ol(class = "omicone-steps",
          lapply(entries, function(e) {
            shiny::tags$li(shiny::tags$b(e$step),
                           shiny::tags$span(class = "omicone-muted",
                                            paste0("  (", e$time, ")")))
          }))
      )
    })

    output$preview <- shiny::renderUI({
      explain_scene("report",
                    paste0("No plot for this step. Tick the sections to include, ",
                             "set a title, then click Download report to save an HTML ",
                             "record of the steps you ran, their parameters, code and ",
                             "the R session."),
                    paste0("这一步没有图表。请勾选要包含的章节、设置标题，",
                             "然后点击“下载报告”，保存一份记录你运行过的步骤、参数、代码与 R 会话的 HTML 报告。"))
    })

    output$download_report <- shiny::downloadHandler(
      filename = function()
        paste0("omicone_report_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".html"),
      content = function(file) {
        title   <- if (nzchar(trimws(input$title %||% ""))) input$title
                   else "OmicOne analysis report"
        entries <- report_entries()

        # Prefer a proper rmarkdown render; fall back to a self-contained HTML
        # string if rmarkdown/pandoc is unavailable or rendering fails.
        rendered <- FALSE
        if (has_pkg("rmarkdown") && rmarkdown::pandoc_available()) {
          rendered <- tryCatch({
            rmd <- tempfile(fileext = ".Rmd")
            writeLines(build_rmd(title, entries), rmd)
            out <- rmarkdown::render(rmd, output_format = "html_document",
                                     output_file = basename(tempfile(fileext = ".html")),
                                     output_dir = dirname(rmd), quiet = TRUE)
            file.copy(out, file, overwrite = TRUE)
            TRUE
          }, error = function(e) FALSE)
        }
        if (!rendered) {
          writeLines(build_html(title, entries), file)
        }
        mark_done(rv, "report")
      }
    )
  })
}
