#' Top-level UI: omics landing page + per-omics plot-first pipeline
#'
#' On launch a splash animation plays, then a landing page lets the user pick an
#' omics type. Choosing one routes into that pipeline (a grouped left stepper +
#' the step's plot-first workspace). The single-cell pipeline is fully
#' implemented; other omics show their planned roadmap via placeholders.
#'
#' @return A [bslib::page_sidebar()] UI wrapped with the splash overlay.
#' @keywords internal
app_ui <- function() {
  # Inter from Google Fonts, but this app is meant to work offline: if the font
  # cannot be resolved, fall back to the system UI stack instead of failing.
  ui_font <- tryCatch(
    bslib::font_google("Inter", local = FALSE),
    error = function(e) bslib::font_collection(
      "Inter", "system-ui", "-apple-system", "Segoe UI", "Roboto",
      "Helvetica Neue", "Arial", "sans-serif"))

  theme <- bslib::bs_theme(
    version = 5, preset = "shiny",
    primary = "#2f81c7",
    base_font = ui_font,
    heading_font = ui_font
  )

  page <- bslib::page_sidebar(
    theme = theme,
    title = shiny::div(
      class = "omicone-topbar",
      shiny::span(class = "omicone-brand",
                  shiny::strong("OmicOne"),
                  shiny::span(class = "omicone-sub", "multi-omics, locally")),
      shiny::uiOutput("progress_chip", inline = TRUE),
      shiny::div(
        class = "omicone-topright",
        shiny::actionLink("switch_omics", i18n("← Omics", "← 切换组学"),
                          class = "omicone-switch"),
        bslib::popover(
          shiny::actionLink("export_menu", i18n("⤓ Export", "⤓ 导出"),
                            class = "omicone-switch"),
          title = i18n("Export current data", "导出当前数据"),
          shiny::downloadButton("dl_rds",   i18n("Object (.rds)", "对象 (.rds)"), class = "btn-sm w-100 mb-1"),
          shiny::downloadButton("dl_meta",  i18n("Cell metadata (.csv)", "细胞元数据 (.csv)"), class = "btn-sm w-100 mb-1"),
          shiny::downloadButton("dl_matrix",i18n("Counts matrix (.rds)", "表达矩阵 (.rds)"), class = "btn-sm w-100 mb-1"),
          shiny::downloadButton("dl_embed", i18n("Embeddings (.csv)", "降维坐标 (.csv)"), class = "btn-sm w-100")
        ),
        shiny::div(class = "omicone-pydot", id = "py-status",
                   i18n("Python: not set up", "Python：未配置")),
        shiny::tags$div(
          class = "omicone-lang",
          shiny::tags$button(class = "omicone-lang-btn active", `data-lang` = "en",
                             onclick = "OmicOneSetLang('en')", "EN"),
          shiny::tags$button(class = "omicone-lang-btn", `data-lang` = "zh",
                             onclick = "OmicOneSetLang('zh')", "中")
        ),
        bslib::input_dark_mode(id = "dark", mode = "dark")
      )
    ),
    sidebar = bslib::sidebar(
      title = i18n("Workflow", "分析流程"),
      width = 232, open = "open", id = "stepbar",
      shiny::uiOutput("step_nav"),
      shiny::hr(),
      shiny::div(class = "omicone-mini", shiny::uiOutput("global_status"))
    ),
    shiny::tags$head(
      shiny::tags$link(rel = "stylesheet", type = "text/css", href = "omicone/custom.css"),
      shiny::tags$script(src = "omicone/app.js")
    ),
    shiny::uiOutput("main_body"),
    shiny::uiOutput("next_hint")
  )

  shiny::tagList(
    shiny::div(
      id = "omicone-splash",
      shiny::tags$canvas(id = "omicone-splash-canvas"),
      shiny::div(id = "omicone-splash-logo",
                 shiny::div(class = "t", "OmicOne"),
                 shiny::div(class = "s", "multi-omics analysis, on your own machine"))
    ),
    page
  )
}
