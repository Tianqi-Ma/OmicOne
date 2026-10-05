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
  # Inter is linked from Google Fonts, but the app must look right offline and
  # where Google is unreachable: the browser then falls through the system UI
  # stack (CJK-capable on every OS) instead of a default serif face.
  ui_font <- app_font_stack()

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
          shiny::uiOutput("export_items")
        ),
        shiny::tags$div(
          class = "omicone-lang",
          shiny::tags$button(class = "omicone-lang-btn active", `data-lang` = "en",
                             onclick = "OmicOneSetLang('en')", "EN"),
          shiny::tags$button(class = "omicone-lang-btn", `data-lang` = "zh",
                             onclick = "OmicOneSetLang('zh')", "中")
        ),
        # follow the operating system's light/dark preference; the toggle
        # still switches it for the session
        bslib::input_dark_mode(id = "dark")
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
      shiny::tags$script(src = "omicone/app.js"),
      shiny::tags$script(src = "omicone/explain.js"),
      shiny::tags$script(src = "omicone/explain-sc.js"),
      shiny::tags$script(src = "omicone/explain-wes.js")
    ),
    shiny::uiOutput("stale_banner"),
    shiny::uiOutput("main_body")
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
