#' Reusable UI building blocks for OmicOne
#'
#' These helpers implement the uniform per-step interaction pattern used by every
#' analysis module. The layout is plot-first: a persistent narrow control rail on
#' the left and a large preview canvas that dominates the view. The explainer is
#' collapsed by default so it never steals space from the plot, and the result
#' summary is a slim strip in the header.
#'
#' @name utils_ui
#' @keywords internal
NULL

#' Bilingual text span (client-side, single language visible)
#'
#' Renders ONE language at a time. Both strings are carried as `data-en` /
#' `data-zh` attributes; only the active language is shown as the element's
#' content. The header switch calls a small JS helper (see app.js) that swaps
#' every `.i18n` element's innerHTML to the chosen language -- so at any moment
#' only one language is on screen (never both), with no server round-trip.
#'
#' @param en,zh English and Chinese strings (may contain inline HTML).
#' @param tag Wrapper tag function (default `shiny::tags$span`).
#' @keywords internal
i18n <- function(en, zh, tag = shiny::tags$span) {
  tag(class = "i18n", `data-en` = en, `data-zh` = zh, shiny::HTML(en))
}

#' Collapsible beginner explainer card ("What is this step?")
#'
#' Collapsed by default (experts ignore it, beginners expand it). Each facet
#' accepts a bilingual pair.
#'
#' @param title Bilingual list(en=, zh=) or single string: short step title.
#' @param what,why,how,read,example Bilingual list(en=, zh=) pairs (how/read/
#'   example optional). `read` explains how to interpret the resulting figure
#'   or table; `example` is highlighted.
#' @param open Logical; start expanded. Default `FALSE` (plot-first).
#' @keywords internal
explainer_card <- function(title, what, why, how = NULL, read = NULL,
                           example = NULL, open = FALSE) {
  bi <- function(x, prefix_en, prefix_zh) {
    if (is.null(x)) return(NULL)
    if (is.list(x)) {
      shiny::tags$p(i18n(paste0("<b>", prefix_en, ":</b> ", x$en),
                         paste0("<b>", prefix_zh, ":</b> ", x$zh)))
    } else {
      shiny::tags$p(shiny::HTML(paste0("<b>", prefix_en, ":</b> ", x)))
    }
  }
  ttl <- if (is.list(title)) title$en else title
  body <- list(
    bi(what, "What", "作用"),
    bi(why,  "Why",  "为什么"),
    bi(how,  "How",  "怎么用"),
    bi(read, "How to read the result", "如何解读结果"),
    if (!is.null(example)) {
      ex <- if (is.list(example)) example else list(en = example, zh = example)
      shiny::div(class = "omicone-example",
                 shiny::tags$span(class = "omicone-example-tag",
                                  i18n("Example", "示例")),
                 i18n(ex$en, ex$zh))
    }
  )
  bslib::accordion(
    class = "omicone-explainer", open = open,
    bslib::accordion_panel(
      title = i18n(paste0("\U0001F4A1 What is this step?"),
                   paste0("\U0001F4A1 这一步是什么？")),
      value = paste0("explainer-", gsub("[^a-z0-9]+", "-", tolower(ttl))),
      icon = NULL, body
    )
  )
}

#' Inline help tooltip on a `?` icon (bilingual)
#' @param en,zh Tooltip text.
#' @keywords internal
help_tip <- function(en, zh = en) {
  bslib::tooltip(shiny::tags$span(class = "omicone-help", "?"),
                 i18n(en, zh), placement = "right")
}

#' A labelled control with an attached help tooltip (bilingual)
#' @param label_en,label_zh Label text. @param tip_en,tip_zh Tooltip text.
#' @keywords internal
label_with_help <- function(label_en, tip_en, label_zh = label_en, tip_zh = tip_en) {
  shiny::tags$label(class = "omicone-label",
                    i18n(label_en, label_zh), help_tip(tip_en, tip_zh))
}

#' Plot-first step container: narrow control rail + large preview
#'
#' @param title Bilingual list(en=, zh=) or string: the step title (header).
#' @param subtitle Bilingual list(en=, zh=) or string: one always-visible line
#'   under the title saying what this step does. Optional.
#' @param explainer An [explainer_card()] (collapsed). Rendered as a full-width
#'   strip between the header and the workspace, so long prose is never
#'   squeezed into the narrow control rail.
#' @param controls Tag list of inputs + the Run button (goes in the left rail).
#' @param summary A UI output slot for the slim result summary (header strip).
#' @param preview A UI output slot for the large preview plot(s).
#' @param rail_width Control-rail width. Default "280px".
#' @keywords internal
step_container <- function(title, explainer, controls, summary, preview,
                           rail_width = "280px", subtitle = NULL) {
  ttl <- if (is.list(title)) i18n(title$en, title$zh) else title
  sub <- NULL
  if (!is.null(subtitle)) {
    sub <- shiny::div(class = "omicone-subtitle",
                      if (is.list(subtitle)) i18n(subtitle$en, subtitle$zh)
                      else subtitle)
  }
  bslib::card(
    class = "omicone-step", full_screen = TRUE,
    bslib::card_header(
      shiny::div(class = "omicone-headwrap",
                 shiny::div(class = "omicone-stephead",
                            shiny::span(class = "omicone-steptitle", ttl),
                            shiny::div(class = "omicone-summarystrip", summary)),
                 sub)
    ),
    if (!is.null(explainer)) explainer,
    bslib::layout_sidebar(
      sidebar = bslib::sidebar(
        width = rail_width, open = "open", position = "left",
        shiny::div(class = "omicone-rail", controls)
      ),
      shiny::div(class = "omicone-plotwrap", preview)
    )
  )
}

#' Primary "Run" button that shows a busy spinner while the step runs
#'
#' Uses [bslib::input_task_button()] so clicking instantly shows a spinner and
#' disables the button until the step finishes -- fixing the "looks like nothing
#' happened, so I clicked again" problem.
#'
#' @param id Namespaced input id.
#' @param en,zh Idle label text.
#' @keywords internal
run_button <- function(id, en = "Run this step", zh = "运行此步骤") {
  bslib::input_task_button(
    id,
    label = i18n(en, zh),
    label_busy = i18n("Running…", "运行中…"),
    icon = shiny::icon("play"),
    class = "btn-primary omicone-run w-100"
  )
}

#' Compact summary pill (used in the slim header strip)
#' @param title,value Character/numeric.
#' @keywords internal
stat_tile <- function(title, value, showcase = NULL) {
  shiny::div(class = "omicone-pill",
             shiny::span(class = "omicone-pill-label", title),
             shiny::span(class = "omicone-pill-value", value))
}
