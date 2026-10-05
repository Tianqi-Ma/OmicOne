#' Input helpers shared by every module
#'
#' Small, pure helpers for the recurring control patterns: parsing a typed
#' gene list, keeping a user's choice when a selector is rebuilt, and picking
#' which metadata columns can sensibly act as a grouping / batch variable.
#'
#' @name utils_inputs
#' @keywords internal
NULL

#' Parse a typed gene list
#'
#' Accepts commas, semicolons, whitespace and new lines as separators, drops
#' empties and duplicates, and keeps the order the user typed.
#' @param x Character scalar (a text input value).
#' @return Character vector of gene symbols (possibly empty).
#' @keywords internal
parse_genes <- function(x) {
  if (is.null(x) || !length(x) || all(is.na(x))) return(character(0))
  g <- trimws(unlist(strsplit(paste(x, collapse = ","), "[,;[:space:]]+")))
  unique(g[nzchar(g)])
}

#' The selection to show when a selector's choices are rebuilt
#'
#' A `renderUI()` that depends on `rv$obj` re-runs after every step; without
#' this the user's choice snaps back to the default each time. Call it with
#' `shiny::isolate(input$x)` as `current`.
#' @param current The current input value (may be NULL).
#' @param choices The new choices (values).
#' @param default Value to fall back to.
#' @keywords internal
keep_selected <- function(current, choices, default = choices[1]) {
  if (!is.null(current) && length(current) && all(current %in% choices)) current else default
}

#' Metadata columns usable as a categorical grouping or batch variable
#'
#' Character / factor / logical columns, plus integer-like numeric columns
#' with few distinct values, that have at least `min_levels` levels and at
#' most `max_levels`. Continuous QC metrics (nCount_RNA, percent.mt, scores)
#' are excluded, so they can never be offered as a batch.
#' @param md A data.frame of cell (or sample) metadata.
#' @param min_levels,max_levels Allowed number of distinct non-missing values.
#' @param max_numeric_levels Numeric columns qualify only up to this many values.
#' @param allow_na Keep columns that contain missing values.
#' @param exclude Regex of column names never offered (QC metrics, scores).
#' @return Character vector of column names.
#' @keywords internal
categorical_cols <- function(md, min_levels = 2, max_levels = 200,
                             max_numeric_levels = 50, allow_na = FALSE,
                             exclude = "^(nCount|nFeature)_|^percent[._]|[._]?[Ss]core$|^doublet_score$") {
  if (is.null(md) || !ncol(md)) return(character(0))
  # on a small object an integer QC metric (nFeature_RNA) can have few enough
  # distinct values to look categorical; never offer metrics as a grouping
  if (length(exclude) && nzchar(exclude)) md <- md[, !grepl(exclude, names(md)), drop = FALSE]
  if (!ncol(md)) return(character(0))
  ok <- vapply(md, function(v) {
    if (!allow_na && anyNA(v)) return(FALSE)
    n <- length(unique(v[!is.na(v)]))
    is_cat <- is.character(v) || is.factor(v) || is.logical(v) ||
      (is.numeric(v) && all(v[!is.na(v)] == round(v[!is.na(v)])) && n <= max_numeric_levels)
    is_cat && n >= min_levels && n <= max_levels
  }, logical(1))
  names(md)[ok]
}
