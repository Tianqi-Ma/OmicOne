#' WES module: Report
#'
#' The report step of the WES pipeline: the shared report module limited to the
#' WES log (import, filters, landscape, signatures, prognosis), plus the
#' runnable R script of the steps that ran.
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_wes_report
NULL

#' @rdname mod_wes_report
#' @keywords internal
mod_wes_report_ui <- function(id) mod_report_ui(id, omics = "wes")

#' @rdname mod_wes_report
#' @keywords internal
mod_wes_report_server <- function(id, rv, log_rv) mod_report_server(id, rv, log_rv, omics = "wes")
