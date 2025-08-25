# helper functions for common bed file stuff

#' Helper to read in bed file
#' Contributes not much except setting first 3 column names
#'
#' @description
#' `bread` reads in a bed file and sets the first three column names.
#'
#' @details
#' The returned object is in data.table format because it's better that way
#'
#' @param bedfile A path to a bed file
#' @param ... Additional arguments passed to data.table::fread()
#' @returns A data table object
#'
#' @export
bread <- function(bedfile, ...) {
  bdt <- data.table::fread(file = bedfile, ...)
  if (colnames(bdt)[1] == "V1") {
    data.table::setnames(bdt,
                         old = c("V1", "V2", "V3"),
                         new = c("chr", "start", "end"))
  }
  return(bdt)
}

#' @description
#' `bintersect` takes two bed-like data tables and returns intersecting
#' ranges as another data table
#' @details
#' Options mirror bedtools options (e.g. -wa), retains
#' properties (additional columns) of first bed file (unless specified)
#'
#' @param bedA path to the first bed file
#' @param bedB path to the second bed file
#' @returns A data table object
#'
#' @export
bintersect <- function(bedA, bedB) {
  nf <- ncol(bed1)
  overlaps <- foverlaps(bedA, bedB, type = "any", nomatch = 0)
  return(overlaps[, `:=` (start = pmax(start, i.start), end = pmax(end, i.end))][, 1:nf])
}
