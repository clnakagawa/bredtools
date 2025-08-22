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
