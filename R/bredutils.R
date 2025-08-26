# helper functions for common bed file stuff

#' @description
#' `bsum` gives the sum of all interval lengths
#'
#'  @param beddt A bed-like data table
#'  @returns A numeric value
#'
#'  @export
bsum <- function(beddt) {
  return(sum(rowSums(beddt[,2:3])))
}

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
#' `bwrite` is just fwrite with some options preset
#'
#' @param beddt Data stored in a data.table object
#' @param bedfile A path to a bed file
#' @param ... Additional arguments passed to data.table::fwrite()
#' @returns nothing
#'
#' @export
bwrite <- function(beddt, bedfile, ...) {
  data.table::fwrite(beddt, bedfile, quote = F, col.names = F, sep = '\t')
}

#' @description
#' `bmerge` takes one bed-like data table and condenses overlapping ranges
#' @details
#' keeps properties of first row in overlapping group by default
#'
#' @param beddt bed-like data.table object
#' @returns A data table object
#'
#' @export
bmerge <- function(beddt) {
  # make copy of dt so changes don't carry over
  tmp <- copy(beddt)

  # sort data.table
  setorder(tmp, chr, start, end)

  # label groups of overlapping regions
  tmp[, group := cumsum(cummax(data.table::shift(end, fill = start[1])) < start),
      by = chr]

  # handle additional columns
  nf <- ncol(tmp)
  otherCols <- names(tmp)[4:nf]

  # collapse overlapping rows
  tmp <- tmp[, c(.(start = min(start), end = max(end)),
                 lapply(.SD, first)),
             by = .(chr, group),
             .SDcols = otherCols][, !'group']
  return(tmp)
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
  nf <- ncol(bedA)
  data.table::setkey(bedA, chr, start, end)
  data.table::setkey(bedB, chr, start, end)
  overlaps <- data.table::foverlaps(bedA, bedB, type = "any", nomatch = 0)
  return(unique(overlaps[, `:=` (start = pmax(start, i.start),
                                 end = pmax(end, i.end))][, 1:nf],
                by = c("chr", "start", "end")))
}

#' @description
#' `bunion` takes two bed-like data tables and gives the union of all regions
#' @details
#' Preserved properties of the first entry in a union "group"
#'
#' @param bedA path to the first bed file
#' @param bedB path to the second bed file
#' @returns A data table object
#'
#' @export
bunion <- function(bedA, bedB) {
  bedAll <- rbindlist(list(bedA, bedB), use.names = T)
  return(bmerge(bedAll))
}

#' @description
#' `bjaccard` takes two bed-like data tables and gives the jaccard "distance"
#' Note that we calculate this as 1 - (length of intersection / length of union)
#'
#' @param bedA path to the first bed file
#' @param bedB path to the second bed file
#' @returns A numeric value between 0 and 1
#'
#' @export
bjaccard <- function(bedA, bedB) {
  bedU <- bsum(bunion(bedA, bedB))
  bedI <- bsum(bintersect(bedA, bedB))
  return(1 - bedI / bedU)
}
