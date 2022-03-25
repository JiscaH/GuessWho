# utility functions

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title wrapper for read.table() with better defaults
#'
#' @description has the defaults \code{\link{read.table}} should have had.
#'
#' @param file  the name of the file which the data are to be read from.
#' @param sep  the field separator character.
#' @param header a logical value indicating whether the file contains the names
#'   of the variables as its first line.
#' @param na.strings a character vector of strings which are to be interpreted
#'   as NA values.
#' @param ... further arguments to \code{\link{read.table}}
#'
#' @details  stringsAsFactors=FALSE.
#'
#' @return A data frame (data.frame) containing a representation of the data in the file.
#'
#' @export

ReadTable <- function(file, sep="\t", header=TRUE, na.strings=c("", NA), ...)
  utils::read.table(file,
             header=header,
             stringsAsFactors=FALSE,
             sep=sep,
             na.strings=na.strings,
             ...)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Max <- function(x, m=NA)  ifelse(all(is.na(x)), m, max(x, na.rm=TRUE))

Min <- function(x, m=NA)  ifelse(all(is.na(x)), m, min(x, na.rm=TRUE))

Max2 <- function(x, m=NA) Max(x[-which.max(x)], m=m)   # 2nd largest value

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
