#' @title split sequence
#' @description splits string sequence into vectors
#' @param x string
#' @export

split_seq <- function(x) unlist(strsplit(x, split = ''))

#' @title extract starts
#' @param x named list
#' @export
extract_starts <- function(x) {
  return(as.numeric(as.matrix(do.call(rbind, strsplit(names(x),split='_')))[,1]))
}

#' @title extract stops
#' @param x named list
#' @export
extract_stops <- function(x) {
  return(as.numeric(as.matrix(do.call(rbind, strsplit(names(x),split='_')))[,2]))
}

