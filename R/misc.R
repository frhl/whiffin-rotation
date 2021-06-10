#' @title split sequence
#' @description splits string sequence into vectors
#' @param x string
#' @export

split_seq <- function(x) unlist(strsplit(x, split = ''))
