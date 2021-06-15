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


#' @title log
# listen in terminal at: nc -l 4000
# rember to set log.socket <- make.socket(port=4000)


pprint <- function(text, ...) {
  msg <- sprintf(paste0(as.character(Sys.time()), ": ", text, "\n"), ...)
  cat(msg)
  write.socket(log.socket, msg)
}
