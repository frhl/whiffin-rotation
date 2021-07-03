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


#' @title print to port
#' @description log output to port
#'
#' @details 
#' In the following order do:
#' * listen in terminal at: nc -l 4000
#' * set log.socket <- make.socket(port=4000)
#' * start script
#' @export

pprint <- function(text, ...) {
  msg <- sprintf(paste0(as.character(Sys.time()), ": ", text, "\n"), ...)
  cat(msg)
  write.socket(log.socket, msg)
}



#' @title generate codons
#' @export

generate_codons <- function(x=c('A','T','C','G')){
  codons <- c()
  for (i in x) for (j in x) for (l in x) codons <- c(codons, (paste0(i,j,l)))
  return(codons)
}


#' @title get by colname
#' @export
get <- function(w,x) return( (1:ncol(x))[grepl(w,colnames(x))])


#' @title index specific sub-strings in vectors
#' @export
indexsplit <- function(x, i, split = '\\.') unlist(lapply(x, function(y) unlist(strsplit(y, split = split))[i] ))

