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
#' @param text what text should be printed
#' @param ... extra text
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
#' @param x letters available
#' @export

generate_codons <- function(x=c('A','T','C','G')){
  codons <- c()
  for (i in x) for (j in x) for (l in x) codons <- c(codons, (paste0(i,j,l)))
  return(codons)
}


#' @title get by colname
#' @param w what data.frame?
#' @param x what colnames
#' @export
get <- function(w,x) return( (1:ncol(x))[grepl(w,colnames(x))])


#' @title index specific sub-strings in vectors
#' @param x vector of strings
#' @param i index
#' @param split split on what?
#' @export
indexsplit <- function(x, i, split = '\\.') unlist(lapply(x, function(y) unlist(strsplit(y, split = split))[i] ))

#' @title remove version
#' @param x vector of strings
#' @export
wo_version <- function(x) unlist(lapply(strsplit(x, split = '\\.'), function(x) x[1]))

#' @title standard error of mean
sem <- function(x) sd(x)/sqrt(length(x))

#' @title not in
#' @description returns true for x not in y
#' @param x value x
#' @param y value y
#' @family misc
#' @export
'%nin%' <- function(x, y) !(x %in% y)

get_time <- function(){
  Sys.time()
}






