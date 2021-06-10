#' @title find codon
#' @description use regex to return position of codon
#' @param x string (dna sequence)
#' @param codon string, regex.
#' @export

find_codon <- function(x, codon = 'ATG'){
  codon <- as.numeric(gregexpr(codon, x)[[1]])+2
  if (all(codon == 1)){
    return(NULL)
  } else {
    return(codon)
  }
}


