#' @title find codons
#' @description cound codons
#' @param x string, sequence
#' @param codon string
#' @export

count_codon <- function(x, codon = 'ATG'){
  return(max(0, length(unlist(strsplit(x, split = codon)))-1))
}


