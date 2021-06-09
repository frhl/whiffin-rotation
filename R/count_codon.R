#' @title find codons


count_codon <- function(x, codon = 'ATG'){
  return(max(0, length(unlist(strsplit(x, split = codon)))-1))
}


