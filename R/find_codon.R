

find_codon <- function(x, codon = 'ATG'){

  split_codon <- unlist(strsplit(x, split = codon))
  split_codon_n <- unlist(lapply(split_codon, nchar))+3 
  codons <- cumsum(split_codon_n) 
  return(codons)
  
}


