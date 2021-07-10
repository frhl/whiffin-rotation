#' @title get pre codon
#' @param codon string
#' @param nts nucleotides
#' @export

get_pre_codons <- function(codon, nts = c('A','T','C','G')){
  
  current <- split_seq(codon)
  result <- do.call(rbind, lapply(1:length(current), function(i){
    new_codon <- current
    do.call(rbind, lapply(nts, function(nt){
      new_codon[i] <- nt
      data.frame(
        codon = paste0(new_codon, collapse = ''),
        position = i,
        before = new_codon[i],
        after = current[i]
      )
    }))
  }))
  
  result <- result[result$before != result$after,]
  return(result)
  
  
  
}