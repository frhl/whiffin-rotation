#' @title simulate codon frequency
#' @description simulates codon frequency
#' @param seqs vector of sequences
#' @param codon what codon should be determined
#' @param export

sim_expected_codons <- function(seqs, codons = c('ATG'), k = 2, iter = 1000, verbose = T){
  
  count <- 0
  
  # simulate each sequence 
  res <- lapply(seqs, function(seq) {
    sim_seq <- shuffle_utrs(seq, k, as.integer(iter))
    
    # cound occurence of individual codons
    mat <- do.call(cbind, lapply(codons, function(cur_codon){
      expt <- length(unlist(lapply(sim_seq, function(x) find_codon(x, codon = cur_codon)))) / iter
      return(expt)
    }))
    colnames(mat) <- codons
    
    # verbose/ logging
    count <<- count + 1
    if (verbose & count %% 100 == 0) print(paste0(count,' / ',length(seqs)))
    return(mat)
  })
  return(res)
}