



#' @title simulate codon frequency
#' @description simulates codon frequency
#' @param seqs vector of sequences
#' @param codon what codon should be determined
#' @param export

sim_seq <- function(seq, f, k = 2, iter = 1000, verbose = T){

  #use_condaenv('r-reticulate')
  ushuffle <- reticulate::import('ushuffle')
  source_python('python/shuffle_utrs.py')
  seqs <- as.vector(seq)
  simulated_sequences <- shuffle_utrs(seq, k, as.integer(iter))
  res <- lapply(simulated_sequences, function(x) f(x))
  return(res)
}