#' @title simulate codon frequency
#' @description simulates codon frequency
#' @param seq vector of sequences
#' @param f what function should be applied to simulate sequences? e.g. find_codon()
#' @param k what kmer-frequency should be preserved?
#' @param iter how many itereations
#' @param verbose printn status
#' @export

sim_seq <- function(seq, f, k = 2, iter = 1000, verbose = T){

  #use_condaenv('r-reticulate')
  ushuffle <- reticulate::import('ushuffle')
  source_python('python/shuffle_utrs.py')
  seqs <- as.vector(seq)
  simulated_sequences <- shuffle_utrs(seq, k, as.integer(iter))
  res <- lapply(simulated_sequences, function(x) f(x))
  return(res)
}