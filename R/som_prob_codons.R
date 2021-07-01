#' @title simulate codon frequency
#' @description simulates codon frequency
#' @param seqs vector of sequences
#' @param codon what codon should be determined
#' @param export

sim_prob_codons <- function(seqs, codons = 'ATG', k = 2, iter = 1000, verbose = T){
  
  #use_condaenv('r-reticulate')
  ushuffle <- reticulate::import('ushuffle')
  source_python('python/shuffle_utrs.py')
  count <- 0
  seqs <- as.vector(seqs)
  
  # simulate each sequence 
  res <- lapply(seqs, function(seq) {
    sim_seq <- shuffle_utrs(seq, k, as.integer(iter))
    
    # We calculated the median and 2.5% and 97.5% quantiles of the number of AUGs in the shuffled 5′ UTRs 
    mat <- do.call(cbind, lapply(codons, function(cur_codon){
      counts <- unlist(lapply(sim_seq, function(x) (count_codon(x, codon = cur_codon) > 0)))
      prob <- sum(counts) / length(counts)
      return(prob)

    }))
    colnames(mat) <- codons
    
    # verbose/ logging
    count <<- count + 1
    if (verbose & count %% 100 == 0) print(paste0(count,' / ',length(seqs)))
    return(mat)
  })
  return(res)
}