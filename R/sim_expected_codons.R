#' @title simulate codon frequency
#' @description simulates codon frequency
#' @param seqs vector of sequences
#' @param codons what codon should be determined
#' @param k what k-mer frequency should be preserved?
#' @param iter how many simulations?
#' @param verbose print status?
#' @param parallel run parallel
#' @param cores how many cores should be used when in parallel?
#' @export

sim_expected_codons <- function(seqs, codons = 'ATG', k = 2, iter = 1000, verbose = T, parallel = F, cores = detectCores(), seed = 1){
  
  stopifnot(length(seqs) > 0)
  #use_condaenv('r-reticulate')
  ushuffle <- reticulate::import('ushuffle')
  source_python('python/shuffle_utrs.py')
  count <- 0
  seqs <- as.vector(seqs)
  
  
  ## parallelize the work
  if (parallel){
    
    #browser()
    
    # setup compute environment
    library(parallel)
    library(doParallel)
    cores <- cores
    registerDoParallel(cores)
    nseqs <- length(seqs)
    
    write(paste('Starting parallel using',cores,'cores...'),stdout())
     
    res <- (foreach (i=1:nseqs) %dopar% {
      
      # load ushuffle
      #use_condaenv('r-reticulate')
      ushuffle <- reticulate::import('ushuffle')
      source_python('python/shuffle_utrs.py')
      
      # prepare sequence
      cur_seq <- seqs[i]
      sim_seq <- shuffle_utrs(cur_seq, k, as.integer(iter), seed)
      
      # We calculated the median and 2.5% and 97.5% quantiles of the number of AUGs in the shuffled 5′ UTRs 
      mat <- do.call(cbind, lapply(codons, function(cur_codon){
        counts <- unlist(lapply(sim_seq, function(x) (count_codon(x, codon = cur_codon))))
        probs <- quantile(counts, probs = c(0.025, 0.5, 0.975))
        expt <- paste0(probs[2],';[',probs[1],'~',probs[3],'];',sum(counts)/iter)
        return(expt)
        
      }))
      colnames(mat) <- codons
      return(mat)
    })
      
  # no parallel
  } else {
    
    # simulate each sequence 
    res <- lapply(seqs, function(seq) {
      sim_seq <- shuffle_utrs(seq, k, as.integer(iter))
      
      # We calculated the median and 2.5% and 97.5% quantiles of the number of AUGs in the shuffled 5′ UTRs 
      mat <- do.call(cbind, lapply(codons, function(cur_codon){
        counts <- unlist(lapply(sim_seq, function(x) (count_codon(x, codon = cur_codon))))
        probs <- quantile(counts, probs = c(0.025, 0.5, 0.975))
        expt <- paste0(probs[2],';[',probs[1],'~',probs[3],'];',sum(counts)/iter)
        return(expt)
        
      }))
      colnames(mat) <- codons
      
      # verbose/ logging
      count <<- count + 1
      if (verbose & count %% 100 == 0) print(paste0(count,' / ',length(seqs)))
      return(mat)
    })
    
  }
  

  return(res)
}
