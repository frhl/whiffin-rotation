#' @title find open reading frames
#' @description takes in a string and finds open reading frames (ORFs)
#' @param x a string
#' @param start what start codon should be used
#' @param stop what stop codons should be used
#' @note This function will use the first stoped codon matching a start codon.
#' The rest of the stop codons (in-frame) will not be considered for the corresponding start codon.
#' @export


find_orfs <- function(x, start = 'ATG', stop = c('TAG', 'TAA', 'TGA')){
  
  seq <- unlist(strsplit(x, split = ''))
  
  # get start codons
  split_start <- unlist(strsplit(x, split = start))
  split_start_n <- unlist(lapply(split_start, nchar))+3 
  starts <- cumsum(split_start_n) # start codon = i-2 : i .. e.g. seq[(1610-2):1610]
  
  if (length(starts) > 0) {
    
    # get stop codons
    split_stop <- unlist(strsplit(x, split = stop))
    split_stop_n <- unlist(lapply(split_stop, nchar))+3 
    stops <- cumsum(split_stop_n) # start codon = i-2 : i .. e.g. seq[(1927-2):1927]
    
    # find matching stop-start codons (ORFs)
    orfs <- lapply(starts, function(cur_start){
      
      # look through stops 
      valid_stops <- stops[stops > cur_start & stops <= length(seq)]
      inframe_stops_bool <- (valid_stops-cur_start) %% 3 == 0
      
      if (any(inframe_stops_bool)){
        #print(cur_start)
        inframe_stops <- valid_stops[inframe_stops_bool] 
        inframe_stop <- inframe_stops[1]
        
        # assertions
        kmers <- seq[(cur_start-2):(inframe_stop)]
        kmers <- unlist(lapply(seq(3, length(kmers), by = 3), function(i){
          paste0(kmers[(i-2):i], collapse = '')
        }))
        stopifnot(kmers[1] == start)
        stopifnot(kmers[length(kmers)] %in% stop)
        orf <- paste0(kmers, collapse = '')
        return(orf)
      }
      
      return(NULL)
    })
    
    # get end of sequence
    ends <- as.numeric(unlist(lapply(orfs, function(x) if (!is.null(x)) return(nchar(x)) else return(NA))))
    names(orfs) <- paste0(starts,'_',starts+ends)
    orfs <- null_omit(orfs)
    return(orfs)
  }
  
  return(NA)

}
  
  
  