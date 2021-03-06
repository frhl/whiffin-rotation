#' @title find open reading frames
#' @description takes in a string and finds open reading frames (ORFs)
#' @param x a string
#' @param start what start codon should be used
#' @param stop what stop codons should be used
#' @param share_stops boolean. If true (default), then ALL open reading frames are returned. However,
#' if false, when multiple start codons have the same stop codon, only the one with the strongest kozak 
#' will be returned.
#' @note This function will use the first stoped codon matching a start codon.
#' The rest of the stop codons (in-frame) will not be considered for the corresponding start codon.
#' @export

get_orf <- function(x, start = 'ATG', stop = '(TAG)|(TAA)|(TGA)', share_stops = F){
  
  # Check for in-frame codons
  start_pos <- find_codon(x, start)
  stop_pos <- find_codon(x, stop)
  outlist <- list()
  if (!is.null(start_pos) & !is.null(stop_pos)){
    
    # get data.frame of starts->stops
    mat <- expand.grid(start_pos, stop_pos)
    mat <- mat[mat$Var2 > mat$Var1,]
    mat <- mat[(mat$Var2 - mat$Var1) %% 3 == 0, ] # keep in frame
    mat <- mat[!duplicated(mat$Var1),] # remove many matching stop codons
    
    # * remove shared stops by selecting strongest kozak
    # * a tie between kozak results in the first one being returned
    if (!share_stops){
      
      # find shared stops and select them
      stop_shared <- mat$Var2[duplicated(mat$Var2)]
      if (length(stop_shared) > 0){
        indexes <- select_shared_stops(x, mat, stop_shared)
        mat <- mat[indexes,]
      }

      
    }
    
    # return sequences if open reading frame exists
    if (nrow(mat) > 0){
      seq_x <- split_seq(x)
      positions <- as.character(apply(mat, 1, function(x) paste(x, collapse = '_')))
      outlist <- lapply(1:nrow(mat), function(i){
        paste0(seq_x[(mat$Var1[i]-2):(mat$Var2[i])], collapse = '')
      })
      names(outlist) <- positions
    }
  }
  return(outlist)
}



