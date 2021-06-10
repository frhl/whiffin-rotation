#' @title find open reading frames
#' @description takes in a string and finds open reading frames (ORFs)
#' @param x a string
#' @param start what start codon should be used
#' @param stop what stop codons should be used
#' @note This function will use the first stoped codon matching a start codon.
#' The rest of the stop codons (in-frame) will not be considered for the corresponding start codon.
#' @export

find_orfs <- function(x, start = 'ATG', stop = '(TAG)|(TAA)|(TGA)'){
  
  # Check for in-frame codons
  start_pos <- find_codon(x, start)
  stop_pos <- find_codon(x, stop)
  outlist <- list()
  if (!is.null(start_pos) & !is.null(stop_pos)){
    
    mat <- expand.grid(start_pos, stop_pos)
    mat <- mat[mat$Var2 > mat$Var1,]
    mat <- mat[(mat$Var2 - mat$Var1) %% 3 == 0, ] # keep in frame
    mat <- mat[!duplicated(mat$Var1),] # remove many matching stop codons
    
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



