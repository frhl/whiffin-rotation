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
      
      # find shared stops
      stop_shared <- mat$Var2[duplicated(mat$Var2)]
      if (length(stop_shared) > 0){
      
        kozak <- stack(get_kozak_strength(x))
        colnames(kozak) <- c('str','Var1')
        kozak_mat <- merge(mat, kozak)
        kozak_mat$index <- 1:nrow(kozak_mat)
        
        keep <- unlist(lapply(unique(stop_shared), function(i){
          tmp_mat <- kozak_mat[kozak_mat$Var2 == i,]
          tmp_mat <- tmp_mat[tmp_mat$str == max(tmp_mat$str),]
          tmp_mat <- tmp_mat[tmp_mat$Var1 == min(tmp_mat$Var1),]
          return(tmp_mat$index)
        }))
        
        mat <- mat[keep,]
        
      }
      
      #kozak_pos <- select_kozak(x)
      
      
      #mat <- mat[mat$Var1 %in% kozak_pos,]
      
      #mat <- mat[!duplicated(mat$Var2),]
      #kozak <- get_kozak_strength(x)
      #d <- data.frame(Var1 = extract_starts(kozak), Var2 = mat$Var2, strength = unlist(kozak))
      #keep_mat <- do.call(rbind, lapply(unique(d$Var2), function(s) {
      #  tmp <- d[d$Var2 == s,]
      #  return(tmp[ min(which(tmp$strength == max(tmp$strength))), ])
      #}))
      #keep_mat$strength <- NULL
      #mat <- keep_mat
      
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



