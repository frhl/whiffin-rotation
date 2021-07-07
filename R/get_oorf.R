#' @title get overlapping open reading frames
#' @description get overlapping open reading frames
#' @param utr utr sequence
#' @param cds cds sequence
#' @param share_stops should all possible sequences be returned? If false, then if two sequences share the same stop,
#' only the one with the strongest kozak will be returned (see \code{?get_kozak_strength()}. If both sequences have the same kozak strength, only the first
#' (5'->3') will be kept.
#' @param inframe boolean should the ORFs be in frame with CDS? Null means that both inframe/out of frame oorfs are returned.
#' @details Some parameteres that could be useful: 
#' 
#' * Out of frame overlapping ORFs (oORFs) - \code{inframe = F} and \code{share_stops = F}
#' * N-terminal extensions (NTE) - \code{inframe = T} and \code{share_stops = F}
#' 
#' 
#' @export

get_oorf <- function(utr, cds, inframe = NULL, share_stops = T){
  
  # check for open reading frames in combined data
  combined <- paste0(utr, cds)
  orfs <- get_orf(combined, share_stops = T) # must always share stops!
  if (length(orfs) > 0){
    
    # check overlap
    utr_end <- nchar(utr)
    cds_end <- nchar(cds)
    mat_orf <- as.data.frame(do.call(rbind, strsplit(names(orfs), split = '_')))
    mat_orf$V1 <- as.numeric(mat_orf$V1)
    mat_orf$V2 <- as.numeric(mat_orf$V2)
    overlap <- mat_orf$V1 < utr_end & mat_orf$V2-3 > utr_end
    
    # in/out of frame to CDS
    if (!is.null(inframe)){
      stopifnot(is.logical(inframe))
      orf_start <- extract_starts(orfs)
      cds_start <- nchar(utr)+3
      inframe_cds <- (orf_start - cds_start) %% 3 == 0 & orf_start < cds_start
      bool <-  inframe_cds & overlap
      return(orfs[unlist(ifelse(inframe, list(bool), list(!bool & overlap)))])
    }

    oorfs <- orfs[overlap]
  
    # deal with shared stops
    if (!share_stops){
      mat_oorf <- as.data.frame(do.call(rbind, strsplit(names(oorfs), split = '_')))
      mat_oorf$V1 <- as.numeric(mat_oorf$V1)
      mat_oorf$V2 <- as.numeric(mat_oorf$V2)
      shared_stops <- mat_oorf$V2[duplicated(mat_oorf$V2)]
      index <- select_shared_stops(x, mat_oorf, shared_stops)
      oorfs <- oorfs[index]
    }
    
    return(oorfs)
  }
  return(NULL)
}




