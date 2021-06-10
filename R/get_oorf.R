#' @title get overlapping open reading frames
#' @description get overlapping open reading frames
#' @param utr utr sequence
#' @param cds cds sequence
#' @param inframe boolean should the ORFs be in frame with CDS?
#' @export

get_oorf <- function(utr, cds, inframe = NULL){
  
  # check for open reading frames in combined data
  combined <- paste0(utr, cds)
  orfs <- get_orf(combined)
  
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
    return(orfs[overlap])
  }
  return(NULL)
}




