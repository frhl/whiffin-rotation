#' @title get overlapping open reading frames kozak strength
#' @description find kozak strength of overlapping reading frames.
#' @return a value of kozak strength (3 strong, 2 moderate and 1 weak)
#' @param utr string. UTR sequence.
#' @param cds string. CDS sequence
#' @param f function to be used that takes two arguments, 
#' Ideally, it should be \code{get_altered_cds} or \code{get_oorf}.
#' @export


get_oorf_kozak <- function(utr, cds, f = get_altered_cds){
  
  kozak_strength <- list()
  combined <- paste0(utr, cds)
  kozak <- get_kozak_strength(combined)
  oorf <- f(utr, cds)
  if (length(oorf) > 0){
    kozak_start <- extract_starts(kozak)
    oorf_start <- extract_starts(oorf)
    kozak_strength <- kozak[kozak_start %in% oorf_start]
  }
  return(kozak_strength)
}

  
  



