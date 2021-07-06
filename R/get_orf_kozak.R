#' @title get open reading frames kozak strength
#' @description find kozak strength of overlapping reading frames.
#' @return a value of kozak strength (3 strong, 2 moderate and 1 weak)
#' @param x string. UTR sequence.
#' @param f function to be used that takes two arguments, 
#' Ideally, it should be \code{get_orf}.
#' @export


get_orf_kozak <- function(x, f = function(y) get_orf(y, share_stops = T)){

  kozak_strength <- list()
  kozak <- get_kozak_strength(x)
  orf <- f(x)
  if (length(orf) > 0){
    kozak_start <- extract_starts(kozak)
    orf_start <- extract_starts(orf)
    kozak_strength <- kozak[kozak_start %in% orf_start]
  }
  return(kozak_strength)
}




