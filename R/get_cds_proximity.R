#' @title get proximity to CDS
#' @param x sequence
#' @param fun what function should be appliedd
#' @description finds the distance to the CDS (assuming that it will be the end of x)
#' @export

get_cds_proximity <- function(x, fun = function(x) get_orf(x, share_stops = F)){
  fun_res <- fun(x)
  if (length(fun_res) > 0) return(nchar(x)-extract_stops(fun_res)) else return(NA)
}

#' @title get proximity to leader/cap
#' @param x sequence
#' @param fun what function should be appliedd
#' @export

get_leader_proximity <- function(x, fun = function(x) get_orf(x, share_stops = F)){
  fun_res <- fun(x)
  if (length(fun_res) > 0) return(extract_starts(fun_res)-3) else return(NA)
  #starts <- extract_starts(fun(x))
  #if (length(starts) > 0)  return(starts-3) else return(NA)
}