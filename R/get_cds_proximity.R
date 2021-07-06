#' @title get proximity to CDS
#' @description finds the distance to the CDS (assuming that it will be the end of x)
#' @export

get_cds_proximity <- function(x, fun = function(x) get_orf(x, share_stops = F)){
  stops <- extract_stops(fun(x))
  return(nchar(x)-stops)
}

#' @title get proximity to leader/cap
#' @export

get_leader_proximity <- function(x, fun = function(x) get_orf(x, share_stops = F)){
  starts <- extract_starts(fun(x))
  return(starts-3)
}