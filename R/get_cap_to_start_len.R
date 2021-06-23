#' @title get length to start position
#' @description get length to start position
#' @param x string. a sequence
#' @param f a function that returns a named list with names indicating starts.
#' E.g. \code{f = get_orf} will return all open reading frames, and this function
#' will return the position of the first start codon in an ORF.
#' @export

get_cap_to_start_len <- function(x, f = get_orf){
  start_pos <- f(x)
  if (length(start_pos) > 0) extract_starts(start_pos)[[1]]-3 else return(NA)
}
