#' @title select kozak
#' @description a sub function used to select among multiple sequences,
#' the kozak that has the highest kozak score (strongest). Ties will be 
#' resovled by selecting the one furthest away from the CDS start site.
#' @export

select_kozak <- function(x){
  kozak <- get_kozak_strength(x)
  strength <- unlist(kozak)
  kozak_pos <- min(extract_starts(kozak[strength == max(strength)]))
  return(kozak_pos)
}
