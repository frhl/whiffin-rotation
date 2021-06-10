#' @title get altered CDS
#' @description a function to get CDS that are either 1) elongated or 2) truncated or 3) both.
#' @param utr string. UTR sequence
#' @param cds string. CDS sequence
#' @export

get_altered_cds <- function(utr,cds){
  
  browser()
  
  truncation <- get_truncating_augs(utr)
  if (!is.null(truncation)) trunc_start <- extract_starts(truncation)
  
  elongation <- get_oorf(utr, cds)
  if (!is.null(elongation)) elong_start <- extract_starts(elongation)
  
  
  
}