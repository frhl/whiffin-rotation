#' @title get altered CDS
#' @description a function to get CDS that are either 1) elongated or 2) truncated or 3) both.
#' @param utr string. UTR sequence
#' @param cds string. CDS sequence
#' @export

get_altered_cds <- function(utr,cds){
  
  # without matching stop
  elongation <- get_oorf(utr, cds)
  truncation <- lapply(get_truncating_augs(utr), function(x) paste0(x, cds))
  combined <- c(elongation, truncation)
  combined <- combined[!duplicated(unlist(combined))]
  return(combined)
  
}