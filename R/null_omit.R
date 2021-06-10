#' @title omit nulls from list
#' @param list
#' @export

null_omit <- function(x) return(x[!sapply(x,is.null)])