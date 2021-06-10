#' @title omit nulls from list
#' @param x list
#' @export

null_omit <- function(x) return(x[!sapply(x,is.null)])