#' @title Melt iterable triangle
#' @description creates a long matrix that contains the paired coordinates
#' for non-redudant iteration.
#' @param cols integer. 
#' @param rows integer. Default is NULL, i.e. square matrix.
#' @param how string. "upper" or "lower".
#' @param diagonal boolean. Should diagonal be included?
#' @export

melt_tri <- function(cols, rows = NULL, how = 'upper', diagonal = FALSE){
  
  # check input
  require(reshape2)
  stopifnot(how %in% c('upper','lower'))
  stopifnot(is.numeric(cols) & length(cols) == 1)
  if (is.null(rows)) rows <- cols
  
  # get long matrix for allowed iterations
  m <- matrix(cols*rows, ncol = cols, nrow = rows)
  triangle <- (ifelse(how == 'upper', 
                      list(upper.tri(m, diag = diagonal)),
                      list(lower.tri(m, diag = diagonal))))
  melted <- melt(triangle[[1]])
  colnames(melted) <- c('i','j','iter')
  return(melted)
  
}