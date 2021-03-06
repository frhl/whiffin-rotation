#' @title get kozak strength
#' @description get kozak strength (Assumes AUG as start codon)
#' @param x sequence
#' @param only_orf Boolean. only return open reading frames.
#' @note currently this only check for kozaks in open reading frames.
#' @details The kozak sequence hierarchy has been used (Strong to weak):
#'  
#' \itemize{
#'  \item \code{strong_1 = '(A|G)..ATGG'}
#'  \item \code{moderate_1 =   '(A|G)..ATG.'} 
#'  \item \code{moderate_2 =   '...ATGG'} 
#'  \item \code{weak_1 =   '...ATG.'} 
#'  \item \code{weak_2 = 'ATG'}
#' }
#' 
#' @references Whiffin et al, 2020
#' @export

get_kozak_strength <- function(x, only_orf = F){
  
  # strengths
  kozaks <- list(
    strong_1 = '(A|G)..ATGG',
    moderate_1 =   '(A|G)..ATG.',
    moderate_2 =   '...ATGG',
    weak_1 =   '...ATG.',
    weak_2 = 'ATG'
  )
  
  # specify whether ORF is needed or not
  starts <- find_codon(x)
  if (only_orf){
    orfs <- get_orf(x)
    if (length(orfs) > 0){
      starts <- extract_starts(orfs)
    } else {
      return(NULL)
    }
  }
  
  # process kozaks
  if (!is.null(starts)){
    sequence <- split_seq(x)
    strengths <- lapply(starts, function(i){
      i1 <- max(1, i-5)
      i2 <- i+1
      context <- paste0(sequence[i1:i2], collapse = '')
      kozak_match <- unlist(lapply(kozaks, function(k){
        grepl(k, context)
      }))
      kozak_number <- max(c(3,2,2,1,1)[kozak_match])
      return(kozak_number)
    })
    
    names(strengths) <- starts
    return(strengths)
    
  }
 
}
