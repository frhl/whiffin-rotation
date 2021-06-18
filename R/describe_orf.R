#' @title describe orf
#' @description describe all orfs in a sequence, by converting them into a data.frame
#' @param x string, the sequence.
#' @param how vector/string. Takes "kozak" and "".
#' @param start start codon
#' @param stop stop codon
#' @details depends on \code{get_orf} and \code{get_kozak_strength}.
#' @export

describe_orf <- function(x, how = c('kozak'), start = 'ATG', stop = '(TAG)|(TAA)|(TGA)'){
 
  # get standard sequence stats
  orfs <- get_orf(x, start, stop)
  df <- data.frame(stack(orfs), start = extract_starts(orfs), stop = extract_stops(orfs))
  colnames(df)[1:2] <- c('seq','bp')
  
  # kozak descriptor
  if (tolower(how) %in% 'kozak'){
    kozak <- get_kozak_strength(x, only_orf = T)
    d_kozak <- stack(kozak)
    colnames(d_kozak) <- c('kozak_strength', 'start')
    df <- merge(df, d_kozak, by = 'start')
    df <- df[,c('seq','start','stop','bp','kozak_strength')]
  }
  
  df$bp <- NULL
  return(df)  
}



