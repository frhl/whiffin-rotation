#' @title map geneomic coordiantes from RNA to DNA
#' @param df a data.frame with the column \code{bp}.
#' @param coords coordinates


map_mane <- function(coords, bp, adjust_pos = -2){
  
 #if (length(coords) != 0) stop('Expected a data.frame with bp column as input')
 mapping <- create_mane_mapping(coords)
 result <- lapply(bp + adjust_pos, function(x){
   interval <- mapping[mapping$rna_start <= x & mapping$rna_stop >= x,]
   dna_position <- interval$dna_start + (x-interval$rna_start)
   return(dna_position)
 })
 return(unlist(result))
}

