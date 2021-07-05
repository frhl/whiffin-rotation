#' @title create mane mapping data.frame for RNA to DNA
#' @description a function that creates a mapping between genomic coordinate systems.
#' @param df a data.frame with the column \code{bp}
#' @export

create_mane_mapping <- function(df){
  
  stopifnot('bp' %in% colnames(df))
  splitted <- strsplit(df$bp,split = ';')
  stopifnot(length(splitted) == 1) # only handle one sequence for now
  result <- lapply(splitted, function(x) {
    
    # initial prep
    dna <- strsplit(x, split = '-')
    rna <- unlist(lapply(dna, function(y) as.numeric(y[2]) - as.numeric(y[1])))
    
    # dna level
    dna_start = as.numeric(unlist(lapply(dna, function(x) x[1])))
    dna_stop = as.numeric(unlist(lapply(dna, function(x) x[2])))
    dna_len = dna_stop - dna_start
    
    # rna level
    rnas <- c(0, cumsum(rna+1)) # R is one based
    rna_start_stop <- do.call(rbind, lapply(1:(length(rnas)-1), function(i) data.frame(rnas[i]+1, rnas[i+1])))
    colnames(rna_start_stop) <- c('rna_start','rna_stop')
    rna_len <- rna_start_stop$rna_stop - rna_start_stop$rna_start
    mapping <- cbind(dna_start, dna_stop, dna_len, rna_start_stop, rna_len)
    stopifnot(mapping$rna_len == mapping$dna_len)
    return(mapping)
  })[[1]]
  
  return(result)
}

