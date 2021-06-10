#' @title kmer
#' @description convers sequence into kmers
#' @param x string, seqeuence
#' @param k kmer length
#' @export

kmer <- function(x, k = 3){
  split_x <- unlist(strsplit(x, split = ''))
  kmers <- unlist(lapply(seq(k,nchar(x)+k, by = k), function(i){
    return(paste0(na.omit(split_x[(i-k+1):i]), collapse =''))
  }))
  return(kmers[kmers != ''])
}
