#' @title make cdna to grch38 mapping
#' @description useds information from mane to make a mapping
#' @param seq srting. RNA seq 
#' @param cdna string cdna postions, e.g. "1-5;6;10"
#' @param bp string grch38 positions, e.g. "105-109;204-208"
#' @param export

make_mapping <- function(seq, cdna, bp){
  
  stopifnot(length(cdna) == 1)
  stopifnot(length(bp) == 1)
  stopifnot(length(seq) == 1)
  
  # setup matrix of cdna and grch38 coordinates
  mat_cdna <- lapply(unlist(strsplit(cdna, split = ';')), function(x) unlist(strsplit(x, split = '-')))
  mat_cdna <- as.data.frame(do.call(rbind, mat_cdna))
  mat_bp <- lapply(unlist(strsplit(bp, split = ';')), function(x) unlist(strsplit(x, split = '-')))
  mat_bp <- as.data.frame(do.call(rbind, mat_bp))
  stopifnot(nrow(mat_cdna) == nrow(mat_bp))
  
  # split seq and do mappings exon-wise
  splitted_seq <- split_seq(seq)
  mapping_lst <- lapply(1:nrow(mat_bp), function(i){
    
    # get interval
    grch_interval <- mat_bp[i,]
    grch_chunk <- grch_interval$V1:grch_interval$V2
    cdna_interval <- mat_cdna[i,]
    cdna_chunk <- cdna_interval$V1:cdna_interval$V2
    exon_start <- cdna_chunk * 0
    exon_start[1] <- 1
    
    return(data.frame(cdna = cdna_chunk,
               bp = grch_chunk,
               ref = splitted_seq[cdna_chunk],
               exon_start = exon_start))
    
    
  })
  
  mat_out <- as.data.frame(do.call(rbind, mapping_lst))
  return(mat_out)
  
}
