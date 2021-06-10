#' @title find truncating AUGs
#' @description find AUGs that are not in an open reading frame
#' @param x string, sequence
#' @param start string, regex
#' @param stop string, regex
#' @export

get_truncating_augs <- function(x, start = 'ATG', stop = '(TAG)|(TAA)|(TGA)'){
  
  # truncating AUGs are AUGs that are not in an ORF
  truncations <- find_codon(x)
  orfs_seq <- get_orf(x, start, stop)
  if (length(orfs_seq) > 0){
    orf_names <- unlist(lapply(strsplit(names(orfs_seq), split = '_'), function(x) x[1]))
    orfs <- as.numeric(orf_names)
    truncations <- truncations[!truncations %in% orfs]
  }

  # split sequence and get truncations
  seq_splitted <- split_seq(x)
  seq_out <- lapply(truncations, function(i){
    paste0(seq_splitted[(i-2):nchar(x)], collapse = '')
  })
  names(seq_out) <- paste0(truncations,'_0')
  return(seq_out)
  
}
