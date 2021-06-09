




get_oorf <- function(utr, cds){
  
  combined <- paste0(utr, cds)
  orfs <- find_orfs(combined)
  if (length(orfs) > 0){
    utr_end <- nchar(utr)
    cds_end <- nchar(cds)
    mat_orf <- as.data.frame(do.call(rbind, strsplit(names(orfs), split = '_')))
    overlap <- mat_orf$V1 < utr_end & mat_orf$V2 > utr_end
    return(orfs[overlap])
  }
  return(NA)
}


