
devtools::load_all()


d1 <- fread('~/Projects/08_genesets/genesets/data/MANE/210708_MANE.GRCh38.v0.95.combined-table.txt')
d1$bp[d1$bp == '-'] <- NA
d1 <- d1[d1$type == 'five_prime_UTR']
d2 <- fread('derived/tables/210708_MANE.v0.95.UTR_features.txt', sep = '\t')
mrg1 <- merge(d1, d2)
#mrg1 <- mrg[mrg$chr == 1,]
#mrg1 <- mrg[mrg$u5_lenx > 50 & mrg$u5_AUG > 1 & mrg$enstid_version %in% minus,]


variants <- fread('extdata/clinvar/clinvar_20210626_chr1_22.txt')
colnames(variants) <- c('chr','bp','ref','alt')

result <- do.call(rbind, lapply(1:nrow(mrg1), function(i){
  row <- mrg1[i]
  if (!is.na(row$bp)){
    bps <- as.data.frame(do.call(rbind, lapply(unlist(strsplit(row$bp, split = ';')), function(x) unlist(strsplit(x, split = '-')))))
    res <- do.call(rbind, lapply(1:nrow(bps), function(j){
      selected <- as.numeric(unlist(bps[j,]))
      vars <- variants[variants$bp > min(selected) & variants$bp < max(selected)]
      if (nrow(vars) == 0) return(NULL)
      if (row$bp == '-') return(NULL)
      vars$enstid <- row$enstid_version
      vars$interval <- row$bp
      vars$strand <- row$strand
      vars$bp_start <- min(selected)
      vars$bp_end <- max(selected)
      #vars$bp<- row$bp
      vars$cdna_bp <- row$bp_cdna
      
      return(vars)
    }))
    return(res)
  }
  
}))









mrg1[grepl('80654490', mrg1$bp)]



## WORKS FOR DIRECTION = +1

# this matches!!
setup_mapping <- function(x, reference){
  sequence_len <- nchar(x)
  mapping_structure = data.frame(
    index = sequence_len:1,
    bp = reference:(reference-sequence_len+1),
    ref = rev(unlist(strsplit(x, split = '')))
  )
  return(mapping_structure)
  
}


setup_mapping_minus <- function(x, reference){
  sequence_len <- nchar(x)
  mapping_structure = data.frame(
    index = rev(1:sequence_len),
    bp = rev(reference:(reference-sequence_len+1)),
    ref = rev(unlist(strsplit(x, split = '')))
  )
  return(mapping_structure)
  
}


# mrg1[mrg1$enstid_version %in% 'ENST00000370225.4']
res <- result[!is.na(result$bp),]
res <- res[nchar(res$ref) == 1]
res <- res[res$strand== '+',]



bool <- (lapply(1:nrow(res), function(i){
  
  row <- res[i,]
  s <- mrg1[mrg1$enstid_version %in% row$enstid]
  q <- make_mapping(s$seq, s$bp_cdna, s$bp)
  ret <- q[q$bp == row$bp,]$ref %in% row$ref
  
  
  mapped <- q[q$bp == row$bp,]
  colnames(mapped)
  original <- row
  
  
  colnames(original)[1:4] <- paste0('clinvar.',colnames(original)[1:4])
  colnames(original)[5:10] <- paste0('MANE.',colnames(original)[5:10])
  colnames(mapped) <- c('mymapping.cDNA','mymapping.bp','mymapping.ref')
  mapped$mymapping.correctly_mapped <- ret
  result <- cbind(mapped,original)
  
  return(result)
  #row <- res[i,]
  #s <- mrg1[mrg1$enstid_version %in% row$enstid]
  #q <- setup_mapping_minus((s$seq), row$bp_end)
  
  #ret <- q[q$bp == row$bp,]$ref %in% row$ref
  
  
  q
  row
  q[q$bp == row$bp,]
  q[q$bp == row$bp-1,]
  
  #q[q$bp == row$bp,]
  
  #return(ret)
  
}))

status <- as.data.frame(do.call(rbind, bool))
status

#fwrite(status, 'derived/210709_clinvar_mane_mapping_status.txt', sep = '\t')

