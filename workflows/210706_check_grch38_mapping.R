
devtools::load_all()


d1 <- fread('~/Projects/08_genesets/genesets/data/MANE/210709_MANE.GRCh38.v0.95.combined-table.txt')
d1$bp[d1$bp == '-'] <- NA
d1 <- d1[d1$type == 'five_prime_UTR']
d2 <- fread('derived/tables/210708_MANE.v0.95.UTR_features.txt', sep = '\t')
mrg <- merge(d1, d2)
mrg1 <- mrg[mrg$chr == 1]


variants <- fread('extdata/clinvar/clinvar_20210626_chr1.txt')
colnames(variants) <- c('chr','bp','ref','alt')

result_list <- lapply(1:nrow(mrg), function(i){
  row <- mrg[i]
  if (!is.na(row$bp)){
    
    mapping <- make_mapping(row$seq, row$bp_cdna, row$bp)
    vars <- variants[variants$bp %in% mapping$bp]
    if (nrow(vars) > 0 & row$bp != '-'){
      vars$enstid <- row$enstid_version
      vars$interval <- row$bp
      vars$strand <- row$strand
      return(vars)
    } 
    return(NULL)
    
  }
  
})




# mrg1[mrg1$enstid_version %in% 'ENST00000370225.4']
result <- do.call(rbind, result_list)
res <- result[!is.na(result$bp),]
res <- res[nchar(res$ref) == 1]
res <- res[res$strand == '-',]


bool <- (lapply(1:nrow(res), function(i){
  
  row <- res[i,]
  
  s <- mrg[mrg$enstid_version %in% row$enstid]
  
  
  q <- make_mapping(s$seq, s$bp_cdna, s$bp)
  ret <- q[q$bp == row$bp,]$ref %in% row$ref
  return(ret)
  
  #mapped <- q[q$bp == row$bp,]
  #colnames(mapped)
  #original <- row
  
  
  #colnames(original)[1:4] <- paste0('clinvar.',colnames(original)[1:4])
  #colnames(original)[5:10] <- paste0('MANE.',colnames(original)[5:10])
  #colnames(mapped) <- c('mymapping.cDNA','mymapping.bp','mymapping.ref')
  #mapped$mymapping.correctly_mapped <- ret
  #result <- cbind(mapped,original)
  
  #return(result)
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

# index 59 + 1 when jumping interval

x <- "GTCAGTCAGTCACCCCACAGTCTCCTCTCTCTTCTTTTCTACTGTGCTATCCTAGAATCAAGGATTTCAGCAACA"

sum(bool)/length(bool)

res[!bool,]
index <- (1:nrow(res))[!bool]

# mrg1[mrg1$gene_symbol %in% 'AMPD2']
s <- mrg1[mrg1$gene_symbol %in% 'AMPD2']
q <- setup_mapping(s$seq, 109621175)
q[q$bp == 109621025,]
q[q$bp == 109621064,]
q[q$bp == 109621101,]
q[q$bp == 109621109,]
q[q$bp == 109621124,]

# mrg1[mrg1$gene_symbol %in% 'ABCC2']
s <- mrg1[mrg1$gene_symbol %in% 'ABCC2']
q <- setup_mapping(s$seq, 99782844)
q[q$bp == 99782756,]
q[q$bp == 99782802,]
q[q$bp == 99782805,]
q[q$bp == 99782821,]
q[q$bp == 99782822,]


## let's check direction = -1
s <- mrg1[grepl('107163510', mrg$bp)]
s <- mrg1[grepl('107867496', mrg$bp)]
make_mapping(s$seq, s$cd)
q <- setup_mapping(s$seq,  107866597) #80654490) # 80654983




q[q$bp == 108207452,]

q[q$bp == 108207498,]  # 13 113671563 A G
q[q$bp == 113671594,]  # 13 113671594 T C

q[q$bp == 114326121,] # 13 114326121 A G
q[q$bp == 114326176,] # 13 114326176 A G





q[q$bp == 80654863,] # G
q[q$bp == 80654723,] # G




q[q$bp == 80654728,] #A# match
q[q$bp == 80654729,] #T# not match
q[q$bp == 80654732,] #C# not match

q[q$bp == 80654737,]  # C


q[q$bp == 80654741,] # A
q[q$bp == 80654742,] # G
q[q$bp == 80654743,] # C

q[q$bp == 80654746,] # G
q[q$bp == 80654747,] # C
q[q$bp == 80654748,] # G
q[q$bp == 80654749,] # T
