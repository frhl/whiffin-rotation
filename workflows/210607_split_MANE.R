

library(data.table)

# load in data
rna <- fread('~/Projects/08_genesets/genesets/data/MANE/210607_MANE.GRCh38.v0.93.select_ensembl_rna_matrix.txt', sep = '\t')
utr <- fread('~/Projects/08_genesets/genesets/data/MANE/210607_MANE.GRCh38.v0.93.select_ensembl_genomic_matrix_UTRs.txt', sep = '\t')
cds <- fread('~/Projects/08_genesets/genesets/data/MANE/210607_MANE.GRCh38.v0.93.select_ensembl_genomic_matrix_cds.txt', sep = '\t')

# combine rna (tranccript data with utr data)
utr <- rbind(utr, cds)


#colnames(cds) <- paste0('cds.',colnames(cds))
colnames(utr) <- paste0('utr.',colnames(utr))
#mrg <- merge(cds[cds$cds.type %in% c('start_codon','stop_codon'),], rna, by.x = 'cds.enstid_version', by.y = 'enstid_version')
mrg <- merge(utr[utr$utr.type %in% c("three_prime_UTR","five_prime_UTR"),], rna, by.x = 'utr.enstid_version', by.y = 'enstid_version')
#mrg <- head(mrg)


# investigate why some ar enot matching
#bool <- ( mrg$utr.bp_start == mrg$bp_start ) | (mrg$utr.bp_end == mrg$bp_end)
#sum(bool)/length(bool)
#mrg[!bool, ]
#mrg_test = mrg[mrg$gene_symbol == 'M6PR',]
#table(mrg_test[bool,]$utr.exon)
#table(mrg_test[!bool,]$utr.exon)
# turns out they different on the exon numbers!)

# looks like 
bool1 <- (( mrg$utr.bp_start == mrg$bp_start ) | (mrg$utr.bp_end == mrg$bp_end)) 
bool2 <- mrg$utr.exon == 1
sum(bool1)/length(bool1) == sum(bool2)/length(bool2) # true when only using five_prime_UTR or three prime UTR
mrg <- mrg[bool2,]

# count uAUGs 
#mrg1 <- 
utrs <- unlist(lapply(1:nrow(mrg), function(i) {
  
  #if (mrg$utr.type[i] == 'five_prime_UTR'){
  direction = mrg$direction[i]
  if (direction == 1){
    start = mrg$bp_start[i]
    start_utr = mrg$utr.bp_start[i] - start
    end_utr = mrg$utr.bp_end[i] - start
    utr <- paste0(strsplit(mrg$seq[i], split = '')[[1]][(start_utr:end_utr)+1], collapse = '')
    return(utr)
  } else {
    end = mrg$bp_end[i]
    start_utr = abs(mrg$utr.bp_start[i] - end)
    end_utr = abs(mrg$utr.bp_end[i] - end)
    utr <- paste0(strsplit(mrg$seq[i], split = '')[[1]][(start_utr:end_utr)+1], collapse = '')
    utr_rev <- paste0(rev(unlist(strsplit(utr, split = ''))), collapse = '')
    return(utr_rev)
  }
  # }
  
  # mathces Nicky's output
  if (mrg$utr.type[i] == 'three_prime_UTR'){
    direction = mrg$direction[i]
    if (direction == 1){
      
      start = mrg$bp_start[i]
      start_utr = mrg$utr.bp_start[i] - start
      end_utr = mrg$utr.bp_end[i] - start
      utr <- paste0(strsplit(mrg$seq[i], split = '')[[1]][(start_utr:end_utr)+1], collapse = '')
      
      return(utr)
    } else {
      end = mrg$bp_end[i]
      start_utr = abs(mrg$utr.bp_start[i] - end)
      end_utr = abs(mrg$utr.bp_end[i] - end)
      utr <- paste0(strsplit(mrg$seq[i], split = '')[[1]][(start_utr:end_utr)+1], collapse = '')
      utr_rev <- paste0(rev(unlist(strsplit(utr, split = ''))), collapse = '')
      return(utr)
    }
  }
  
  
  # looks like all exons are used here.
  #if (mrg$utr.type[i] == 'CDS'){
  #  direction = mrg$direction[i]
  #  if (direction == 1){
  #    
  #    start = mrg$bp_start[i]
  #    start_utr = mrg$utr.bp_start[i] - start
  #    end_utr = mrg$utr.bp_end[i] - start
  #    utr <- paste0(strsplit(mrg$seq[i], split = '')[[1]][(start_utr:end_utr)+1], collapse = '')
  #    
  #    return(utr)
  #  } else {
  #    end = mrg$bp_end[i]
  #    start_utr = abs(mrg$utr.bp_start[i] - end)
  #    end_utr = abs(mrg$utr.bp_end[i] - end)
  #    utr <- paste0(strsplit(mrg$seq[i], split = '')[[1]][(start_utr:end_utr)+1], collapse = '')
  #    utr_rev <- paste0(rev(unlist(strsplit(utr, split = ''))), collapse = '')
  #    return(utr)
  #  }
  #}
  
  
  
}))
mrg$seq <- NULL
mrg$utr <- utrs

mrg$utr.biotype <- NULL

mrg <- mrg[,as.logical(!duplicated(t(mrg))), with = F]
colnames(mrg) <- gsub('utr\\.','',colnames(mrg))
mrg
fwrite(mrg, '~/Projects/08_genesets/genesets/data/MANE/210607_MANE.GRCh38.v0.93.5-3_prime_UTRs.txt', sep = '\t')
