

library(data.table)

# load in data
rna <- fread('~/Projects/08_genesets/genesets/data/MANE/210629_MANE.GRCh38.v0.95.select_ensembl_rna_matrix.txt', sep = '\t')
utr <- fread('~/Projects/08_genesets/genesets/data/MANE/210629_MANE.GRCh38.v0.95.select_ensembl_genomic_matrix_UTRs.txt', sep = '\t')
cds <- fread('~/Projects/08_genesets/genesets/data/MANE/210629_MANE.GRCh38.v0.95.select_ensembl_genomic_matrix_cds.txt', sep = '\t')

# combine rna (tranccript data with utr data)
utr <- rbind(utr, cds)

# merge data
colnames(utr) <- paste0('utr.',colnames(utr))
mrg <- merge(utr[utr$utr.type %in% c('CDS',"exon","three_prime_UTR","five_prime_UTR"),], rna, by.x = 'utr.enstid_version', by.y = 'enstid_version')


# Prep for writing out the data
enstids <- unique(mrg$utr.enstid_version)
enstid <- enstids 

# helper functions
revseq <- function(x) paste0(rev(unlist(strsplit(x, ''))), collapse = '')

# clinvar matching
#clinvar <- fread('extdata/clinvar/clinvar_20210626_chr1.txt')
#colnames(clinvar) <- c('chr','bp','ref','alt')
#mymapping <- list()

sequences <- lapply(enstids, function(enstid){
  
  # enstid <- "ENST00000255082.8"
  bool_transcript <- mrg$utr.enstid_version == enstid
  df <- mrg[bool_transcript, ]
  df <- df[order(df$utr.bp_start,)]
  df <- df[df$utr.type != 'exon',]
  direction <- unique(df$direction)
  seq <- unique(df$seq)
  seq <- ifelse(direction == 1, seq, revseq(seq))
  df$seq <- NULL
  df$newseq <- NA
  df$cdna <- NA
  df$cdna_rev <- NA

  # setup vars
  df$lens <- df$utr.bp_end - df$utr.bp_start
  
  # check that exons + introns matches up to total gene length
  stopifnot(sum(df$utr.bp_end - df$utr.bp_start) == nchar(seq)- nrow(df))
  
  # sequence handling
  start <- unique(df$bp_start)
  end <- unique(df$bp_end)
  lens <- 0

  for (i in 1:nrow(df)){
    
      # current position (cDNA)
      x1 <- df$utr.bp_start[i] - start
      x2 <- df$utr.bp_end[i] - start
      
      # remove introns
      intron_len <- x1-lens
      x1 <- abs(x1 - intron_len)
      x2 <- abs(x2 - intron_len)
      
      # position in gene
      stopifnot(length(df$utr.bp_start[i]:df$utr.bp_end[i]) == length(x1:x2))
      
      # keep track of interval 
      #interval_rev <- (y1:y2)+(i)
      interval <- (x1:x2)+(i)
      newseq <- paste0(unlist(strsplit(seq,split=''))[interval], collapse = '')
      df$newseq[i] <- ifelse(direction == 1, newseq, revseq(newseq))
      df$cdna[i] <- paste0(x1+i,'-',x2+i)
      
      
      #df$cdna_rev[i] <- paste0(y1+i,'-',y2+i)
      
      
      lens <- lens + df$lens[i]
      
  }
  
  
  ## NOTE: THIS IS AN UNGLY HACK TO AVOID DEALING WITH REVERSE STRAND
  # IN A PROPER WAY. TODO, AVOID THIS MESS
  if (direction == -1){
    lens <- 0
    lst <- list()
    for (i in nrow(df):1){
      
      # current position (cDNA)
      x1 <- df$utr.bp_start[i] - start
      x2 <- df$utr.bp_end[i] - start
      
      # remove introns
      intron_len <- x1-lens
      x1 <- abs(x1 - intron_len)
      x2 <- abs(x2 - intron_len)

      # keep track of interval 
      u <- nrow(df)-i+1
      lst[[i]] <- (paste0(x1+u,'-',x2+u))
      lens <- lens + df$lens[i]
      
    }
    
    df$cdna <- unlist(lst)
    
    
  }
  
  return(df)
  
})



# combine data
seq_df <- do.call(rbind, sequences)
seq_df <- seq_df[,as.logical(!duplicated(t(seq_df))), with = F]
colnames(seq_df) <- gsub('utr\\.bp_','exon\\.bp_',colnames(seq_df))
colnames(seq_df) <- gsub('utr\\.','',colnames(seq_df))
#fwrite(seq_df, '~/Projects/08_genesets/genesets/data/MANE/210710_MANE.GRCh38.v0.95.5-3_all_exon_seqs.txt', sep = '\t')
#seq_df <- fread('~/Projects/08_genesets/genesets/data/MANE/210710_MANE.GRCh38.v0.95.5-3_all_exon_seqs.txt', sep = '\t')


sequences_combined <- do.call(rbind, lapply(enstids, function(enstid){
  
  # enstid <- "ENST00000255082.8"
  bool_transcript <- seq_df$enstid_version == enstid
  df <- seq_df[bool_transcript, ]
  direction <- unique(df$direction)
  stopifnot(length(direction) == 1)
  if (direction == '-'){
    df <- df[nrow(df):1]
  }
  
  # get RNA sequence
  UTR_5 <- paste0(df$newseq[df$type == 'five_prime_UTR'], collapse = '')
  CDS <- paste0(df$newseq[df$type == 'CDS'], collapse = '')
  UTR_3 <- paste0(df$newseq[df$type == 'three_prime_UTR'], collapse = '')
  
  # Keep cDNA positions
  UTR_5_cdna <- paste0(df$cdna[df$type == 'five_prime_UTR'], collapse = ';')
  CDS_cdna <- paste0(df$cdna[df$type == 'CDS'], collapse = ';')
  UTR_3_cdna <- paste0(df$cdna[df$type == 'three_prime_UTR'], collapse = ';')  
  
  # Keep Grch38 positions
  UTR_5_bp <- paste0(df$exon.bp_start[df$type == 'five_prime_UTR'], '-',df$exon.bp_end[df$type == 'five_prime_UTR'], collapse = ';')
  CDS_bp <- paste0(df$exon.bp_start[df$type == 'CDS'], '-',df$exon.bp_end[df$type == 'CDS'], collapse = ';')
  UTR_3_bp <- paste0(df$exon.bp_start[df$type == 'three_prime_UTR'], '-',df$exon.bp_end[df$type == 'three_prime_UTR'], collapse = ';')  
  chrom = unique(df$chr)
  stopifnot(length(chrom) == 1)
  
  row <- df[1,c(10,11,16, 1, 12)]
  r1 <- data.frame(row, type = 'five_prime_UTR', chr = chrom, bp = UTR_5_bp, bp_cdna = UTR_5_cdna, strand = direction, seq = UTR_5)
  r2 <- data.frame(row, type = 'three_prime_UTR', chr = chrom, bp = CDS_bp,  bp_cdna = CDS_cdna, strand = direction, seq = UTR_3)
  r3 <- data.frame(row, type = 'CDS', chr = chrom, bp = UTR_3_bp, bp_cdna = UTR_3_cdna, strand = direction, seq = CDS)
  
  res <- rbind(r1, r2, r3)
  return(res)
  
}))

#test <- sequences_combined
#test <- test[test$enstid_version == "ENST00000361915.8" & test$type == 'five_prime_UTR',]
#mapping <- make_mapping(test$seq, test$bp_cdna, test$bp)
#merge(mapping, clinvar)


fwrite(sequences_combined, '~/Projects/08_genesets/genesets/data/MANE/210710_MANE.GRCh38.v0.95.combined-table.txt', sep = '\t')
