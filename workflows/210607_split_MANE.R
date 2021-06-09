

library(data.table)

# load in data
rna <- fread('~/Projects/08_genesets/genesets/data/MANE/210607_MANE.GRCh38.v0.93.select_ensembl_rna_matrix.txt', sep = '\t')
utr <- fread('~/Projects/08_genesets/genesets/data/MANE/210607_MANE.GRCh38.v0.93.select_ensembl_genomic_matrix_UTRs.txt', sep = '\t')
cds <- fread('~/Projects/08_genesets/genesets/data/MANE/210607_MANE.GRCh38.v0.93.select_ensembl_genomic_matrix_cds.txt', sep = '\t')

# combine rna (tranccript data with utr data)
utr <- rbind(utr, cds)

#colnames(cds) <- paste0('cds.',colnames(cds))
colnames(utr) <- paste0('utr.',colnames(utr))
mrg <- merge(utr[utr$utr.type %in% c('CDS',"exon","three_prime_UTR","five_prime_UTR"),], rna, by.x = 'utr.enstid_version', by.y = 'enstid_version')


#tbl <- table(mrg$gene_symbol, mrg$utr.type) # there can be more than one five prime UTR for the same gene
#mrg <- mrg[mrg$gene_symbol == "ACAN", ]#'DDX1',]
#mrg <- mrg[mrg$gene_symbol == 'M6PR', ]
#mrg <- mrg[mrg$gene_symbol == 'DDX1',]
#mrg$seq <- NULL


enstids <- unique(mrg$utr.enstid_version)
enstid <- enstids 

revseq <- function(x) paste0(rev(unlist(strsplit(x, ''))), collapse = '')


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

  
  # setup vars
  df$lens <- df$utr.bp_end - df$utr.bp_start
  
  # check that things matches up
  stopifnot(sum(df$utr.bp_end - df$utr.bp_start) == nchar(seq)- nrow(df))
  
  # sequence handling
  start <- unique(df$bp_start)
  end <- unique(df$bp_end)
  lens <- 0
  for (i in 1:nrow(df)){
    
    # current position with repsect to genomic sequence
    x1 <- df$utr.bp_start[i] - start
    x2 <- df$utr.bp_end[i] - start
    
    # current position with respect to mRNA
    intron_len <- x1-lens
    x1 <- x1 - intron_len
    x2 <- x2 - intron_len
    
    # keep track of interval 
    interval <- (x1:x2)+(i)
    newseq <- paste0(unlist(strsplit(seq,split=''))[interval], collapse = '')
    df$newseq[i] <- ifelse(direction == 1, newseq, revseq(newseq))
    
    lens <- lens + df$lens[i]
  }
  
  return(df)
  
})

# combne data
seq_df <- do.call(rbind, sequences)
seq_df <- seq_df[,as.logical(!duplicated(t(seq_df))), with = F]
colnames(seq_df) <- gsub('utr\\.','',colnames(seq_df))
fwrite(seq_df, '~/Projects/08_genesets/genesets/data/MANE/210607_MANE.GRCh38.v0.93.5-3_all_exon_seqs.txt', sep = '\t')
seq_df <- fread('~/Projects/08_genesets/genesets/data/MANE/210607_MANE.GRCh38.v0.93.5-3_all_exon_seqs.txt', sep = '\t')
# subset data




sequences_combined <- do.call(rbind, lapply(enstids, function(enstid){
  
  # enstid <- "ENST00000255082.8"
  bool_transcript <- seq_df$enstid_version == enstid
  df <- seq_df[bool_transcript, ]
  direction <- unique(df$direction)
  if (direction == '-'){
    df <- df[nrow(df):1]
  }
  
  
  UTR_5 <- paste0(df$newseq[df$type == 'five_prime_UTR'], collapse = '')
  CDS <- paste0(df$newseq[df$type == 'CDS'], collapse = '')
  UTR_3 <- paste0(df$newseq[df$type == 'three_prime_UTR'], collapse = '')
  
  row <- df[1,c(10,11,16, 1, 12)]
  r1 <- data.frame(row, type = 'five_prime_UTR', seq = UTR_5)
  r2 <- data.frame(row, type = 'three_prime_UTR', seq = UTR_3)
  r3 <- data.frame(row, type = 'CDS', seq = CDS)
  
  res <- rbind(r1, r2, r3)
  return(res)
  
}))

fwrite(sequences_combined, '~/Projects/08_genesets/genesets/data/MANE/210607_MANE.GRCh38.v0.93.combined-table.txt', sep = '\t')
