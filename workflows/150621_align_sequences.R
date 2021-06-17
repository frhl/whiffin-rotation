setwd('~/Projects/09_whiffin_rotation/whiffin-rotation/')

# load python environment. Note, do not install minoconda. Just type NO and continue.
library(reticulate)
path <- virtualenv_python(envname = 'sandbox')
use_virtualenv(path)
os <- import("os")
os$listdir(".")
#use_virtualenv('sandbox')
source_python('python/align.py')

# load sequences
library(data.table)
sequences <- fread('~/Projects/08_genesets/genesets/data/MANE/210607_MANE.GRCh38.v0.93.combined-table.txt')
sequences <- sequences[sequences$type =='five_prime_UTR']
all(unique(ensgid) == ensgid)

# load DDG2P data
g2p <- fread('extdata/DDG2P_4_6_2021.csv')
g2p <- g2p[g2p$`DDD category` %in% c('confirmed','probable')]
colnames(g2p) <- gsub(' ','_',colnames(g2p))

# only go over DDG2P sequences
sequences <- sequences[sequences$gene_symbol %in% g2p$gene_symbol]
sequences <- head(sequences, n = 10)
ensgid <- unique(sequences$ensgid)

# parallize it
library(parallel)
library(doParallel)
cores <- detectCores()-1
registerDoParallel(cores)
nrows <- nrow(sequences)

log.socket <- make.socket(port=4002)
on.exit(close.socket(log.socket))
pprint('starting..')

mat <- (foreach (i=1:nrows, .combine=cbind) %dopar% {

  g1 <- ensgid[i]
  pprint(i)
  rows <- lapply(ensgid, function(g2){
    if (g1 != g2) {
      entry1 <- sequences[sequences$ensgid == g1,][1,]
      entry2 <- sequences[sequences$ensgid == g2,][1,]
      alignment <- align(entry1$seq, entry2$seq)
      alignment <- paste0(alignment, collapse = ':')
      return(alignment)
    } else return('')
  })
  
  rows <- do.call(rbind, rows)
  return(rows)
  
})

mat <- as.data.table(mat)
colnames(mat) <- ensgid
rownames(mat) <- ensgid
fwrite(mat, file = 'derived/tables/210615_G2P-confirmedprobable_five_prime_utr_alignment_NW.txt', sep = '\t', row.names = T, col.names = T)



# iterate through sequences and get alignment
compared <- c()
count <- 0
esnsgid <- ensgid[1:200]
mat <- lapply(ensgid[295], function(g1){
  count <<- count + 1
  progress <- paste0(count, '/', nrows)
  #progress <- paste0(round(round(100*(count/nrows)),4),collapse = '%')
  print(progress)
  
  
  
  rows <- lapply(ensgid, function(g2){
    
    comparison <- paste0(sort(c(g1, g2)), collapse = '-')
    if (!comparison %in% compared & g1 != g2) {
      entry1 <- sequences[sequences$ensgid == g1,][1,]
      entry2 <- sequences[sequences$ensgid == g2,][1,]
      alignment <- align(entry1$seq, entry2$seq)
      alignment <- paste0(alignment, collapse = ':')
      compared <<- c(compared, comparison)
      return(alignment)
    } else return('')
    
  })
  rows <- do.call(cbind, rows)
  colnames(rows) <- ensgid
  return(rows)
})

mat <- (as.data.table(do.call(rbind, cols)))
colnames(mat) <- ensgid
rownames(mat) <- ensgid
fwrite(mat, file = 'derived/tables/210615_G2P-confirmedprobable_five_prime_utr_alignment_NW.txt', sep = '\t', row.names = T, col.names = T)

