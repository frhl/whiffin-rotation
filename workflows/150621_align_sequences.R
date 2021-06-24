setwd('~/Projects/09_whiffin_rotation/whiffin-rotation/')

# load python environment. Note, do not install minoconda. Just type NO and continue.
library(reticulate)
path <- virtualenv_python(envname = 'sandbox')
use_virtualenv(path)
use_virtualenv('sandbox')
source_python('python/align.py')
os <- import("os")
os$listdir(".")


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
sequences <- head(sequences, n = 1000)
ensgid <- unique(sequences$ensgid)


# tri(500x500) = 124.750  iter: 1.5h

# setup iter grid
itermat <- melt_tri(nrow(sequences), how = 'upper', diagonal = F)

# parallize it
library(parallel)
library(doParallel)
cores <- detectCores()
registerDoParallel(cores)
nrows <- nrow(sequences)

#log.socket <- make.socket(port=4001)
#on.exit(close.socket(log.socket))
pprint('starting..')
pprint(paste0('iterations:',sum(itermat$iter)))

mat <- (foreach (i=1:nrows, .combine=cbind) %dopar% {

  #g1 <- ensgid[i]
  pprint(i)
  rows <- lapply(1:nrows, function(j){
    
    allow_iter <- itermat$iter[itermat$i == i & itermat$j == j]
    if (allow_iter){
      
      # get genes
      g1 <- ensgid[i]
      g2 <- ensgid[j]
      
      # get sequences
      entry1 <- sequences[sequences$ensgid == g1,][1,]
      entry2 <- sequences[sequences$ensgid == g2,][1,]
      
      # align using python
      alignment <- align(entry1$seq, entry2$seq)
      alignment <- paste0(alignment, collapse = ':')
      return(alignment)
    } else return('')
  })
  
  rows <- do.call(rbind, rows)
  return(rows)
  
})

tmat <- as.data.table(mat)
colnames(tmat) <- ensgid
rownames(tmat) <- ensgid
fwrite(tmat, file = 'derived/tables/210623_1000_G2P-confirmedprobable_five_prime_utr_alignment_NW.txt', sep = '\t', row.names = T, col.names = T)


