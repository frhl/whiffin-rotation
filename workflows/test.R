setwd('~/Projects/09_whiffin_rotation/whiffin-rotation/')
devtools::load_all()

# import 5' UTR data
library(data.table)
d <- fread('~/Projects/08_genesets/genesets/data/MANE/210629_MANE.GRCh38.v0.95.combined-table.txt', sep = '\t')
d <- d[d$type == 'five_prime_UTR']
features <- fread('derived/tables/210629_MANE.v0.95.UTR_features.txt', sep = '\t')
ensgids <- features$ensgid[features$u5_AUG > 0]

# import shuffler
library(reticulate)
use_condaenv('r-reticulate')
ushuffle <- reticulate::import('ushuffle')
source_python('python/shuffle_utrs.py')


sim_expected_codons



count <- 0













