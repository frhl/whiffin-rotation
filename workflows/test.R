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


mrg <- merge(res_mat_obs, mat_split_expt)
mrg$ensgid_version <- NULL
mrg$ensgid <- NULL

obs <- mrg[,3:8] 
expt <- mrg[,9:14] 

sum(obs$obs.ATG) / sum(expt$expt.ATG)
sum(obs$obs.TAG) / sum(expt$expt.TAG)


d <- as.data.frame(colSums(obs) / colSums(expt))
colnames(d) <- 'oe'
d$codon <- unlist(lapply(strsplit(rownames(d), split = '\\.'), function(x) x[2]))

ggplot(d, aes(x=codon, y = oe)) +
  geom_point() +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  ggtitle('Observed versus Expected 5" UTR codons') +
  ylab('Observed / Expected')


oe_mat <- obs / expt


oe_mat[is.na(oe_mat)] <- 0
oe_mat$enstid <- mrg$enstid


library(reshape2)

melted <- melt(oe_mat)
ggplot(melted, aes(x=variable, y = value)) +
  geom_violin() +
  ylim(0,3)


mrg














