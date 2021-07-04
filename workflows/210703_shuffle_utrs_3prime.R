setwd('~/Projects/09_whiffin_rotation/whiffin-rotation/')
devtools::load_all()

# import 5' or 3' UTR data
library(data.table)
d <- fread('~/Projects/08_genesets/genesets/data/MANE/210629_MANE.GRCh38.v0.95.combined-table.txt', sep = '\t')
#d <- d[d$type == 'five_prime_UTR']
d <- d[d$type == 'three_prime_UTR']
features <- fread('derived/tables/210629_MANE.v0.95.UTR_features.txt', sep = '\t')
ensgids <- features$ensgid[features$u5_AUG > 0]

# import shuffler
library(reticulate)
use_condaenv('r-reticulate')
ushuffle <- reticulate::import('ushuffle')
source_python('python/shuffle_utrs.py')


# helper
wo_version <- function(x) unlist(lapply(strsplit(x, split = '\\.'), function(x) x[1]))


all_codons <- generate_codons()

# test

# get all observed codons
codons <- c('ATG')#all_codons

# count observed codons
res_obs <- lapply(d$seq, function(seq){
  mat <- do.call(cbind, lapply(codons, function(codon) count_codon(seq, codon)))
  colnames(mat) <- codons
  return(mat)
})
res_mat_obs <- as.data.frame(do.call(rbind, res_obs))
colnames(res_mat_obs) <- paste0('obs.',codons)
res_mat_obs$ensgid_version <- d$ensgid
res_mat_obs$enstid_version <- d$enstid_version
res_mat_obs$enstid <- wo_version(res_mat_obs$enstid_version)
fwrite(res_mat_obs, 'derived/210701_MANE.GRCh38.v0.95_three_prime_utr_codons_obs.csv', sep = ',')

# simulate expected codons given sequence context
interval = TRUE
res_expt <- sim_expected_codons(d$seq[interval], k = 2, iter = 1000, codons = codons)
res_mat_expt <- as.data.frame(do.call(rbind, res_expt))
colnames(res_mat_expt) <- paste0('expt.',codons)
res_mat_expt$ensgid <- d$ensgid[interval]
res_mat_expt$enstid_version <- d$enstid_version[interval]
res_mat_expt$enstid <- wo_version(res_mat_expt$enstid_version)
fwrite(res_mat_expt, 'derived/210701_MANE.GRCh38.v0.95_three_prime_utr_codons_expt_ci.csv', sep = ',')
mat_split_expt <- as.data.frame(matrixsplit(res_mat_expt, ';', as.numeric, 3))
mat_split_expt$ensgid <- d$ensgid[interval]
mat_split_expt$enstid_version <- d$enstid_version[interval]
mat_split_expt$enstid <- wo_version(mat_split_expt$enstid_version)
fwrite(mat_split_expt, 'derived/210701_MANE.GRCh38.v0.95_three_prime_utr_codons_expt.csv', sep = ',')

# plot all 
#res_prob <- sim_prob_codons(d$seq[interval], k = 2, iter = 1000, codons = codons)
##res_mat_prob <- as.data.frame(do.call(rbind, res_prob))
#colnames(res_mat_prob) <- paste0('prob.',codons)
#res_mat_prob$ensgid <- d$ensgid[interval]
#res_mat_prob$enstid_version <- d$enstid_version[interval]
#res_mat_prob$enstid <- wo_version(res_mat_prob$enstid_version)
#fwrite(res_mat_prob, 'derived/210701_MANE.GRCh38.v0.95_codons_probs.csv', sep = ',')


#res_mat <- merge(res_mat , features)
#res_mat$delta <- abs(res_mat$prob.ATG - as.numeric(res_mat$u5_AUG > 0))
#res_mat$transcript <- wo_version(res_mat$enstid_version)
#mrg <- merge(res_mat, constraints, by = 'transcript')
#mrg <- mrg[mrg$delta != 0,]
#plot(mrg$delta, mrg$pLI)



