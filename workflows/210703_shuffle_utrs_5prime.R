#setwd('~/Projects/09_whiffin_rotation/whiffin-rotation/')
devtools::load_all()
library(data.table)
library(reticulate)

#use_condaenv('r-reticulate')
ushuffle <- reticulate::import('ushuffle')
source_python('python/shuffle_utrs.py')

# import 5' or 3' UTR data
d <- fread('../../210629_MANE.GRCh38.v0.95.combined-table.txt', sep = '\t')
#d <- fread('~/Projects/08_genesets/genesets/data/MANE/210629_MANE.GRCh38.v0.95.combined-table.txt', sep = '\t')
d <- d[d$type == 'five_prime_UTR']
features <- fread('derived/tables/210629_MANE.v0.95.UTR_features.txt', sep = '\t')
#ensgids <- features$ensgid #[features$u5_AUG > 0]


# paramters to evaluate
replicates = 50
iterations = 200
interval = TRUE
all_codons <- generate_codons()
codons <- all_codons

print('counting observed codons..')
# count observed codons
res_obs <- lapply(d$seq, function(seq){
  mat <- do.call(cbind, lapply(codons, function(codon) count_codon(seq, codon)))
  colnames(mat) <- codons
  return(mat)
})
res_mat_obs <- as.data.frame(do.call(rbind, res_obs))
colnames(res_mat_obs) <- paste0('obs.',codons)
res_mat_obs$ensgid_version <- d$ensgid_version
res_mat_obs$enstid_version <- d$enstid_version
fwrite(res_mat_obs, 'derived/210707_MANE.GRCh38.v0.95_five_prime_utr_codons_obs.csv', sep = ',')

# simulate expected codons given sequence context
write('simulating expected codons..',stdout())
write(paste('Running',replicates*iterations,'simulations.'),stdout())

res <- lapply(1:replicates, function(i){
  write(paste0(get_time(), ' - Replicate ',i),stdout())
  res_expt <- sim_expected_codons(d$seq[interval], k = 2, iter = iterations, codons = codons, parallel = T)
  
  # save confidence intervals
  outfile_ci <- paste0('derived/210707_MANE.GRCh38.v0.95_five_prime_utr_codons_expt_ci_rep',i,'.csv')
  res_mat_expt <- as.data.frame(do.call(rbind, res_expt))
  colnames(res_mat_expt) <- paste0('expt.',codons)
  res_mat_expt$ensgid_version <- d$ensgid_version[interval]
  res_mat_expt$enstid_version <- d$enstid_version[interval]
  fwrite(res_mat_expt, outfile_ci, sep = ',')
  
  # save only estimates
  outfile_est <- paste0('derived/210707_MANE.GRCh38.v0.95_five_prime_utr_codons_expt_rep',i,'.csv')
  mat_split_expt <- as.data.frame(matrixsplit(res_mat_expt, ';', as.numeric, 3))
  mat_split_expt$ensgid_version <- d$ensgid_version[interval]
  mat_split_expt$enstid_version <- d$enstid_version[interval]
  fwrite(mat_split_expt, outfile_est, sep = ',')
  
})


