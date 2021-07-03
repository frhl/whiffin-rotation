devtools::load_all()

library(data.table)
library(reticulate)

d <- fread('~/Projects/08_genesets/genesets/data/MANE/210629_MANE.GRCh38.v0.95.combined-table.txt', sep = '\t')
d <- d[d$type == 'five_prime_UTR']
features <- fread('derived/tables/210629_MANE.v0.95.UTR_features.txt', sep = '\t')
ensgids <- features$ensgid[features$u5_ORF_kozak > 2]


f <- function(x) { 
  strengths <- unlist(get_kozak_strength(x))
  if (is.null(strengths)) return(0) else return(strengths)
  }

f()

enstids <- d$enstid_version

#mean_kozak <- lapply(enstids, function(cur_ensgid){print(cur_ensgid); sim_seq(d$seq[d$ensgid == cur_ensgid], f)})


table_kozak <- lapply(enstids, function(cur_enstid) table(unlist(sim_seq(d$seq[d$enstid_version == cur_enstid], f))))


mean_kozak <- lapply(enstids, function(cur_enstid) mean(unlist(sim_seq(d$seq[d$enstid_version == cur_enstid], f))))


kozaks <- list(
  strong_1 = '(A|G)..ATGG',
  moderate_1 =   '(A|G)..ATG.',
  moderate_2 =   '...ATGG',
  weak_1 =   '...ATG.',
  weak_2 = 'ATG'
)

# get all observed codons
#all_codons <- generate_codons()
codons <- unlist(kozaks)

# count observed codons
res_obs <- lapply(d$seq, function(seq){
  mat <- do.call(cbind, lapply(codons, function(codon) count_codon(seq, codon)))
  colnames(mat) <- codons
  return(mat)
})
res_mat_obs <- as.data.frame(do.call(rbind, res_obs))
colnames(res_mat_obs) <-  paste0('obs.',names(kozaks))  #paste0('obs.',codons)
res_mat_obs$ensgid_version <- d$ensgid
res_mat_obs$enstid_version <- d$enstid_version
res_mat_obs$enstid <- wo_version(res_mat_obs$enstid_version)
fwrite(res_mat_obs, 'derived/210702_MANE.GRCh38.v0.95_kozak_codons_obs.csv', sep = ',')

# simulate expected codons given sequence context
interval = TRUE
res_expt <- sim_expected_codons(d$seq[interval], k = 2, iter = 1000, codons = codons)
res_mat_expt <- as.data.frame(do.call(rbind, res_expt))
colnames(res_mat_expt) <- paste0('expt.',names(kozaks)) #codons) #paste0('expt.',codons)
res_mat_expt$ensgid <- d$ensgid[interval]
res_mat_expt$enstid_version <- d$enstid_version[interval]
res_mat_expt$enstid <- wo_version(res_mat_expt$enstid_version)
fwrite(res_mat_expt, 'derived/210702_MANE.GRCh38.v0.95_kozak_codons_expt_ci.csv', sep = ',')
mat_split_expt <- as.data.frame(matrixsplit(res_mat_expt, ';', as.numeric, 3))
mat_split_expt$ensgid <- d$ensgid[interval]
mat_split_expt$enstid_version <- d$enstid_version[interval]
mat_split_expt$enstid <- wo_version(mat_split_expt$enstid_version)
fwrite(mat_split_expt, 'derived/210702_MANE.GRCh38.v0.95_kozak_codons_expt.csv', sep = ',')







