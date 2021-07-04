

library(data.table)

d <- fread('~/Projects/08_genesets/genesets/data/MANE/210629_MANE.GRCh38.v0.95.combined-table.txt', sep = '\t')
d <- d[d$type == 'five_prime_UTR']
features <- fread('derived/tables/210629_MANE.v0.95.UTR_features.txt', sep = '\t')
ensgids <- features$ensgid[features$u5_ORF_kozak > 2]


# extract all uORFS
d$ensgid
d_uorf <- lapply(d$ensgid, function(x) get_orf(d$seq[d$ensgid == x], share_stops = F))
names(d_uorf) <- d$ensgid
d_mat <- stack(d_uorf)
colnames(d_mat) <- c('seq','ensgid')
d_copy <- d
d_copy$seq <- NULL
d_mat <- merge(d_mat, d_copy)

# get all observed codons
all_codons <- generate_codons()
codons <- all_codons #c('TAA')

# count observed codons
res_obs <- lapply(d_mat$seq, function(seq){
  mat <- do.call(cbind, lapply(codons, function(codon) count_codon(seq, codon)))
  colnames(mat) <- codons
  return(mat)
})
res_mat_obs <- as.data.frame(do.call(rbind, res_obs))
colnames(res_mat_obs) <- paste0('obs.',codons)
res_mat_obs$ensgid_version <- d_mat$ensgid
res_mat_obs$enstid_version <- d_mat$enstid_version
res_mat_obs$enstid <- wo_version(res_mat_obs$enstid_version)
res_mat_obs$ensgid <- wo_version(as.character(res_mat_obs$ensgid_version))
fwrite(res_mat_obs, 'derived/210702_MANE.GRCh38.v0.95_u5orf_codons_obs.csv', sep = ',')

# simulate expected codons given sequence context
interval = TRUE
res_expt <- sim_expected_codons(d_mat$seq[interval], k = 2, iter = 1000, codons = codons)
res_mat_expt <- as.data.frame(do.call(rbind, res_expt))
colnames(res_mat_expt) <- paste0('expt.',codons)
res_mat_expt$ensgid_version <- d_mat$ensgid_version[interval]
res_mat_expt$enstid_version <- d_mat$enstid_version[interval]
res_mat_expt$enstid <- wo_version(res_mat_expt$enstid_version)
fwrite(res_mat_expt, 'derived/210702_MANE.GRCh38.v0.95_u5orf_codons_expt_ci.csv', sep = ',')
mat_split_expt <- as.data.frame(matrixsplit(res_mat_expt, ';', as.numeric, 3))
mat_split_expt$ensgid <- d_mat$ensgid[interval]
mat_split_expt$enstid_version <- d_mat$enstid_version[interval]
mat_split_expt$enstid <- wo_version(mat_split_expt$enstid_version)
fwrite(mat_split_expt, 'derived/210702_MANE.GRCh38.v0.95_u5orf_codons_expt.csv', sep = ',')
