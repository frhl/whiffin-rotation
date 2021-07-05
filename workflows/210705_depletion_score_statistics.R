library(data.table)

# compare with LOEUF
constraints <- fread('~/Projects/08_genesets/genesets/data/gnomad/karczewski2020/supplementary_dataset_11_full_constraint_metrics.tsv') #transcript-level 
constraints$loeuf <- constraints$oe_lof_upper 
constraints <- constraints[constraints$canonical == TRUE,]
constraints <- constraints[,c('gene','transcript','loeuf'), with = F]

# obs / expected
d_obs_ci <- fread('derived/210701_MANE.GRCh38.v0.95_codons_expt_ci.csv', sep = ',')
d_expt <- fread('derived/210701_MANE.GRCh38.v0.95_codons_expt.csv', sep = ',')
d_obs <- fread('derived/210701_MANE.GRCh38.v0.95_codons_obs.csv', sep = ',')

# merge observed / expected
mrg <- merge(d_expt, d_obs)
mrg$ensgid_version <- NULL
mrg$ensgid <- NULL

# observed versus expected
index_obs <- 3:66
index_expt <- 67:(67+63)

obs <- mrg[,index_expt, with = F] 
expt <- mrg[,index_obs, with = F] 

d <- (as.data.frame(colSums(obs) / colSums(expt)))
colnames(d) <- 'oe'


for (s1 in seq(0,1, by = 0.1)){
  
  print(s1)
  
  # setup score
  triplets <- rowSums(obs)
  score <- colSums(apply(obs, 1, function(x) exp(1 / (x * d$oe)))  ) / (triplets) #**s1
  
  #score <- colSums(apply(obs, 1, function(x) x * d$oe)) / (triplets) #**s1
  ds <- data.frame(ensgid = wo_version(d_obs$ensgid_version), transcript = wo_version(d_obs$enstid_version), score = (score), triplets = triplets)
  #plot(triplets, triplets**0.2)
  
  # deciles for score
  deciles_seq <- seq(0,1,by = 0.1)
  deciles <- quantile(ds$score, probs = deciles_seq, na.rm = T)
  ds$decile <- cut(ds$score, deciles)
  levels(ds$decile) <- deciles_seq*100
  
  # decile for triplets
  deciles_len_seq <- seq(0,1,by = 0.05)
  deciles_len <- quantile(ds$triplets, probs = deciles_len_seq, na.rm = T)
  ds$decile_len <- cut(ds$triplets, deciles_len)
  levels(ds$decile_len) <- deciles_len_seq*100
  compare <- merge(ds, constraints)
  compare <- compare[!is.na(compare$score) & !is.na(compare$decile),]

  
  fit <- summary(lm(log(triplets) ~ score, data = compare))$coefficients  
  print(fit)
  
  fit <- summary(lm(loeuf ~ score, data = compare))$coefficients
  print(fit)
  
  fit <- summary(lm(loeuf ~ log(triplets), data = compare))$coefficients
  print(fit)
  
  fit <- summary(lm(loeuf ~ triplets + score, data = compare))$coefficients
  print(fit)
  Sys.sleep(2)
  
}








# decile for loeuf
deciles_loeuf_seq <- seq(0,1,by = 0.1)
deciles_loeuf <- quantile(compare$loeuf, probs = deciles_loeuf_seq, na.rm = T)
compare$decile_loeuf <- cut(compare$loeuf, deciles_loeuf)
levels(compare$decile_loeuf) <- deciles_loeuf_seq*100

## plotting
compare <- compare[!is.na(compare$score) & !is.na(compare$decile),]

# depletion score
ggplot(compare, aes(x=score)) +
  geom_histogram() +
  xlab('Distribution of depletion score (Decile)') +
  ylab('Frequency') +
  ggtitle('Depletion Score') +
  geom_vline(xintercept = deciles[5], linetype = 'dashed')

# triplets versus loeuf (i.e. it's not a gene length score)
#ggplot(compare, aes(x=decile_loeuf, y=log(triplets))) +
#  geom_boxplot() +
#  xlab('gnomAD LOEUF') +
#  ylab('Triplets') +
#  ggtitle('LOEUF vs Triplet Count') 



