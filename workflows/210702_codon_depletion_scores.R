

library(data.table)
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


# setup score
triplets <- rowSums(obs)
score <- colSums(apply(obs, 1, function(x) x * d$oe)) #/ sqrt(triplets)
ds <- data.frame(ensgid = wo_version(d_obs$ensgid_version), transcript = wo_version(d_obs$enstid_version), score = log(score), triplets = triplets)

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

# compare with LOEUF
constraints <- fread('~/Projects/08_genesets/genesets/data/gnomad/karczewski2020/supplementary_dataset_11_full_constraint_metrics.tsv') #transcript-level 
constraints$loeuf <- constraints$oe_lof_upper 
constraints <- constraints[constraints$canonical == TRUE,]
constraints <- constraints[,c('gene','transcript','loeuf'), with = F]
compare <- merge(ds, constraints)


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

# depletion vs triplets (more triplets result in higher score)
ggplot(compare, aes(x=decile, y=triplets)) +
  geom_boxplot() +
  xlab('Depletion score (Decile)') +
  ylab('Triplets') +
  ggtitle('LOEUF vs Triplet Count') 

# triplets versus loeuf (i.e. it's not a gene length score)
ggplot(compare, aes(x=decile_loeuf, y=log(triplets))) +
  geom_boxplot() +
  xlab('gnomAD LOEUF') +
  ylab('Triplets') +
  ggtitle('LOEUF vs Triplet Count') 

# compare the resulting scores
ggplot(compare, aes(x=decile_len, y = decile_loeuf, fill = decile)) +
  geom_tile() +
  xlab('Depletion score (Decile)') +
  ylab('gnomAD LOEUF') +
  ggtitle('Depletion Score versus gnomAD')


# compare the resulting scores (boxplot)
ggplot(compare, aes(x=decile, y = loeuf)) +
  geom_boxplot()

ggplot(compare, aes(x=score, y = loeuf)) +
  geom_point() +
  geom_smooth(method = 'lm')


summary(lm(loeuf ~ triplets, data = compare))
summary(lm(loeuf ~ score, data = compare))
summary(lm(loeuf ~ score + triplets, data = compare)) # significant after conditioning on triplets

#ggplot(compare, aes(x=score/triplets, y = loeuf)) +
#  geom_point() +
#  geom_smooth(method = 'lm')

  #xlab('Depletion score (Decile)') +
  #ylab('gnomAD LOEUF') +
  #ggtitle('Depletion Score versus gnomAD')



# compare the resulting scores
#compare1 <- compare[,c("transcript", "loeuf","score")]
#compare2 <- compare1
#compare2$score <- compare$score / compare$triplets
#compare2$variable <- 'w triplet'
#compare1$variable <- 'w/o triplet'
#compare12 <- rbind(compare1, compare2)

#summary(lm(loeuf ~ score, data = compare1))
#summary(lm(loeuf ~ score, data = compare2))


#ggplot(compare12, aes(x = loeuf, y=score, fill = variable)) +
#  geom_boxplot() +
#  xlab('Depletion score (Decile)') +
#  ylab('gnomAD LOEUF') +
#  ggtitle('Depletion Score versus gnomAD') +
#  geom_smooth()




# 






