# constraint modelling in 5' UTR

library(data.table)
d_obs_ci <- fread('derived/210701_MANE.GRCh38.v0.95_codons_expt_ci.csv', sep = ',')
d_obs <- fread('derived/210701_MANE.GRCh38.v0.95_codons_expt.csv', sep = ',')
d_expt <- fread('derived/210701_MANE.GRCh38.v0.95_codons_obs.csv', sep = ',')

# merge observed / expected
mrg <- merge(d_expt, d_obs)
mrg$ensgid_version <- NULL
mrg$ensgid <- NULL

# observed versus expected
obs <- mrg[,3:66] 
expt <- mrg[,67:(67+63)] 

# plot data
d <- as.data.frame(colSums(obs) / colSums(expt))
colnames(d) <- 'oe'
d$codon <- unlist(lapply(strsplit(rownames(d), split = '\\.'), function(x) x[2]))
ggplot(d, aes(x=reorder(codon, oe), y = oe)) +
  geom_point() +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  ggtitle('Observed versus Expected 5" UTR codons',
          '1000 simulations preservering di-nt frequency for each sequnce') +
  ylab('Observed / Expected') +
  xlab('Codon') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# mong CG-containing codons, CGA codons for arginine are unique due to 
# their ability to create stop codons TGA (UGA in mRNA) upon epigenetic-mediated mutation


























##### ------ #####

# setup confidence interval
d_ci <- as.data.frame(matrixsplit(d_obs_ci, '\\;', as.character, index = 2))
d_ci <- as.data.frame(apply(d_ci, 2, function(x)  gsub('(\\[)|(\\])','',x)))
d_ci_expanded <- as.data.frame(do.call(cbind, lapply(colnames(d_ci)[1:64], function(col) {
  mat <- do.call(rbind, strsplit(d_ci[[col]], split = '\\~'))
  mat <- as.data.frame(apply(mat, 2, as.numeric))
  colnames(mat) <- paste0(col,c('.upper','.lower'))
  return(mat)
})))



# with a 95% confidence interval of [nobs/n0.975, nobs/n0.025]. 
sums_expt_ci <- colSums(d_ci_expanded)
sums_obs_ci <- colSums(obs[,rep(1:64, each = 2), with = F])
d_ci <- as.data.frame(sums_obs_ci / sums_expt_ci)
colnames(d_ci) <- 'oe'
d_ci$id <- names(sums_expt_ci)
d_ci$codon <- as.factor(unlist(lapply(strsplit(d_ci$id, split = '\\.'), function(x) x[2])))
d_ci$level <- as.factor(unlist(lapply(strsplit(d_ci$id, split = '\\.'), function(x) x[3])))
d_ci$id <- NULL

# 
c1 <- d_ci[d_ci$level == 'upper',]
colnames(c1) <- c('oe_upper','codon','level')
c1$level <- NULL
c2 <- d_ci[d_ci$level == 'lower',]
colnames(c2) <- c('oe_lower','codonB','level')
c2$level <- NULL
d_ci_combined <- cbind(c1, c2)
d_ci_combined$codonB <- NULL
d_ci_full <- merge(d_ci_combined, d)


ggplot(d_ci_full, aes(x=reorder(codon, oe), y = oe, ymax = oe_upper, ymin = oe_lower)) +
  geom_point() +
  geom_errorbar() +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  ggtitle('Observed versus Expected 5" UTR codons',
          '1000 simulations preservering di-nt frequency for each sequnce') +
  ylab('Observed / Expected') +
  xlab('Codon') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



#d_ci <- as.data.frame(matrixsplit(d_ci, '\\;', index = 2))


