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
index_obs <- 3:66
index_expt <- 67:(67+63)

obs <- mrg[,index_obs, with = F] 
expt <- mrg[,index_expt, with = F] 

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

#########################################
# ATG/GCA depletion in genes w/wo uORFS #
#########################################

features <- fread('derived/tables/210629_MANE.v0.95.UTR_features.txt', sep = '\t')

genes_w_uorf <- features$enstid_version[features$u5_ORF > 0]
genes_wo_uorf <- features$enstid_version[features$u5_ORF == 0]

# observed versus expected
obs_w_uorf <- mrg[mrg$enstid_version %in% genes_w_uorf, 3:66] 
expt_w_uorf <- mrg[mrg$enstid_version %in% genes_w_uorf, 67:(67+63)] 

obs_wo_uorf <- mrg[mrg$enstid_version %in% genes_wo_uorf, 3:66] 
expt_wo_uorf <- mrg[mrg$enstid_version %in% genes_wo_uorf, 67:(67+63)]


# get observed versus expected stratified by uORF status
d_w <- data.frame(oe = colSums(obs_w_uorf) / colSums(expt_w_uorf), uORF = 'Present')
d_w$codon <- unlist(lapply(strsplit(rownames(d_w), split = '\\.'), function(x) x[2]))
d_wo <- data.frame(oe = colSums(obs_wo_uorf) / colSums(expt_wo_uorf), uORF = 'Not present')
d_wo$codon <- unlist(lapply(strsplit(rownames(d_wo), split = '\\.'), function(x) x[2]))
d_wwo <- rbind(d_w, d_wo)

# plot data
ggplot(d_wwo, aes(x=reorder(codon, oe), y = oe, group = uORF, color = uORF)) +
  geom_point() +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  ggtitle('Observed versus Expected 5" UTR codons',
          '1000 simulations preservering di-nt frequency for each sequnce') +
  ylab('Observed / Expected') +
  xlab('Codon') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



###########################################
# ATG depletion in LOF genes after decile #
###########################################

# get LOF deciles
constraints <- fread('~/Projects/08_genesets/genesets/data/gnomad/karczewski2020/supplementary_dataset_11_full_constraint_metrics.tsv') #transcript-level 

# subset and merge
constraints$loeuf <- constraints$oe_lof_upper 
constraints <- constraints[constraints$canonical == TRUE,]

# create deciles
probs <- seq(0,1, by = 0.1)
deciles <- quantile(na.omit(constraints$loeuf),probs = probs) # only 13k/ 18k transcripts match
constraints$decile_cut <- cut(constraints$loeuf, deciles)
constraints$decile <- constraints$decile_cut
levels(constraints$decile) <- probs[-11]*100+10
constraints <- constraints[,c('gene','transcript','decile')]

m <- merge(mrg, constraints, by.x = 'enstid', by.y = 'transcript')
obs <- melt(m[,c(3:66,131:132)], id.vars = c('gene','decile'))
colnames(obs)[3:4] <- c('codon','obs')
obs$codon <- as.character(obs$codon)
obs$codon <- unlist(lapply(strsplit(obs$codon, split = '\\.'), function(x) x[2]))

expt <- melt(m[,c(67:(67+63),131:132)], id.vars = c('gene','decile'))
colnames(expt)[3:4] <- c('codon','expt')
expt$codon <- as.character(expt$codon)
expt$codon <- unlist(lapply(strsplit(expt$codon, split = '\\.'), function(x) x[2]))

# finally merge data and check
mcombi <- merge(expt, obs, by = c('gene','decile','codon'))
nrow(mcombi) == nrow(expt)
nrow(mcombi) == nrow(obs)


res <- do.call(rbind, lapply(codons, function(codon){
  m_codon <- mcombi[mcombi$codon == codon,]
  aggr1 <- aggregate(expt ~ decile, data = m_codon, FUN = sum)
  aggr2 <- aggregate(obs ~ decile, data = m_codon, FUN = sum)
  aggr12 <- merge(aggr1, aggr2)
  aggr12$codon <- codon
  return(aggr12)
}))



# plot deciles versus codon depletion
selected_res <- res[res$codon %in% c('ATG','CGA','TAA','AGT'),]
selected_res <- res[res$codon %in% head(d$codon[order(d$oe)], n = 15),]
ggplot(selected_res, aes(x=decile, y = obs/expt, color = codon, group = codon)) +
  geom_point(alpha=1) +
  geom_line(alpha=1) + 
  geom_hline(yintercept = 1, linetype = 'dashed') +
  ggtitle('Top 15 most depleted codons as a function of LOEUF deciles') +
  ylab('Observed / Expected') +
  xlab('LOEUF decile')

selected_res <- res[res$codon %in% head(d$codon[order(d$oe)], n = 15),]
ggplot(selected_res, aes(x=decile, y = obs/expt, fill = decile, group = decile)) +
  geom_boxplot() +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  ggtitle('Top 15 most depleted codons aggregated within LOEUF deciles') +
  ylab('Observed / Expected') +
  xlab('LOEUF decile')

############################################
# Mono-allelic/Bi-allelic vs ATG depletion #
############################################

# DDG2P genes
g2p <- fread('extdata/DDG2P_4_6_2021.csv')
g2p <- g2p[g2p$`DDD category` %in% c('confirmed','probable')]
colnames(g2p) <- gsub(' ','_',colnames(g2p))

genes_mono_allelic <- features$enstid_version[features$gene_symbol %in% g2p$gene_symbol[g2p$allelic_requirement %in% "monoallelic"]]
genes_bi_allelic <- features$enstid_version[features$gene_symbol %in% g2p$gene_symbol[g2p$allelic_requirement %in% "biallelic"]]

# observed versus expected
obs_mono_allelic <- mrg[mrg$enstid_version %in% genes_mono_allelic, 3:66] 
expt_mono_allelioc <- mrg[mrg$enstid_version %in% genes_mono_allelic, 67:(67+63)] 

obs_bi_allelic <- mrg[mrg$enstid_version %in% genes_bi_allelic, 3:66] 
expt_bi_allelic <- mrg[mrg$enstid_version %in% genes_bi_allelic, 67:(67+63)]

# get observed versus expected stratified by uORF status
d_mono <- data.frame(oe = colSums(obs_mono_allelic) / colSums(expt_mono_allelioc), requirement = 'mono allelic')
d_mono$codon <- unlist(lapply(strsplit(rownames(d_mono), split = '\\.'), function(x) x[2]))
d_bi <- data.frame(oe = colSums(obs_bi_allelic) / colSums(expt_bi_allelic), requirement = 'bi allelic')
d_bi$codon <- unlist(lapply(strsplit(rownames(d_bi), split = '\\.'), function(x) x[2]))
d_allelic <- rbind(d_mono, d_bi)

# plot data
ggplot(d_allelic, aes(x=reorder(codon, oe), y = oe, group = requirement, color = requirement)) +
  geom_point() +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  ggtitle('Observed versus Expected 5" UTR codons',
          '1000 simulations preservering di-nt frequency for each sequnce') +
  ylab('Observed / Expected') +
  xlab('Codon') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

###################################################### 
# Are there higher depletion in very expressed genes #
######################################################

# get protein / RNA 
expression <- fread('derived/tables/210609_prt_rna_numerical.txt', sep = '\t')
aggr_rna <- aggregate(rna ~ gene.id, data = expression, FUN = function(x) sd(x, na.rm = T))
seq_quantile <- seq(0,1, by = 0.05)
quantile_rna <- quantile(aggr_rna$rna, probs = seq_quantile)
aggr_rna$percentile <- cut(aggr_rna$rna, quantile_rna)
levels(aggr_rna$percentile) <- seq_quantile*100


merge(aggr_rna, )

obs <- mrg[,index_obs, with = F] 
expt <- mrg[,index_expt, with = F] 








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


