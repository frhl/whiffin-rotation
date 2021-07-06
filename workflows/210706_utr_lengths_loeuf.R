

library(data.table)
#d_obs_ci <- fread('derived/210701_MANE.GRCh38.v0.95_codons_expt_ci.csv', sep = ',')
#d_expt <- fread('derived/210701_MANE.GRCh38.v0.95_codons_expt.csv', sep = ',')
#d_obs <- fread('derived/210701_MANE.GRCh38.v0.95_codons_obs.csv', sep = ',')

df <- fread('derived/tables/210629_MANE.v0.95.UTR_features.txt', sep = '\t')

# merge observed / expected
mrg <- merge(d_expt, d_obs)
mrg$ensgid_version <- NULL
mrg$ensgid <- NULL
mrg

# observed versus expected
index_obs <- 3:66
index_expt <- 67:(67+63)

obs <- mrg[,index_expt, with = F] 
expt <- mrg[,index_obs, with = F] 


d <- data.frame(features)
d <- (as.data.frame(colSums(obs) / colSums(expt)))
colnames(d) <- 'oe'


ds <- data.frame(ensgid = wo_version(df$ensgid_version), transcript = wo_version(df$enstid_version), score = df$u5_len, triplets = df$u5_len)


#plot(ds$triplets, (ds$score), xlab = 'Triplets (Frame = 1)', ylab = 'Depletion score')

# deciles for score
deciles_seq <- seq(0,1,by = 0.1)
deciles <- quantile(ds$score, probs = deciles_seq, na.rm = T)
ds$decile <- cut(ds$score, deciles)
levels(ds$decile) <- deciles_seq*100

# decile for triplets
#deciles_len_seq <- seq(0,1,by = 0.1)
#deciles_len <- quantile(ds$triplets, probs = deciles_len_seq, na.rm = T)
#ds$decile_len <- cut(ds$triplets, deciles_len)
#levels(ds$decile_len) <- deciles_len_seq*100

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
compare <- compare[!is.na(compare$score) & !is.na(compare$decile) & !is.na(compare$decile_loeuf),]

pdf('derived/plots/210706_LENGTH_ONLY_score_summary.pdf', width = 8, height = 5)

annotate_codon <- function(oe_df, codon){
  x <- geom_vline(xintercept = oe_df$oe[oe_df$codon == codon], linetype = 'dashed') +
    annotate(geom = 'text', x = oe_df$oe[oe_df$codon == codon]+0.02, y = 3, label = codon)
  return(x)
}

# depletion score
oe_df <- data.frame(oe = oe_weights, codon = indexsplit(colnames(obs), i = 2))
p0 <- ggplot(oe_df, aes(x = oe)) +
  geom_density(fill = 'grey', color = 'black', alpha = 0.7, size = 1) +
  geom_vline(xintercept = oe_df$oe[oe_df$codon == 'ATG'], linetype = 'dashed') +
  annotate(geom = 'text', x = oe_df$oe[oe_df$codon == 'ATG']+0.02, y = 3, label = 'ATG') +
  geom_vline(xintercept = oe_df$oe[oe_df$codon == 'GCG'], linetype = 'dashed') +
  annotate(geom = 'text', x = oe_df$oe[oe_df$codon == 'GCG']-0.02, y = 3, label = 'GCG') +
  ggtitle('O/E weight distribution') + 
  theme_bw()
print(p0)

# depletion score
p1 <- ggplot(compare, aes(x=score)) +
  geom_histogram() +
  xlab('Depletion score') +
  ylab('Frequency') +
  ggtitle('Depletion score distribution and deciles') +
  geom_vline(xintercept = deciles, linetype = 'dashed', alpha = 0.3)
print(p1)

# depletion vs triplets (more triplets result in higher score)
p2 <- ggplot(compare, aes(x=decile, y=log10(triplets))) +
  geom_boxplot() +
  xlab('Depletion score (decile)') +
  ylab('log10 (triplets)') +
  ggtitle('Depletion score versus 5" UTR triplets') 
print(p2)

# triplets versus loeuf (i.e. it's not a gene length score)
p3 <- ggplot(compare, aes(x=decile_loeuf, y=log10(triplets))) +
  geom_boxplot() +
  xlab('gnomAD LOEUF') +
  ylab('log10 (triplets)') +
  ggtitle('LOEUF (gnomAD) vs Triplet Count') 
print(p3)


# gnomad Loeuf versus depletion score
test <- cor.test(compare$loeuf, compare$score, method = 'pearson')
est <- round(test$estimate, 4); pval <- formatC(test$p.value, format = "e", digits = 2)
p4 <- ggplot(compare, aes(x=score, y = loeuf)) +
  geom_point() +
  geom_density2d(color = 'grey') +
  #geom_point() +
  geom_smooth(method = 'lm') + 
  xlab('Depletion score') + 
  ylab('LOEUF') +
  ggtitle('Depletion Score versus LOEUF',
          paste0('Pearson correlation = ',est,' (P-value = ', pval,')'))
print(p4)

# gnomad Loeuf versus depletion score
test <- cor.test(log10(compare$triplets), compare$score, method = 'pearson')
est <- round(test$estimate, 4); pval <- formatC(test$p.value, format = "e", digits = 2)
p5 <- ggplot(compare, aes(x=score, y = log10(triplets))) +
  geom_point() +
  geom_density2d(color = 'grey') +
  geom_smooth(method = 'lm') + 
  xlab('Depletion score') + 
  ylab('log10 (triplets)') +
  ggtitle('Depletion Score versus triplets',
          paste0('Pearson correlation = ',est,' (P-value = ', pval,')'))
print(p5)

# triplets versus loeuf (i.e. it's not a gene length score)
p6 <- ggplot(compare, aes(x=decile_loeuf, y=score)) +
  geom_boxplot() +
  xlab('LOEUF deciles') +
  ylab('Depletion score') +
  ggtitle('LOEUF deciles versus depletion score')
print(p6)


plotit(load_mcarthur_list_loeuf('mgi_essential.tsv','mice essential genes', compare))
plotit(load_mcarthur_list_loeuf('blekhman_ar.tsv','Autosomal Rescessive (Blekhman) ', compare))
plotit(load_mcarthur_list_loeuf('clingen_level3_genes_2018_09_13.tsv','Evidence for dosage pathogenencity (Clingen, 2018)', compare))
plotit(load_mcarthur_list_loeuf('homozygous_lof_tolerant_twohit.tsv','Genes with at least two different high-confidence\nLoF variants found in a homozygous state in at least one individual in ExAC', compare))
plotit(load_mcarthur_list_loeuf('haploinsufficiency_severe_curated_2016.tsv','haploinsufficiency (Severe)', compare))
plotit(load_mcarthur_list_loeuf('NEGv1_subset_universe.tsv','Genes deemed non-essential in multiple cultured cell lines\nbased on CRISPR/Cas screen data', compare))
plotit(load_mcarthur_list_loeuf('gwascatalog.tsv','Closest gene to GWAS hits with P < 5-e8 in the NHGRI GWAS catalog', compare))
plotit(load_mcarthur_list_loeuf('olfactory_receptors.tsv', 'olfactory receptors', compare))
plotit(load_mcarthur_list_loeuf('fmrp_list_gencode.tsv','FMRP interactors', compare))

# are uAUG creating variants enriched?
features <- fread('derived/tables/210629_MANE.v0.95.UTR_features.txt', sep = '\t')
features <- features[,c('gene_symbol','ensgid_version','enstid_version')] # for mapping
aug <- fread('extdata/uAUG-creating_all_possible_annotated.txt')
aug <- aug[aug$effect == 'uORF_created',]
aug_table <- as.data.frame(table(aug$Gene, aug$Kozak_strength))
colnames(aug_table) <- c('Gene','kozak_strength','uORF_count')

# format data for plotting
mrg <- merge(features, aug_table, by.x = 'gene_symbol', by.y = 'Gene')
mrg$transcript <- wo_version(mrg$enstid_version)
dt <- merge(mrg, ds, by ='transcript')
dt1 <- aggregate(uORF_count ~ decile + kozak_strength, data = dt, FUN = sum)

ggplot(dt1, aes(x=decile, y=uORF_count, fill = kozak_strength)) +
  geom_bar(stat='identity', position = 'dodge', color = 'black') + 
  xlab('Depletion score (Decile)') +
  ylab('Frequency of uORF creating variants') +
  ggtitle('uORF creating variants (canonical Gencode)')


# split clinvar codons to see what they were before
table_aug2 <- as.data.frame(table(aug$ref, aug$alt))
colnames(table_aug2) <- c('ref','alt','freq')

aug[aug$alt == 'C',]

graphics.off()

