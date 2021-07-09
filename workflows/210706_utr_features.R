# Do a GO enrichment analysis
devtools::load_all()
library(genoppi)
library(cowplot)
library(readxl)
library(ggsignif)
library(scales)

d <- fread('derived/tables/210708_MANE.v0.95.UTR_features.txt', sep = '\t')
dseq <- fread('~/Projects/08_genesets/genesets/data/MANE/210708_MANE.GRCh38.v0.95.combined-table.txt', sep = '\t')
expr <- fread('extdata/210709_table_s3_s2_combined.csv')

####################
# overall features #
####################

# " do UTRs have introns, and if so, how many"?
dseq$introns <- unlist(lapply(strsplit(dseq$bp, split = ';'), length))
count_above_9 <- sum(dseq$introns[dseq$type == 'five_prime_UTR'] > 5)
ggplot(dseq[dseq$introns < 10 & dseq$type == 'five_prime_UTR'], aes(x=introns)) +
  geom_bar(fill = 'lightgreen', color = 'black') +
  scale_x_continuous(breaks=1:9) +
  geom_vline(xintercept = 9, linetype = 'dashed') +
  annotate(geom = 'text', x = 7, y = 10000, label = '3 Genes with > 9 Introns') +
  xlab("Introns in 5' UTRs") +
  ylab('Count') +
  theme_classic()

# "How long are 5' UTRs?"
quantile(d$u5_len, probs = c(0.025, 0.975))
sum(d$u5_len > 11 & d$u5_len < 747)
ggplot(d, aes(x=(u5_len))) +
  geom_bar(fill = 'black', color = 'black', width = 0.02) +
  scale_x_log10(breaks = log_breaks()) +
  geom_vline(xintercept = 11, linetype = 'dashed', color = 'red') +
  geom_vline(xintercept = 747, linetype = 'dashed', color = 'red') +
  annotate(geom = 'text', x = 100, y = 120, label = '17.697 Genes', color = 'red') +
  annotate(geom = 'text', x = 7, y = 120, label = '2.5%', color = 'red') +
  annotate(geom = 'text', x = 1300, y = 120, label = '97.5%', color = 'red') +
  xlab("5' UTR lengths") +
  ylab('count') +
  theme_classic()

# " Are genes with uORF more liekly to be triplosensitive haploinsufficeint?"
pTS_threshold <- 0.993
pHI_threshold <- 0.84
collins2020 <- setDT(read_xlsx('~/Projects/08_genesets/genesets/data/dosage/Collins2021medrxvic_supplementary_data.xlsx', 13))
strict <- collins2020$gene[collins2020$pHI > pHI_threshold & collins2020$pTS > pTS_threshold]

dcollins <- merge(collins2020, d, by.x = 'gene', by.y = 'gene_symbol')
dcollins <- dcollins[dcollins$u5_ORF < 5]
dcollins$u5_ORF <- as.factor(dcollins$u5_ORF)
ggplot(dcollins, aes(x=u5_ORF, y = pHI, group = u5_ORF)) +
  geom_jitter(alpha = 0.05) +
  geom_boxplot(fill = 'lightblue') +
  geom_signif(comparisons = list(c("0", "1")), test = 'wilcox.test') +
  geom_hline(yintercept = pHI_threshold, linetype = 'dashed') +
  xlab('Number of uORFs') +
  ylab('Probability of haploinsufficiency (Collins 2021)') +
  #ggtitle('') +
  theme_bw() 
  

ggplot(dcollins, aes(x=u5_ORF, y = pTS, group = u5_ORF)) +
  geom_jitter(alpha = 0.05) +
  geom_boxplot(fill = 'tomato1') +
  geom_signif(comparisons = list(c("0", "1")), test = 'wilcox.test') +
  geom_hline(yintercept = pTS_threshold, linetype = 'dashed') +
  xlab('Number of uORFs') +
  ylab('Probability of triplosensitivity (Collins 2021)') +
  #ggtitle('') +
  theme_bw() 

# "Are longer UTRs more complex"?
ggplot(dcollins, aes(x=u5_len, y = pTS, group = u5_len)) +
  scale_x_log10(breaks = log_breaks()) +
  geom_point() +
  #geom_boxplot(fill = 'lightblue') +
  geom_signif(comparisons = list(c("0", "1")), test = 'wilcox.test') +
  geom_hline(yintercept = pHI_threshold, linetype = 'dashed') +
  xlab('Number of uORFs') +
  ylab('Probability of haploinsufficiency (Collins 2021)') +
  theme_bw() 


# "DO gwas loci tend to fall in the UTR"? pie-chart



#################
# GO features  #
################

# GO MF enrichment
go_mf <- fread('download/210709_hpc_derived/210709_hypergeom_go_mf_uORF_analysis.txt')
go_mf_wo_orf <- go_mf[go_mf$u5_ORF == 0]
go_mf_wo_orf <- go_mf_wo_orf[go_mf_wo_orf$pvalue < 0.0005]
go_mf_w_orf <- go_mf[go_mf$u5_ORF == 1]
go_mf_w_orf <- go_mf_w_orf[go_mf_w_orf$pvalue < 0.0005]
bonf <- 0.05 / length(unique(go_mf$list_name))

p0 <- gg_bar(go_mf_wo_orf, bonf = bonf)
p1 <- gg_bar(go_mf_w_orf, bonf = bonf)




####################
# expression wide  #
####################

# "Does multiple uORFs correspond to lowered expression?"
dexpr <- merge(d, expr)
dexpr$how <- as.factor(ifelse(grepl('RNA',dexpr$variable), 'RNA','PRT'))
dexpr$variable <- gsub('(RNA\\.)|(PRT\\.)','',dexpr$variable)
dexpr$variable <- as.factor(dexpr$variable)
dexpr <- dexpr[dexpr$u5_ORF < 5]
dexpr$u5_ORF <- factor(dexpr$u5_ORF)

ggplot(dexpr[dexpr$how == 'PRT'], aes(x=u5_ORF, y = value, group = u5_ORF)) +
  #geom_jitter(alpha = 0.05) +
  geom_boxplot(fill = 'lightblue') +
  geom_signif(comparisons = list(c("0", "1"), c("2","3")), test = 't.test') +
  xlab('Number of uORFs') +
  ylab('Standardized Protein expression (Log2 intensity)') +
  coord_flip() +
  #ggtitle('') +
  theme_bw() 

ggplot(dexpr[dexpr$how == 'RNA'], aes(x=u5_ORF, y = value, group = u5_ORF)) +
  #geom_jitter(alpha = 0.05) +
  geom_boxplot(fill = 'lightgreen') +
  geom_signif(comparisons = list(c("0", "1"), c("2","3")), test = 't.test') +
  xlab('Number of uORFs') +
  ylab('Standardized RNA expression (Log2 intensity)') +
  coord_flip() +
  #ggtitle('') +
  theme_bw() 






