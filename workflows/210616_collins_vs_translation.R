

# description: “239/295 triplosensitive genes are also haploinsufficient” according to Collins paper. 
# These genes must be under super strict translational regulation?

collins2020 <- setDT(read_xlsx('~/Projects/08_genesets/genesets/data/dosage/Collins2021medrxvic_supplementary_data.xlsx', 13))
pTS_threshold <- 0.993
pHI_threshold <- 0.84
strict <- collins2020$gene[collins2020$pHI > pHI_threshold & 
                           collins2020$pTS > pTS_threshold]

strict

# what are the features of these genes?
# library(genoppi)
# d <- data.frame(gene = strict, significant = TRUE)
# res <- lapply_calc_hyper(d, gtex_rna, intersectN = F)
# ggbarplot(res, bonf = 0.01) 

# load translation data
expression <- fread('derived/tables/210609_prt_rna_numerical.txt', sep = '\t')
complexity <- fread('derived/tables/210615_MANE.v0.93.UTR_features.txt', sep = '\t')
dt <- merge(complexity, expression, by.x = 'ensgid',by.y = 'gene.id')
dt$transcript <- unlist(lapply(strsplit(dt$enstid_version, split = '\\.'), function(x) x[1]))
dt$rna_std <- (dt$rna - mean(dt$rna, na.rm = T))/sd(dt$rna, na.rm = T)
dt$prt_std <- (dt$prt - mean(dt$prt, na.rm = T))/sd(dt$prt, na.rm = T)
dt$ddg2p <- as.factor(ifelse(dt$gene_symbol %in% g2p$gene_symbol, 'Y','N'))
dt$prt_rna <- dt$prt - dt$rna

# all tissues
ggplot(dt[dt$gene_symbol %in% strict], aes(x=rna_std, y =prt_std)) +
  geom_point() + 
  geom_smooth(method = 'lm', formula = y ~ x) +
  facet_wrap(~tissue)

# look at pearson corss
tissue <- 'Brain_Cortex'
tissue <- 'Brain_Cerebellum'
bool <- dt$tissue %in% tissue
(cor.test(method = 'pearson', x=dt$rna[bool], y=dt$prt[bool])$estimate)
bool <- dt$tissue %in% tissue & dt$gene_symbol %in% strict
(cor.test(method = 'pearson', x=dt$rna[bool], y=dt$prt[bool])$estimate)

# somehow the strict geneset has very similar correlation

# check categorical labels
expression_binary <- fread('derived/tables/210609_prt_rna_binary_concordance.txt', sep = '\t')
colnames(expression_binary)[1] <- 'ensgid'
dt_binary <- merge(dt, expression_binary, by = c('ensgid','tissue'))
dt_binary$strict <- dt_binary$gene_symbol %in% strict
dt_binary <- dt_binary[!duplicated(dt_binary),]

# proteins in cortex are higher abundance than rna
bool <- dt_binary$tissue == 'Brain_Cortex'
(mat <- as.matrix(table(dt_binary$value[bool], dt_binary$strict[bool])))
chisq.test(mat[-2,], correct = F) # p-value = 0.384
chisq.test(mat[-3,], correct = F) # p-value = 0.02068 *

# No real difference in cerebellum
bool <- dt_binary$tissue == 'Brain_Cerebellum'
(mat <- as.matrix(table(dt_binary$value[bool], dt_binary$strict[bool])))
chisq.test(mat[-2,], correct = F) # p-value = 0.1684
chisq.test(mat[-3,], correct = F) # p-value = 0.07661



