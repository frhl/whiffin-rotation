# author: Frederik Lassen (21-06-09)
# description: 

library(ggsignif)
library(gridExtra)

# load DDG2P data
g2p <- fread('extdata/DDG2P_4_6_2021.csv')
g2p <- g2p[g2p$`DDD category` %in% c('confirmed','probable')]
colnames(g2p) <- gsub(' ','_',colnames(g2p))

# load expression and UTR complexity data
expression <- fread('derived/tables/210609_prt_rna_numerical.txt', sep = '\t')
complexity <- fread('derived/tables/210611_MANE.v0.93.UTR_features.txt', sep = '\t')
expression$bottom.percetile <- NULL
dt <- merge(complexity, expression, by.x = 'ensgid',by.y = 'gene.id')
dt$rna_std <- (dt$rna - mean(dt$rna, na.rm = T))/sd(dt$rna, na.rm = T)
dt$prt_std <- (dt$prt - mean(dt$prt, na.rm = T))/sd(dt$prt, na.rm = T)
dt$ddg2p <- as.factor(ifelse(dt$gene_symbol %in% g2p$gene_symbol, 'Y','N'))
dt$prt_rna <- dt$prt - dt$rna

# 0) 

ggplot(dt, aes(x = prt_rna, y =u5_oORF_all)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  facet_wrap(~tissue)


# 1) Do they have a higher brain expression compared to all genes?
brain_tissue <- c('Brain_Cortex','Brain_Cerebellum')
melted_expression <- melt(dt[,c('gene_symbol','tissue','ddg2p','rna_std','prt_std')])

# Look at differences in expression for G2P genes
sci <- function(x, n = 2) formatC(x, format = "e", digits = n)
cortex <- melted_expression[melted_expression$tissue == 'Brain_Cortex',]
p.values.cortex <- sci(sapply(split(cortex, cortex$variable), function(x){wilcox.test(value~ddg2p, x)$p.value}))
cerebellum <- melted_expression[melted_expression$tissue == 'Brain_Cerebellum',]
p.values.cerebellum <- sci(sapply(split(cerebellum, cortex$variable), function(x){wilcox.test(value~ddg2p, x)$p.value}))

## Brain Cortex
ggplot(melted_expression[melted_expression$tissue %in% 'Brain_Cortex',], aes(x = variable, y = value, fill = ddg2p)) +
  geom_boxplot() +
  geom_signif(y_position = c(3,5.9), xmin = c(0.8,1.8), xmax = c(1.2,2.2), annotation=c(p.values.cortex), tip_length=0.01) +
  ggtitle('Brain_Cortex RNA/Protein expression versus DDG2P genes', 'Wilcox.test') +
  ylab('Standardized Protein Expression or RNA expression') +
  xlab('Group') +
  theme_bw() +
  facet_wrap(~tissue)

## Cerebellum
ggplot(melted_expression[melted_expression$tissue %in% 'Brain_Cerebellum',], aes(x = variable, y = value, fill = ddg2p)) +
  geom_boxplot() +
  geom_signif(y_position = c(3,5.9), xmin = c(0.8,1.8), xmax = c(1.2,2.2), annotation=c(p.values.cerebellum), tip_length=0.01) +
  ggtitle('Brain_Cerebellum RNA/Protein expression versus DDG2P genes', 'Wilcox.test') +
  ylab('Standardized Protein Expression or RNA expression') +
  xlab('Group') +
  theme_bw() +
  facet_wrap(~tissue)

# slightly higher expression of DDG2P genes in the brain compared to other genes

# 2) Are there any with surprisingly low RNA/Protein expression?
dt <- merge(dt, g2p, by = 'gene_symbol')
dt <- dt[grepl("Brain/Cognition", dt$organ_specificity_list),]







cur_tissue <- 'Brain_Cortex'
cur_tissue <- 'Brain_Cerebellum'
#cur_tissue <- "Artery_Aorta"


# What is the difference of RNA/Protein distribution of mono-allelic versus bi-allelic genes? 
title = 'allelic requirement versus standardized log RNA/Protein expression'
data = dt[dt$tissue %in% cur_tissue & dt$allelic_requirement  %in% c('monoallelic','biallelic'),]
mapping = aes(x=rna_std, y = prt_std, color = allelic_requirement, label = gene_symbol)
mapping_dens = aes(x=rna_std, fill = allelic_requirement)
gg_scatter_dens(data, mapping, mapping_dens, title, cur_tissue)

# what about missense versus LoF genes?
title2 = 'missense/LoF status versus standardized log RNA/Protein expression'
data2 = dt[dt$tissue %in% cur_tissue & dt$mutation_consequence %in% c('loss of function','all missense/in frame'),]
mapping2 = aes(x=rna_std, y = prt_std, color = mutation_consequence, label = gene_symbol)
mapping_dens2 = aes(x=rna_std, fill = mutation_consequence)
gg_scatter_dens(data2, mapping2, mapping_dens2, title2, cur_tissue)

#Are there any with discordant RNA vs protein expression? Are these more likely to have UTR regulatory elements?
dt$u5_reg <- factor(ifelse(as.logical(dt$u5_oORF_inframe | dt$u5_oORF_altered_cds | dt$u5_ORF), "5-UTR Regulatory Element", "N/A"))
title3 = 'U5 regulatory elements versus standardized log RNA/Protein expression'
data3 = dt[dt$tissue %in% cur_tissue,]
mapping3 = aes(x=rna_std, y = prt_std, color = u5_reg, label = gene_symbol)
mapping_dens3 = aes(x=rna_std, fill = u5_reg)
plt <- gg_scatter_dens(data3, mapping3, mapping_dens3, title3, cur_tissue)
#ggplotly(plt)

# chisq test to compare
expression_binary <- fread('derived/tables/210609_prt_rna_binary_concordance.txt', sep = '\t')
colnames(expression_binary)[1] <- 'ensgid'
dt_binary <- merge(dt, expression_binary, by = c('ensgid','tissue'))

# RNA first
dt_rna <- dt_binary[dt_binary$value %in% c("concordance",'discordance_rna_high')]
mat_rna <- as.matrix(table(dt_rna$u5_reg, dt_rna$value))
chisq.test(mat_rna, correct = F) # 0.0009751

dt_prt <- dt_binary[dt_binary$value %in% c("concordance",'discordance_prt_high')]
mat_prt <- as.matrix(table(dt_prt$u5_reg, dt_prt$value))
chisq.test(mat_prt, correct = F) # 0.04172


## what about for only oORFs
dt_binary$u5_reg_oORF <- factor(ifelse(as.logical(dt_binary$u5_oORF_inframe | dt_binary$u5_oORF_altered_cds | dt_binary$u5_oORF_all), "5-UTR Regulatory oORF", "N/A"))

# RNA first
dt_rna <- dt_binary[dt_binary$value %in% c("concordance",'discordance_rna_high')]
mat_rna <- as.matrix(table(dt_rna$u5_reg_oORF, dt_rna$value))
chisq.test(mat_rna, correct = F) # 0.6481

dt_prt <- dt_binary[dt_binary$value %in% c("concordance",'discordance_prt_high')]
mat_prt <- as.matrix(table(dt_prt$u5_reg_oORF, dt_prt$value))
chisq.test(mat_prt, correct = F) # 0.2384


## what about for oORFs stratified by Kozak?
dt_binary$u5_reg_oORF <- factor(ifelse(as.logical(dt_binary$u5_oORF_inframe | dt_binary$u5_oORF_altered_cds | dt_binary$u5_oORF_all), "5-UTR Regulatory oORF", "N/A"))

# RNA first
dt_rna <- dt_binary[dt_binary$value %in% c("concordance",'discordance_rna_high')]
mat_rna <- as.matrix(table(dt_rna$u5_reg_oORF, dt_rna$value))
chisq.test(mat_rna, correct = F) # 0.6481

dt_prt <- dt_binary[dt_binary$value %in% c("concordance",'discordance_prt_high')]
mat_prt <- as.matrix(table(dt_prt$u5_reg_oORF, dt_prt$value))
chisq.test(mat_prt, correct = F) # 0.2384



#Do these have de novo variants in GEL NDD patients?

