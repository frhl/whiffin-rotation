# author: Frederik Lassen (21-06-09)
# description: 

library(ggsignif)
library(gridExtra)
library(plotly)

# load DDG2P data
g2p <- fread('extdata/DDG2P_4_6_2021.csv')
g2p <- g2p[g2p$`DDD category` %in% c('confirmed','probable')]
colnames(g2p) <- gsub(' ','_',colnames(g2p))

# load constraint data
#constraint <- fread('~/Projects/08_genesets/genesets/data/gnomad/karczewski2020/supplementary_dataset_11_full_constraint_metrics.tsv')

# load expression and UTR complexity data
expression <- fread('derived/tables/210609_prt_rna_numerical.txt', sep = '\t')
complexity <- fread('derived/tables/210615_MANE.v0.93.UTR_features.txt', sep = '\t')
dt <- merge(complexity, expression, by.x = 'ensgid',by.y = 'gene.id')
dt$transcript <- unlist(lapply(strsplit(dt$enstid_version, split = '\\.'), function(x) x[1]))
dt$rna_std <- (dt$rna - mean(dt$rna, na.rm = T))/sd(dt$rna, na.rm = T)
dt$prt_std <- (dt$prt - mean(dt$prt, na.rm = T))/sd(dt$prt, na.rm = T)
dt$ddg2p <- as.factor(ifelse(dt$gene_symbol %in% g2p$gene_symbol, 'Y','N'))
dt$prt_rna <- dt$prt - dt$rna

# 0) 
#ggplot(dt, aes(x = prt_rna, y =u5_oORF_all)) +
#  geom_point() +
#  geom_vline(xintercept = 0, linetype = 'dashed') +
#  facet_wrap(~tissue)


# 1) Do they have a higher brain expression compared to all genes?
brain_tissue <- c('Brain_Cortex','Brain_Cerebellum')
melted_expression <- melt(dt[,c('gene_symbol','tissue','ddg2p','rna_std','prt_std')])

# Look at differences in expression for G2P genes
sci <- function(x, n = 2) formatC(x, format = "e", digits = n)
cortex <- melted_expression[melted_expression$tissue == 'Brain_Cortex',]
p.values.cortex <- sci(sapply(split(cortex, cortex$variable), function(x){wilcox.test(value~ddg2p, x)$p.value}))
cerebellum <- melted_expression[melted_expression$tissue == 'Brain_Cerebellum',]
p.values.cerebellum <- sci(sapply(split(cerebellum, cerebellum$variable), function(x){wilcox.test(value~ddg2p, x)$p.value}))
intestine <- melted_expression[melted_expression$tissue == 'Small_Intestine',]
p.values.intestine <- sci(sapply(split(intestine, intestine$variable), function(x){wilcox.test(value~ddg2p, x)$p.value}))


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
ggplot(dt[dt$tissue %in% brain_tissue,], aes(x = rna_std, y = prt_std, color = ddg2p)) +
  geom_point(data = dt[dt$tissue %in% brain_tissue & dt$ddg2p == 'N',]) +
  geom_point(data = dt[dt$tissue %in% brain_tissue & dt$ddg2p == 'Y',]) +
  ggtitle('Brain_Cerebellum/Cortex RNA/Protein expression of DDG2P genes') +
  ylab('Standardized Protein expression') +
  xlab('Standardized RNA expression') +
  geom_smooth() +
  theme_bw() +
  facet_wrap(~tissue)

data = dt[dt$tissue %in% 'Brain_Cortex',]
mapping = aes(x=rna_std, y = prt_std, color = ddg2p, label = gene_symbol)
mapping_dens = aes(x=rna_std, fill = ddg2p)
gg_scatter_dens(data, mapping, mapping_dens, title, cur_tissue, debug = T)



dt <- merge(dt, g2p, by = 'gene_symbol')
#dt <- merge(dt, constraint, by = 'transcript') 
dt <- dt[grepl("Brain/Cognition", dt$organ_specificity_list),]


cur_tissue <- 'Brain_Cortex'
cur_tissue <- 'Brain_Cerebellum'
#cur_tissue <- "Artery_Aorta"
#cur_tissue <- "Breast"

pdf('derived/plots/210617_scatter_density_g2p_versus_rna-prt_expression.pdf', width = 10, height = 12)
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
ggplotly(plt)
graphics.off()

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

# RNA first
dt_rna <- dt_binary[dt_binary$value %in% c("concordance",'discordance_rna_high') ]
mat_rna <- as.matrix(table(dt_rna$u5_oORF_kozak, dt_rna$value))
chisq.test(mat_rna, correct = F) # 0.6481
DescTools::CochranArmitageTest(t(mat_rna), 'increasing')
DescTools::Desc(t(mat_rna))

dt_prt <- dt_binary[dt_binary$value %in% c("concordance",'discordance_prt_high') ]
mat_prt <- as.matrix(table(dt_prt$u5_oORF_kozak, dt_prt$value))
DescTools::CochranArmitageTest(t(mat_prt), 'increasing')
chisq.test(mat_prt, correct = F) # 0.2384

## what about for uORFs stratified by Kozak?

# RNA first 5' UTR (armitage trend test)
dt_rna <- dt_binary[dt_binary$value %in% c("concordance",'discordance_rna_high') ]
mat_rna <- t(as.matrix(table(dt_rna$u5_ORF_kozak, dt_rna$value)))
mat_rna[1,]/colSums(mat_rna)
DescTools::CochranArmitageTest(mat_rna, 'increasing')
#chisq.test(mat_rna, correct = F) # 0.00038 ***

dt_prt <- dt_binary[dt_binary$value %in% c("concordance",'discordance_prt_high') ]
mat_prt <- t(as.matrix(table(dt_prt$u5_ORF_kozak, dt_prt$value)))
DescTools::CochranArmitageTest(mat_prt, 'increasing')
#chisq.test(mat_prt, correct = F) # 0.01555 *

# RNA first 3' UTR (armitage trend test)
#dt_rna <- dt_binary[dt_binary$value %in% c("concordance",'discordance_rna_high')] #& dt$u3_ORF_kozak %in% c(0,3)]
#mat_rna <- t(as.matrix(table(dt_rna$u3_ORF_kozak, dt_rna$value)))
#DescTools::CochranArmitageTest(mat_rna, 'increasing')
#DescTools::Desc(mat_rna)
#hisq.test(mat_rna, correct = F) # 0.00038 ***

#dt_prt <- dt_binary[dt_binary$value %in% c("concordance",'discordance_prt_high') ]
#mat_prt <- t(as.matrix(table(dt_prt$u3_ORF_kozak, dt_prt$value)))
#DescTools::CochranArmitageTest(mat_prt, 'increasing')
#DescTools::Desc(mat_prt)
#chisq.test(mat_prt, correct = F) # 0.01555 *


#Do these have de novo variants in GEL NDD patients?




Desc(mat)






#Can organ involvement be predicted by expression? If so, is RNA or protein better for this? (would be neat to show that this is better with protein for e.g.)

