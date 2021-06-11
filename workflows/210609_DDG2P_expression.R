# author: Frederik Lassen (21-06-09)
# description: 

library(ggsignif)

# load DDG2P data
g2p <- fread('extdata/DDG2P_4_6_2021.csv')
g2p <- g2p[g2p$`DDD category` %in% c('confirmed','probable')]
g2p$gene_symbol <- g2p$`gene symbol`

# load expression and UTR complexity data
expression <- fread('derived/tables/210609_prt_rna_numerical.txt', sep = '\t')
complexity <- fread('derived/tables/210611_MANE.v0.93.UTR_features.txt', sep = '\t')
expression$bottom.percetile <- NULL
dt <- merge(complexity, expression, by.x = 'ensgid',by.y = 'gene.id')
dt$rna_std <- (dt$rna - mean(dt$rna, na.rm = T))/sd(dt$rna, na.rm = T)
dt$prt_std <- (dt$prt - mean(dt$prt, na.rm = T))/sd(dt$prt, na.rm = T)
dt$ddg2p <- as.factor(ifelse(dt$gene_symbol %in% g2p$gene_symbol, 'Y','N'))



#Do they have a higher brain expression compared to all genes?
brain_tissue <- c('Brain_Cortex','Brain_Cerebellum')
melted_expression <- melt(dt[,c('gene_symbol','tissue','ddg2p','rna_std','prt_std')])


cortex <- melted_expression[melted_expression$tissue == 'Brain_Cortex',]
p.values <- sapply(split(cortex, cortex$variable), function(x){wilcox.test(ddg2p~value, x)$p.value})


ggplot(melted_expression[melted_expression$tissue %in% brain_tissue,], aes(x = variable, y = value, fill = ddg2p)) +
  geom_boxplot() +
  geom_signif(comparisons = ) +
  ylab('Standardized Protein Expression or RNA expression') +
  xlab('Group') +
  theme_bw() +
  facet_wrap(~tissue)

#Are there any with surprisingly low protein or RNA expression?

#Are there any with discordant RNA vs protein expression? Are these more likely to have UTR regulatory elements?
  
#Do these have de novo variants in GEL NDD patients?

