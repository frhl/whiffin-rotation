# author: frederik heymann 21-06-16
# description Are genes with other organ involvement less likely to be brain specific in 
# expression? Also, are they more or less likely to have discordant RNA/protein 
# levels?

devtools::load_all()

# load DDG2P data
g2p <- fread('extdata/DDG2P_4_6_2021.csv')
g2p <- g2p[g2p$`DDD category` %in% c('confirmed','probable')] # 1809 genes
colnames(g2p) <- gsub(' ','_',colnames(g2p))

# load organ data
organ <- fread('extdata/210616_DDG2P_organ_matrix.csv')
organ <- data.frame(gene_symbol = organ$gene_symbol, 
                    organ_count = apply(organ[,-1],1 ,sum),
                    organ_pct = apply(organ[,-1],1 ,mean))

# merge data
dt_binary2 <- fread('derived/tables/210609_prt_rna_binary_concordance.txt', sep = '\t')
dt_binary2 <- merge(dt_binary2, g2p, by.x = 'hgnc_symbol', by.y = 'gene_symbol')
dt_organ <- merge(dt_binary2, organ, by.x = 'hgnc_symbol', by.y = 'gene_symbol')
dt_organ <- dt_organ[!duplicated(dt_organ),]

# Are genes with other organ involvement less likely to be brain specific in expression? 
tissues <- unique(dt_organ$tissue)


chisq <- lapply(tissues, function(x){
  mat <- as.matrix(table(dt_organ$tissue %in% x, dt_organ$organ_count > 1))
  
  sum(dt_organ$tissue %in% x & dt_organ$organ_count == 1)
  sum(dt_organ$tissue %in% x & dt_organ$organ_count > 1)
  rownames(mat) <- c('Other Tissue','Target Tissue')
  colnames(mat) <- c('1 organ','1+ organ')
  mat <- t(mat)
  mat[1, ] / colSums(mat) # 
  test <- chisq.test(mat, correct = F) # p-value < 2.2e-16
  return(test$p.value)
})
names(chisq) <- tissues
chisq <- stack(chisq)[,c(2,1)]
colnames(chisq) <- c('list_name','pvalue')
chisq$significant <- chisq$pvalue < (0.05 / 32)
chisq$logpvalue <- -log10(chisq$pvalue)

ggbarplot(chisq, 0.05/32) + 
  ylab('Tissue') + 
  ggtitle('Tissue-specificity versus DDG2P organ count',
          paste0('Probability of tissue being involved in more than 1 tissue\n',
          'chisq(table(data$tissue %in% tissue, data$organ_count > 1))'))



# Also, are they more or less likely to have discordant RNA/protein levels?
concordance <- dt_organ$organ_pct[dt_organ$value == "concordance"]
discordance_rna <- dt_organ$organ_pct[dt_organ$value == "discordance_rna_high"]
discordance_prt <- dt_organ$organ_pct[dt_organ$value == "discordance_prt_high"]
t.test(concordance, discordance_rna) # p-value = 0.5747
t.test(concordance, discordance_prt) # p-value = 0.1733

# evaluate chi square tests, are top 90% more likely to be discordant
perc90 <- quantile(dt_organ$organ_count, prob = 0.9)
dt_organ$organ_perc90 <- dt_organ$organ_count >= perc90
table(dt_organ$organ_perc90)

mat_prt <- as.matrix(table(dt_organ$value, dt_organ$organ_perc90))[-3,]
chisq.test(mat_prt, correct = F); mat_prt # p-value = 0.7613 <- this value is very sensitive to percentile
mat_rna <- as.matrix(table(dt_organ$value, dt_organ$organ_perc90))[-2,]
chisq.test(mat_rna, correct = F); mat_rna # p-value = 0.01638

# Can organ involvement be predicted by expression? 
# If so, is RNA or protein better for this? (would be neat to show that 
# this is better with protein for e.g.)




