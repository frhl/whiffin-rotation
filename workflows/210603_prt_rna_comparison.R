# author: Frederik Heymann Lassen, 03-June-2021
# description: Investigate whether the correlation between protein and RNA expression 
# are a result of driving 

library(readxl)
library(data.table)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(ggpubr)
library(biomaRt)

# load package scripts
source('R/read_workbook.R')

# get UTR lengths
utr <- fread('extdata/biomart-5prime-utr-canonical-protein_coding.txt')
utr <- utr[!is.na(utr$`5' UTR start`) & !is.na(utr$`5' UTR end`),]
utr$UTR_length <- utr$`5' UTR end`- utr$`5' UTR start`
utr <- utr[,c(1,8)]
colnames(utr)[1] <- 'gene.id'

# load data
table_rna <- read_excel('extdata/prt_rna/table_s3.xlsx', sheet = 'B RNA tissue median', skip = 1) # C RNA standard z-score
table_prt <- read_excel('extdata/prt_rna/table_s2.xlsx', sheet = 'E protein tissue median', skip = 1) # F protein standard z-score
table_cor <- read_excel('extdata/prt_rna/table_s4.xlsx', sheet = 'B concordance comparison', skip = 1) # F protein standard z-score

# ensure that data is numeric
table_rna <- cbind(table_rna[,1], apply(table_rna[,-1],2, as.numeric))
table_prt <- cbind(table_prt[,1], apply(table_prt[,-1],2, as.numeric))

# add columns
colnames(table_rna)[-1] <- paste0('RNA.', unlist(lapply(strsplit(colnames(table_rna)[-1], split = ' '), function(x) paste0(x,collapse = '_'))))
colnames(table_prt)[-1] <- paste0('PRT.', unlist(lapply(strsplit(colnames(table_prt)[-1], split = ' '), function(x) paste0(x,collapse = '_'))))
df <- merge(table_rna, table_prt, all = T)
long <- melt(setDT(df), id.vars = 1)

# pairwise comparisons
n = 32 
pairs <- lapply(2:(n+1), function(i){
  cur_df <- df[,c(1, i, i+n), with = F]
  cur_df[[2]] <- as.numeric(cur_df[[2]])
  cur_df[[3]] <- as.numeric(cur_df[[3]])
  cols <- colnames(cur_df)[2:3]
  tissue <- unlist(strsplit(cols, split = '\\.')[2])[2]
  cur_df$tissue <- tissue
  colnames(cur_df)[2:3] <- c('rna','prt')
  cur_df <- cur_df[,c(1,4,2,3)]
  return(cur_df)
})

# check what genes are consistently found in the bottom 5th percetile
df_check <- unlist(lapply(pairs, function(x) x$gene.id[x$rna < quantile(x$rna, probs = 0.01, na.rm = T)]))
counts <- rev(sort(table(df_check)))
hist(counts)
exclude <- counts[counts > median(counts)]

# plot pearson correlations in scatter plot
pairs_df <- do.call(rbind, pairs)
pairs_df$bottom.percetile <- pairs_df$gene.id %in% names(counts)
p <- ggplot(pairs_df[pairs_df$rna != 0 & pairs_df$prt != 0,], aes(x=prt, y=rna, color = bottom.percetile), ) +
  geom_point(size = 0.5) + 
  stat_cor(method = "pearson", label.x = -2, label.y = 15, vjust = -1) + 
  geom_smooth(method='lm', formula= y~x) +
  xlab('Protein Tissue median') +
  ylab('RNA Tissue median') +
  facet_wrap(~tissue)

ggsave('derived/plots/210603_rna_prt_median_correlations_5thpercentile.pdf',width = 16, height = 16)

# check bottom quantile
for (tissue in)

cat <- unique()
cortex <- pairs_df[pairs_df$tissue %in% 'Brain_Cortex',]
cortex <- merge(cortex, utr)
cortex_narm <- cortex[!is.na(cortex$rna) & !is.na(cortex$prt),]
cor(cortex_narm$rna, cortex_narm$prt)
fit <- lm(prt ~ rna + UTR_length, data = cortex_narm)


summary(fit)
summary(lm(UTR_length ~ prt, data = cortex_narm))
summary(lm(UTR_length ~ rna, data = cortex_narm))

plot(cortex_narm$rna, log(cortex_narm$UTR_length))
summary(fit)

plot(cortex_narm$rna, log(cortex_narm$UTR_length))
plot(cortex_narm$prt, log(cortex_narm$UTR_length))


plot(residuals(fit), (cortex_narm$UTR_length), xlab = 'Resdiuals', ylab = '5"-UTR length')



# check table cor instead
head(table_cor)
dis_rna <- table_cor == 'discordance_rna_high'; sum(dis_rna, na.rm = T)
dis_prt <- table_cor == 'discordance_prt_high'; sum(dis_prt, na.rm = T)
concordance <- table_cor == 'concordance'; sum(concordance, na.rm = T)

# extract genes x categories that are disconcordant
cols <- colnames(table_cor)[grepl('category',colnames(table_cor))]
for (col in cols){
  
  
}


# compare 5' UTR lengths bween concordant / discordant across genes stratified by tissue








