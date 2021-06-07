# author: Frederik Heymann Lassen, 02-June-2021
# description: investigate upper-bound LoF fraction of O/E protein trunctaing variation
# and haploinsifficient/triplosensitive genes

library(readxl)
library(data.table)
library(ggplot2)
library(ggrepel)
library(cowplot)

# load package scripts
source('R/read_workbook.R')

# load data
collins2020 <- setDT(read_xlsx('~/Projects/08_genesets/genesets/data/dosage/Collins2021medrxvic_supplementary_data.xlsx', 13))
constraints <- fread('~/Projects/08_genesets/genesets/data/gnomad/karczewski2020/supplementary_dataset_11_full_constraint_metrics.tsv') #transcript-level 
olink <- read_workbook('extdata/Olink - Gene List Comparison - Tissue specific essential genes.xlsx')

# subset and merge
constraints$loeuf <- constraints$oe_lof_upper 
constraints <- constraints[constraints$canonical == TRUE,]
constraints <- merge(constraints, collins2020, all = T)

# create deciles
probs <- seq(0,1, by = 0.1)
deciles <- quantile(na.omit(constraints$loeuf),probs = probs)
constraints$decile_cut <- cut(constraints$loeuf, deciles)
constraints$decile <- constraints$decile_cut
constraints$decile_plot <- constraints$decile_cut
levels(constraints$decile) <- probs[-11]*100+10

# create pHI groups and pTS group
constraints$pTSgroup <- cut(constraints$pTS, probs)
constraints$pHIgroup <- cut(constraints$pHI, probs)
levels(constraints$pTSgroup) <- paste0(probs[-11]*100,'-',probs[-11]*100+10)
levels(constraints$pHIgroup) <- paste0(probs[-11]*100,'-',probs[-11]*100+10)
keep <- c('gene','pLI','loeuf','decile','pTS','pHI')

# Combine olink data
df <- olink$`Tissue-specific essential genes`
colnames(df)[1] <- c('gene')
mergedf <- merge(df, constraints[,keep, with = F], all.x = T, all.y=F, by = 'gene')
mergedf$numdecile <- as.numeric(as.character(mergedf$decile))
deciles <- as.numeric(levels(mergedf$decile))
deciles_new <- paste0(deciles-10, '-', deciles)
levels(mergedf$decile) <- deciles_new

# diagnostic plot
p1 <- ggplot(mergedf, aes(x=decile)) +
  geom_bar() + 
  xlab('LOEUF deciles') + 
  ylab('Frequency') +
  theme_bw() +
  ggtitle('LOUEF deciles vs Tissue-specific essential genes (OLINK)')

# deal with instances with only one tissue first
simple <- mergedf[!grepl('\\,',mergedf$`t statistic`),]
simple$p.val <- as.numeric(simple$p.val)
simple$p.adj <- as.numeric(simple$p.adj)
simple$logFC <- as.numeric(simple$logFC)
simple$label <- paste0(simple$gene, ' - ',simple$Tissues, ' (LOEUF = ',simple$loeuf,')')
simple$ADS <- as.numeric(simple$`Averege Dependency Score`)

p2 <- ggplot(simple, aes(x=logFC, y = -log10(p.adj), color = decile, size = 1/numdecile, label = gene)) +
  geom_point() +
  theme_bw() +
  labs(color = 'LOEUF decile', size = '1/LOUEF decile') + 
  geom_text_repel(data = simple[simple$numdecile <= 20 & simple$logFC < -0.5,], color = 'black', box.padding = 1, size = 3) +
  ggtitle('Tissue-specific essential genes (OLINK)',paste('Showing',nrow(simple),'rows of',nrow(mergedf)))

p3 <- ggplot(simple, aes(x=logFC, y = ADS, color = decile, size = 1/numdecile, label = gene)) +
  geom_point() +
  theme_bw() +
  ylab('Average Dependency Score') +
  labs(color = 'LOEUF decile', size = '1/LOUEF decile') + 
  geom_text_repel(data = simple[simple$numdecile <= 20 & simple$ADS < quantile(simple$ADS, probs = 0.05, na.rm = T),], color = 'black', box.padding = 1, size = 3, max.overlaps = Inf) +
  ggtitle('Tissue-specific essential genes (OLINK)',paste('Showing',nrow(simple),'rows of',nrow(mergedf)))

p4 <- ggplot(simple, aes(x=pTS, y = pHI, color = ADS, size = 1/numdecile, label = gene)) +
  geom_point() +
  theme_bw() +
  ggtitle('Tissue-specific essential genes (OLINK)',paste('Showing',nrow(simple),'rows of',nrow(mergedf)))


pdf('derived/plots/210603_todd_olink_genes.pdf', width = 8, height = 6)
print(p1)
print(p2)
print(p3)
print(p4)
graphics.off()

#  
olink[['LOUEF deciles']] <- mergedf
olink$Sheet <- NULL
writexl::write_xlsx(olink, 'derived/tables/210603_todd_olink_genes.xlsx')




