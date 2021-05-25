# author: Frederik Heymann Lassen, 25-May-2021
# description: investigate upper-bound LoF fraction of O/E protein trunctaing variation
# and haploinsifficient/triplosensitive genes

library(readxl)
library(data.table)
library(ggplot2)
library(ggrepel)
library(cowplot)


# load data
collins2020 <- setDT(read_xlsx('~/Projects/08_genesets/genesets/data/dosage/Collins2021medrxvic_supplementary_data.xlsx', 13))
constraints <- fread('~/Projects/08_genesets/genesets/data/gnomad/karczewski2020/supplementary_dataset_11_full_constraint_metrics.tsv') #transcript-level 

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

# set gene of interest
constraints$goi_hi <- constraints$pHIgroup == '90-100' & constraints$loeuf > deciles[8]
constraints$goi_ts <- constraints$pTSgroup == '90-100' & constraints$loeuf > deciles[8]
constraints$type <- ifelse(constraints$goi_hi, 'HI', ifelse(constraints$goi_ts,'TS',NA))

# condense data so it can be plotted on top
plot_df <- constraints[,c('gene','decile','type','loeuf','pHI','pTS')]
plot_df <- melt(plot_df, id.vars = 1:4, measure.vars = 5:6)
plot_df <- plot_df[!is.na(plot_df$loeuf),]
plot_df <- plot_df[!is.na(plot_df$value),]

# plot it
pdf('derived/210525_karczewski_vs_collins.pdf',width = 12, height = 5)
ggplot(plot_df, aes(x=decile, y = value, label = gene)) +
  geom_jitter() +
  xlab('LOEUF decile (%) [Karczewski 2020]') +
  ylab('pHI / pTS value [Collins 2021]') +
  facet_grid(~variable) +
  #geom_text_repel(data = plot_df[!is.na(plot_df$type),], show.legend = ) +
  geom_hline(yintercept = 0.75, linetype = 'dashed', color = 'grey') +
  theme_bw()

# plot counts 
mat_pt <- as.data.frame(table(constraints$decile, constraints$pTSgroup))
mat_hi <- as.data.frame(table(constraints$decile, constraints$pHIgroup))
p1 <- ggplot(mat_pt, aes(x=Var1, y=Var2, fill=Freq, label = Freq)) +
  geom_tile(show.legend = F) + geom_text() + theme_minimal() + 
  xlab('LOEUF decile (%) [Karczewski 2020]') + 
  ylab('pTS bins [Collins 2021]')
p2 <- ggplot(mat_hi, aes(x=Var1, y=Var2, fill=Freq, label = Freq)) +
  geom_tile(show.legend = F) + geom_text() + theme_minimal() + 
  xlab('LOEUF decile (%) [Karczewski 2020]') + 
  ylab('pHI bins [Collins 2021]')

plot_grid(p1, p2, labels = 'AUTO')

# plot scatter plot relation between the two
p3 <- ggplot(constraints, aes(x=loeuf, y=pTS, label = gene, color = type)) +
  geom_point(show.legend = F) +
  geom_vline(xintercept = deciles, linetype = 'dashed', color = 'red') + 
  geom_hline(yintercept = c(0.5, 0.75, 0.9), linetype = 'dashed', color = 'red') +
  geom_text_repel(data=constraints[!is.na(constraints$type),],show.legend = F) +
  xlab('LOEUF decile (Karczewski 2020)') +
  ylab('pTS (Collins 2021)') +
  theme_bw()
p4 <- ggplot(constraints, aes(x=loeuf, y=pHI, label = gene, color = type)) +
  geom_point(show.legend = F) +
  geom_vline(xintercept = deciles, linetype = 'dashed', color = 'red') + 
  geom_hline(yintercept = c(0.5, 0.75, 0.9), linetype = 'dashed', color = 'red') +
  geom_text_repel(data=constraints[!is.na(constraints$type),],show.legend = F) +
  xlab('LOEUF (Karczewski 2020)') +
  ylab('pHI (Collins 2021)') +
  theme_bw()

plot_grid(p3, p4, labels = 'AUTO')
graphics.off()

## save data
out <- constraints
out <- out[!(is.na(out$pHI) | is.na(out$pTS) | is.na(out$loeuf)),]
out <- out[with(out, order(-decile, pTSgroup, pHIgroup)),]
out <- out[out$pHI > 0.5 | out$pTS > 0.5,]
out <- out[,c('gene','transcript','canonical','loeuf','decile','pHI','pHIgroup','pTS','pTSgroup')]
writexl::write_xlsx(out,path = 'derived/210525_karczewski_vs_collins_0.5cutoff.xlsx')


