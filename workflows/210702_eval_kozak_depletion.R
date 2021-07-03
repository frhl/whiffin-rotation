# Evaluate kozak depletion

library(data.table)
d_obs_ci <- fread('derived/210702_MANE.GRCh38.v0.95_kozak_codons_expt_ci.csv', sep = ',')
d_obs <- fread('derived/210702_MANE.GRCh38.v0.95_kozak_codons_expt.csv', sep = ',')
d_expt <- fread('derived/210702_MANE.GRCh38.v0.95_kozak_codons_obs.csv', sep = ',')


expression <- fread('derived/tables/210609_prt_rna_numerical.txt', sep = '\t')
#aggr_rna <- aggregate(rna ~ gene.id, data = expression, FUN = function(x) quantilef(x))
aggr_rna <- aggregate(rna ~ gene.id, data = expression, FUN = function(x) max(x, na.rm = T))
colnames(aggr_rna) <- c("gene.id","rna")
seq_quantile <- seq(0,1, by = 0.1)
quantile_rna <- quantile(aggr_rna$rna, probs = seq_quantile)
aggr_rna$percentile <- cut(aggr_rna$rna, quantile_rna)
levels(aggr_rna$percentile) <- seq_quantile*100

mrg <- merge(d_expt, d_obs)
aggr_mrg <- merge(mrg, aggr_rna, by.x = 'ensgid', by.y = 'gene.id')
aggr_mrg$ensgid_version <- NULL

obs_aggr <- aggr_mrg[,get('(enstid)|(obs)|(perc)',aggr_mrg), with = F]
expt_aggr <- aggr_mrg[,get('(enstid)|(expt)|(perc)',aggr_mrg), with = F]

expr1 <- do.call(rbind, lapply(seq_quantile*100, function(perc){
  
  selected_obs <- obs_aggr[obs_aggr$percentile == perc, get('obs',obs_aggr), with = F]
  selected_expt <- expt_aggr[expt_aggr$percentile == perc, get('expt',expt_aggr), with = F]
  d_out <- data.frame(perc = perc, obs = colSums(selected_obs), expt = colSums(selected_expt))  
  d_out$codon <- indexsplit(rownames(d_out), 2)
  rownames(d_out) <- NULL
  d_out$oe <- d_out$obs / d_out$expt
  return(d_out)
  
}))


expr1$codon <- factor(expr1$codon, levels = rev(unique(expr1$codon)))
ggplot(expr1, aes(x=codon, y = oe, group = perc, color = perc)) +
  geom_point(size = 2) +
  geom_line(size = 1) +
  labs(color = 'Percentile (Max RNA expression across 32 tissues)') +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  ggtitle('Observed versus Expected Kozak Sequence as a function of RNA expression',
          '1000 simulations preservering di-nt frequency for each sequnce') +
  ylab('Observed / Expected') +
  xlab('Kozak Strength') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = 'bottom')
