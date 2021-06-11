# author: frederik heymann (21-06-09)
# description: compare complexity of UTR across different tissue/protein expression tables

devtools::load_all()

# load expression and UTR data
expression <- fread('derived/tables/210609_prt_rna_numerical.txt', sep = '\t')
complexity <- fread('derived/tables/210611_MANE.v0.93.UTR_features.txt', sep = '\t')
dt <- merge(complexity, expression, by.x = 'ensgid',by.y = 'gene.id')
dt$rna_std <- (dt$rna - mean(dt$rna, na.rm = T))/sd(dt$rna, na.rm = T)
dt$prt_std <- (dt$prt - mean(dt$prt, na.rm = T))/sd(dt$prt, na.rm = T)
dt <- dt[dt$prt != 0, ]


# setup ggplot
pdf('derived/210611_UTR_complexity_vs_continious.pdf', width = 12, height = 12)

ggplot(dt[dt$u5_ORF > 1], aes(x=rna_std, y = prt_std, color = u5_ORF)) +
  geom_point() +
  facet_wrap(~tissue) +
  xlab('Normalized Log RNA expression') +
  ylab('Normalized Log Protein Expression') +
  ggtitle('U5_ORF (u5_ORF > 1)','2nd order polynominal fit') +
  geom_hline(yintercept = 0, alpha = 0.2) +
  geom_vline(xintercept = 0, alpha = 0.2) +
  geom_smooth(color = 'red', linetype = 'dashed', formula = y ~ poly(x, 2)) +
  theme_bw() 

text = paste0(sum(complexity$u5_oORF_inframe > 0),'/',nrow(complexity),' genes')
ggplot(dt[dt$u5_oORF_inframe > 0], aes(x=rna_std, y = prt_std, color = u5_oORF_inframe)) +
  geom_point() +
  facet_wrap(~tissue) +
  xlab('Normalized Log RNA expression') +
  ylab('Normalized Log Protein Expression') +
  ggtitle('U5_oORF_cds_inframe (uORF_inframe > 0)', text) +
  geom_hline(yintercept = 0, alpha = 0.2) +
  geom_vline(xintercept = 0, alpha = 0.2) +
  ylim(-5,5) +
  geom_smooth(color = 'red', linetype = 'dashed', formula = y ~ x) +
  theme_bw() 

text = paste0(sum(complexity$u5_oORF_outframe > 0),'/',nrow(complexity),' genes')
ggplot(dt[dt$u5_oORF_outframe > 0], aes(x=rna_std, y = prt_std, color = u5_oORF_outframe)) +
  geom_point() +
  facet_wrap(~tissue) +
  xlab('Normalized Log RNA expression') +
  ylab('Normalized Log Protein Expression') +
  ggtitle('U5_oORF_cds_outframe (uORF_outframe > 0)', text) +
  geom_hline(yintercept = 0, alpha = 0.2) +
  geom_vline(xintercept = 0, alpha = 0.2) +
  geom_smooth(color = 'red', linetype = 'dashed', formula = y ~ x) +
  theme_bw()

text = paste0(sum(complexity$u5_oORF_altered_cds > 0),'/',nrow(complexity),' genes')
ggplot(dt[dt$u5_oORF_altered_cds > 0], aes(x=rna_std, y = prt_std, color = u5_oORF_outframe)) +
  geom_point() +
  facet_wrap(~tissue) +
  xlab('Normalized Log RNA expression') +
  ylab('Normalized Log Protein Expression') +
  ggtitle('u5_oORF_altered_cds i.e. truncations/elongation of CDS (u5_oORF_altered_cds > 0)', text) +
  geom_hline(yintercept = 0, alpha = 0.2) +
  geom_vline(xintercept = 0, alpha = 0.2) +
  geom_smooth(color = 'red', linetype = 'dashed', formula = y ~ x) +
  theme_bw()

text = paste0(sum(complexity$u5_oORF_kozak > 0),'/',nrow(complexity),' genes')
ggplot(dt[dt$u5_oORF_kozak != 0], aes(x=rna_std, y = prt_std, color = u5_oORF_kozak)) +
  geom_point(data = dt[dt$u5_oORF_kozak %in% c(1)], alpha = 0.9) +
  geom_point(data = dt[dt$u5_oORF_kozak %in% c(2)], alpha = 0.9) +
  geom_point(data = dt[dt$u5_oORF_kozak %in% c(3)], alpha = 0.9) +
  facet_wrap(~tissue) +
  xlab('Normalized Log RNA expression') +
  ylab('Normalized Log Protein Expression') +
  ggtitle('u5_oORF_kozak i.e. kozak strength truncations/elongation of CDS (u5_oORF_kozak > 0) (3=strong, 2=moderate, 1=weak)', text) +
  geom_hline(yintercept = 0, alpha = 0.2) +
  geom_vline(xintercept = 0, alpha = 0.2) +
  geom_smooth(color = 'red', linetype = 'dashed', formula = y ~ x) +
  theme_bw()

text = paste0(sum(complexity$u5_AUG > 0),'/',nrow(complexity),' genes')
ggplot(dt[u5_AUG > 0], aes(x=rna_std, y = prt_std, color = log2(u5_AUG))) +
  geom_point() +
  facet_wrap(~tissue) +
  xlab('Normalized Log RNA expression') +
  ylab('Normalized Log Protein Expression') +
  geom_hline(yintercept = 0, alpha = 0.2) +
  geom_vline(xintercept = 0, alpha = 0.2) +
  ggtitle('U5_AUG (U5_AUG > 0)', text) +
  theme_bw() 

graphics.off()


# compare top 10 % versus bottom 10 %
tissue <- unique(dt$tissue)
hows <- c('u5_AUG','u5_ORF','u5_oORF_outframe','u5_oORF_inframe')
how = 'u5_oORF_inframe'
perc <- c(0.2, 0.8)

lst <- list()

for (how in hows){ #hows, function(how){
  
  comp.rna <- do.call(rbind, lapply(tissue, function(x){
    
    df <- dt[dt$tissue %in% x]
    df$delta <- df$rna_std - df$prt_std
    quant <- quantile(df$delta, probs = perc, na.rm = T)
    lower <- quant[1]
    upper <- quant[2]
    lower_count <- df[df$delta < lower, ][[how]]
    upper_count <- df[df$delta > upper, ][[how]]
    
    print(paste0(length(lower_count) - length(upper_count)))
    
    d1 <- data.frame(tissue = x, percentile = perc[1], how = how, table(lower_count))
    d2 <- data.frame(tissue = x, percentile = perc[2], how = how, table(upper_count))
    
    colnames(d1)[4] <- 'count'
    colnames(d2)[4] <- 'count'
    return(rbind(d1, d2))
    
  }))
  comp.rna$perc <- factor(comp.rna$perc)
  
  
  lst[[how]] <- ggplot(comp.rna[comp.rna$count != 0,], aes(x = count, y = Freq, group = percentile, fill = perc)) +
    geom_bar(stat="identity") +
    ggtitle(how) +
    theme_bw()
    #facet_wrap(~tissue)
  
  
}




hist(dt$delta)
hist(dt$rna)

# initial peaking
graphics.off()


# understand the relationshop between concordance and not
expression1 <- fread( 'derived/tables/210609_prt_rna_binary_concordance.txt', sep = '\t')
expression1 <- expression1[!duplicated(expression1),]
expression2 <- fread('derived/tables/210609_prt_rna_numerical.txt', sep = '\t')
expression2 <- expression2[!duplicated(expression2),]
colnames(expression2)[1] <- 'ensembl_id'
dt <- merge(expression1, expression2)
dt$rna_std <- (dt$rna - mean(dt$rna, na.rm = T))/sd(dt$rna, na.rm = T)
dt$prt_std <- (dt$prt - mean(dt$prt, na.rm = T))/sd(dt$prt, na.rm = T)

ggplot(dt, aes(x=rna_std, y = prt_std, color = value)) +
  geom_point() +
  facet_wrap(~tissue) +
  xlab('Normalized Log RNA expression') +
  ylab('Normalized Log Protein Expression') +
  geom_abline(linetype = 'dashed') +
  theme_bw() 

# Groups and values
expression <- fread( 'derived/tables/210609_prt_rna_binary_concordance.txt', sep = '\t')
expression <- expression[!duplicated(expression),]
dt <- merge(complexity, expression, by.x = 'ensgid',by.y = 'ensembl_id')



p1 <- ggplot(dt, aes(x=value,y=u5_oORF_inframe, fill = value)) +
  geom_jitter(size=0.5, alpha = 0.1) +
  geom_violin(alpha = 0.9) +
  facet_wrap(~tissue) +
  theme_bw() +
  xlab('Concordance Group') +
  ggtitle('oORF_inframe')
p1

p1 <- ggplot(dt, aes(x=value,y=u5_oORF_outframe, fill = value)) +
  geom_jitter(size=0.5, alpha = 0.1) +
  geom_violin(alpha = 0.9) +
  facet_wrap(~tissue) +
  theme_bw() +
  xlab('Concordance Group') +
  ggtitle('oORF_outframe')
p1

p1 <- ggplot(dt, aes(x=value,y=u5_AUG, fill = value)) +
  geom_jitter(size=0.5, alpha = 0.1) +
  geom_violin(alpha = 0.9) +
  facet_wrap(~tissue) +
  theme_bw() +
  xlab('Concordance Group') +
  ggtitle('uAUG')
p1







