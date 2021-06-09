

expression <- fread('derived/tables/210609_prt_rna_numerical.txt', sep = '\t')
complexity <- fread('derived/tables/210609_MANE.v0.93.UTR_features.txt', sep = '\t')
dt <- merge(complexity, expression, by.x = 'ensgid',by.y = 'gene.id')
dt <- dt[dt$tissue == 'Brain_Cortex',]
dt$rna <- (dt$rna - mean(dt$rna, na.rm = T))/sd(dt$rna, na.rm = T)
dt$prt <- (dt$prt - mean(dt$prt, na.rm = T))/sd(dt$prt, na.rm = T)
dt$delta <- dt$rna - dt$prt

hist(dt$delta)
hist(dt$rna)

# initial peaking
graphics.off()

pdf('210609_numeric_concordance_vs_utrs.pdf', width = 12, height = 9)
par(mfrow=c(4,3))

# U5
plot(dt$rna, log(dt$u5_len), main = 'RNA')
plot(dt$prt, log(dt$u5_len), main = 'Protein')
plot(dt$delta, log(dt$u5_len), main = 'Delta')

plot(dt$rna, dt$u5_ORF, main = 'RNA')
plot(dt$prt, dt$u5_ORF, main = 'Protein')
plot(dt$delta, dt$u5_ORF, main = 'Delta')

plot(dt$rna, dt$u5_GC, main = 'RNA')
plot(dt$prt, dt$u5_GC, main = 'Protein')
plot(dt$delta, dt$u5_GC, main = 'Delta')

plot(dt$rna, dt$u5_oORF, main = 'RNA')
plot(dt$prt, dt$u5_oORF, main = 'Protein')
plot(dt$delta, dt$u5_oORF, main = 'Delta')

## U3
plot(dt$rna, log(dt$u3_len), main = 'RNA')
plot(dt$prt, log(dt$u3_len), main = 'Protein')
plot(dt$delta, log(dt$u3_len), main = 'Delta')

plot(dt$rna, dt$u3_ORF, main = 'RNA')
plot(dt$prt, dt$u3_ORF, main = 'Protein')
plot(dt$delta, dt$u3_ORF, main = 'Delta')

plot(dt$rna, dt$u3_GC, main = 'RNA')
plot(dt$prt, dt$u3_GC, main = 'Protein')
plot(dt$delta, dt$u3_GC, main = 'Delta')

plot(dt$rna, dt$u3_oORF, main = 'RNA')
plot(dt$prt, dt$u3_oORF, main = 'Protein')
plot(dt$delta, dt$u3_oORF, main = 'Delta')
graphics.off()


expression <- fread( 'derived/tables/210609_prt_rna_binary_concordance.txt', sep = '\t')
dt <- merge(complexity, expression, by.x = 'ensgid',by.y = 'ensembl_id')


p1 <- ggplot(dt[dt$u5_oORF != 0], aes(x=value,y=u5_oORF, fill = value)) +
  geom_jitter(size=0.5, alpha = 0.1) +
  geom_violin(alpha = 0.9) +
  facet_wrap(~tissue) +
  theme_bw()
p1


p2 <- ggplot(dt, aes(x=value,y=u5_ORF, fill = value)) +
  geom_jitter(size=0.5, alpha = 0.1) +
  geom_violin(alpha = 0.9) +
  facet_wrap(~tissue) +
  theme_bw()
p2

p3 <- ggplot(dt, aes(x=value,y=u5_ORF, fill = value)) +
  geom_jitter(size=0.5, alpha = 0.1) +
  geom_violin(alpha = 0.9) +
  facet_wrap(~tissue) +
  theme_bw()
p3



