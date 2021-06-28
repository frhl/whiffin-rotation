

expression <- fread('derived/tables/210609_prt_rna_numerical.txt', sep = '\t')
complexity <- fread('derived/tables/210609_MANE.v0.93.UTR_features.txt', sep = '\t')
dt <- merge(complexity, expression, by.x = 'ensgid',by.y = 'gene.id')
