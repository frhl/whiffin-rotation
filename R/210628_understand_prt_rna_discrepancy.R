# 21-06-28: why does it seem like there is more protein per rna?

devtools::load_all()

# get protein / RNA 
expression <- fread('derived/tables/210609_prt_rna_numerical.txt', sep = '\t')
complexity <- fread('derived/tables/210615_MANE.v0.93.UTR_features.txt', sep = '\t')
dt <- merge(complexity, expression, by.x = 'ensgid',by.y = 'gene.id')
dt <- dt[dt$prt != 0, ]
dt$prt_rna_ratio <-  dt$prt-dt$rna
dt$rna_prt_ratio <-  dt$rna-dt$prt