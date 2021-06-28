# 21-06-28: are protein/rna ratios concordant across tissues?

# get protein / RNA 
expression <- fread('derived/tables/210609_prt_rna_numerical.txt', sep = '\t')
complexity <- fread('derived/tables/210609_MANE.v0.93.UTR_features.txt', sep = '\t')
dt <- merge(complexity, expression, by.x = 'ensgid',by.y = 'gene.id')
dt <- dt[dt$prt != 0, ]
dt$prt_rna_ratio <-  dt$prt-dt$rna
dt$rna_prt_ratio <-  dt$rna-dt$prt

# selected DDG2P genes that are dose sensitive
collins2020 <- setDT(read_xlsx('~/Projects/08_genesets/genesets/data/dosage/Collins2021medrxvic_supplementary_data.xlsx', 13))
constraints <- fread('~/Projects/08_genesets/genesets/data/gnomad/karczewski2020/supplementary_dataset_11_full_constraint_metrics.tsv') #transcript-level 
g2p <- fread('extdata/DDG2P_4_6_2021.csv')
g2p <- g2p[g2p$`DDD category` %in% c('confirmed','probable')] # 1809 genes
colnames(g2p) <- gsub(' ','_',colnames(g2p))
g2p_genes_dosage <- g2p$gene_symbol[g2p$gene_symbol %in% collins2020$gene[collins2020$pHI >= 0.84 | collins2020$pTS >= 0.993]]

# sample protein-rna ratios
sim_prt_rna <- function(gene_symbols, reps, n, verbose = T){
  sd_narm <- function(x) sd(x, na.rm = T)
  iter <- 1:reps
  emp <- lapply(iter, function(i){
    if ( (i %% 100 == 0) & verbose) print(paste0(round((i/reps)*100, 10),'%'))
    selected_genes <- sample(unique(gene_symbols), n)
    selected_rows <- dt[dt$gene_symbol %in% selected_genes]
    aggr <- aggregate(prt_rna_ratio ~ gene_symbol, data = selected_rows, FUN = sd_narm)
    return(mean(aggr$prt_rna_ratio, na.rm = T))
  })
  return(as.numeric(na.omit(unlist(emp))))
}


sim_samples <- function(x, k = 1000, f = function(x) sd(x, na.rm = T), replace = T, probs = c(0.025,0.975)){
  simsamples <- replicate(k, sample(dt$prt_rna_ratio, replace = replace))
  simmedians <- apply(simsamples, 2, f)
  quantile(simmedians, probs)
}


## are g2p dosage genes more variable in protein/rna ratio compared to other genes? 
mean_sd <- sim_samples(dt$prt_rna_ratio, 10, function(x) sd(x, na.rm = T))





# plot time series of selected genes
set.seed(101)
selected_genes <- sample(g2p_genes_dosage,5)
ggplot(dt[dt$gene_symbol %in% selected_genes], aes(x=tissue, y = prt_rna_ratio, color = gene_symbol, group = gene_symbol)) +
  geom_point() +
  geom_line() 

# 95% bootstrap confidence interval for SD



# plot emprical mean of standard deviations of all genes
emp1.nums <- sd(sim_prt_rna(dt$gene_symbol, 200, 10))
plot(density(emp1.nums), main = 'Mean of Standard Deviations')

# compare with mean sd of 100 G2P genes
emp.g2p <- sim_prt_rna(g2p$gene_symbol, 1000, 100)
lines(density(emp.g2p), lty = 2, col = 'red')
t.test(emp1.nums, emp.g2p)

# compare with mean of dosage sensitive g2p genes
emp.dosage <- sim_prt_rna(g2p_genes_dosage, 1000, 100)
lines(density(emp.dosage), lty = 2, col = 'blue')
t.test(emp1.nums, emp.dosage)

# compare with mean of dosage sensitive g2p genes
probs <- complexity$u5_ORF/(max(complexity$u5_ORF)+1)
genes_orf <- sample(complexity$gene_symbol, 2000, prob = probs)


emp.orf <- sim_prt_rna(complexity$gene_symbol[complexity$u5_ORF > 5], 1000, 100)
lines(density(emp.orf), lty = 2, col = 'orange')
t.test(emp1.nums, emp.orf)

plot(density(emp1.nums), main = 'Mean of Standard Deviations')
for (i in 1:3){
  emp.orf <- sim_prt_rna(complexity$gene_symbol[complexity$u5_max_kozak == i], 500, 100)
  #emp.orf <- sim_prt_rna(complexity$gene_symbol[complexity$u5_ORF > i], 100, 100)
  lines(density(emp.orf), lty = 2, col = 'blue')
  #t.test(emp1.nums, emp.orf)
  
}



# are genes with uORFs more variable in protein/RNA compared to other genes?













