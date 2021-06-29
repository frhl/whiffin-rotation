# 21-06-28: are protein/rna ratios concordant across tissues?

devtools::load_all()
library(readxl)

# get protein / RNA 
expression <- fread('derived/tables/210609_prt_rna_numerical.txt', sep = '\t')
complexity <- fread('derived/tables/210615_MANE.v0.93.UTR_features.txt', sep = '\t')
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

# simulate samples
sim_samples <- function(x, k = 1000, func = function(x) sd(x, na.rm = T), probs = c(0.025,0.5,0.975)){
  simsamples <- replicate(k, sample(x, replace = T))
  simmedians <- apply(simsamples, 2, func)
  outdf <- t(as.data.frame(quantile(simmedians, probs = probs)))
  rownames(outdf) <- NULL
  return(outdf)
}

# is there a discrepancy of gene/rna ratio across tissues?
run_mean_conf <- function(dt, func, k = 100){
  tissues <- unique(dt$tissue)
  result <- lapply(tissues, function(cur_tissue) sim_samples(dt$prt_rna_ratio[dt$tissue %in% cur_tissue], k = k, func = func))
  mat <- data.frame(do.call(rbind, result))
  colnames(mat) <- c('lower','median','upper')
  mat$tissue <- tissues
  return(mat)
}

sets <- list(
  all = dt,
  g2p = dt[dt$gene_symbol %in% g2p_genes_dosage],
  u5_no_orf = dt[dt$u5_ORF == 0],
  u5_orf = dt[dt$u5_ORF > 0],
  u5_oorf = dt[dt$u5_oORF_altered_cds > 0]
)

g <- function(x) mean(x, na.rm = T)
f <- function(x) sd(x, na.rm = T)


# why is there a higher ratio of protein to mRNA for u5_orf
# means
res_mean <- lapply(names(sets), function(s) data.frame(run_mean_conf(sets[[s]], g, k = 20), how = s))
res_mean1 <- do.call(rbind, res_mean)

ggplot(res_mean1, aes(x=median, xmin=lower, xmax=upper, y = reorder(tissue, median), color = how, group = how)) +
  geom_errorbar() +
  geom_point() +
  xlab('Mean (Protein/RNA)') +
  ylab('Tissues')

# standard deviations
res_sds <- lapply(names(sets), function(s) data.frame(run_mean_conf(sets[[s]], f, k = 100), how = s))
res_sds1 <- do.call(rbind, res_sds)
res_sds1
sd(dt$prt_rna_ratio[dt$tissue %in% 'Brain_Cortex'], na.rm = T)

ggplot(res_sds1, aes(x=median, xmin=lower, xmax=upper, y = reorder(tissue, median), color = how, group = how)) +
  geom_errorbar(position = 'dodge') +
  #geom_point(position = 'dodge') +
  xlab('Standard Deviation (Protein/RNA)') +
  ylab('Tissues')








## are g2p dosage genes more variable in protein/rna ratio compared to other genes? 


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













