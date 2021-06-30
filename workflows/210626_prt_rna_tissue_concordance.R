# 21-06-28: are protein/rna ratios concordant across tissues?

devtools::load_all()
library(readxl)

# get protein / RNA 
expression <- fread('derived/tables/210609_prt_rna_numerical.txt', sep = '\t')
complexity <- fread('derived/tables/210629_MANE.v0.95.UTR_features.txt', sep = '\t')
dt <- merge(complexity, expression, by.x = 'ensgid',by.y = 'gene.id')
dt <- dt[dt$prt != 0, ]
dt$prt_rna_ratio <-  dt$prt-dt$rna
dt$rna_prt_ratio <-  dt$rna-dt$prt

# selected DDG2P genes that are dose sensitive
collins2020 <- setDT(read_xlsx('~/Projects/08_genesets/genesets/data/dosage/Collins2021medrxvic_supplementary_data.xlsx', 13))
constraints <- fread('~/Projects/08_genesets/genesets/data/gnomad/karczewski2020/supplementary_dataset_11_full_constraint_metrics.tsv') #transcript-level 
g2p <- fread('extdata/DDG2P_4_6_2021.csv')
g2p <- g2p[g2p$`DDD category` %in% c('confirmed','probable')] # 1809 genes
g2p <- g2p[g2p$`mutation consequence` %in% 'loss of function']  # 1431 genes
colnames(g2p) <- gsub(' ','_',colnames(g2p))
g2p_genes_dosage <- g2p$gene_symbol[g2p$gene_symbol %in% collins2020$gene[collins2020$pHI >= 0.84 | collins2020$pTS >= 0.993]]
#g2p_genes_dosage <- g2p$gene_symbol

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
run_mean_conf <- function(dt, func, k = 100, param = 'prt', tissues = unique(dt$tissue)){
  result <- lapply(tissues, function(cur_tissue) sim_samples(dt[[param]][dt$tissue %in% cur_tissue], k = k, func = func))
  mat <- data.frame(do.call(rbind, result))
  colnames(mat) <- c('lower','median','upper')
  mat$tissue <- tissues
  return(mat)
}

# stability -> more protein? 
length(unique(dt$gene_symbol[dt$u5_ORF == 0]))
len_no_orf_median <- median(dt[dt$u5_ORF == 0]$u5_len)
len_orf_median <- median(dt[dt$u5_ORF != 0]$u5_len)


sets <- list(
  all = dt,
  #secreted = dt[dt$ensgid %in% x],
  #g2p = dt[dt$gene_symbol %in% g2p_genes_dosage],
  #g2p_uorf = dt[dt$gene_symbol %in% g2p_genes_dosage & dt$u5_ORF > 0],
  #g2p_no_uorf = dt[dt$gene_symbol %in% g2p_genes_dosage &  dt$u5_ORF == 0],
  
  
  #u5_orf_1_50 = dt[dt$u5_len < 50],
  #u5_orf_50_100 = dt[dt$u5_len >= 50 & dt$u5_len < 100],
  #u5_orf_100_200 = dt[dt$u5_len >= 100 & dt$u5_len < 200],
  #u5_orf_200_300 = dt[dt$u5_len >= 100 & dt$u5_len < 200],
  #u5_orf_300_upper = dt[dt$u5_len >= 300]

  #u5_no_orf_secreted = dt[dt$u5_ORF == 0 & dt$ensgid %in% x],
  #u5_orf_secreted = dt[dt$u5_ORF > 0 & dt$ensgid %in% x],
  #u5_oorf_secreted = dt[dt$u5_oORF_altered_cds > 0 & dt$ensgid %in% x],
    
  u5_no_orf = dt[dt$u5_ORF == 0],
  u5_orf = dt[dt$u5_ORF > 0 ],
  u5_oorf = dt[dt$u5_oORF_altered_cds > 0]
)


#sets <- list(
  #all = dt,
  #g2p = dt[dt$gene_symbol %in% g2p_genes_dosage],
  #g2p_uorf = dt[dt$gene_symbol %in% g2p_genes_dosage & dt$u5_ORF > 0],
  #g2p_no_uorf = dt[dt$gene_symbol %in% g2p_genes_dosage &  dt$u5_ORF == 0],
  #u5_no_orf = dt[dt$u5_ORF == 0],
  #u5_orf = dt[dt$u5_ORF > 0]
  #u5_oorf = dt[dt$u5_oORF_altered_cds > 0]
#)

g <- function(x) mean(x, na.rm = T)
f <- function(x) sd(x, na.rm = T)









summary(lm(u5_ORF ~ u5_GC, data = dt))
summary(lm(prt ~ u5_ORF + u5_len + u5_GC, data = dt))
summary(lm(rna ~ u5_ORF + u5_len + u5_GC, data = dt))



fit <- lm(u5_len ~ prt + tissue, data = dt)
summary(fit)


#tissue <- unique(dt$tissue)
#dq <- dq[!is.na(dq$prt) & !is.na(dq$rna)]
#dq <- dq[!duplicated(dq[,c('rna','prt')])]
#lapply(tissue, function(t) cor(dq$prt[dq$tissue %in% t], dq$rna[dq$tissue %in% t])  )


tissues <- unique(dt$tissue)

# 1) comparing mean expression across all tissues
rna_mean <- lapply(names(sets), function(s) data.frame(run_mean_conf(sets[[s]], g, k = 20, 'rna', ), how = s))
rna_mean1 <- do.call(rbind, rna_mean)

p1 <- ggplot(rna_mean1, aes(x=median, xmin=lower, xmax=upper, y = reorder(how, median), fill = how, group = how)) +
  geom_boxplot(color = 'black', show.legend = F) +
  xlab('Mean RNA (Log2 TPM)') +
  ylab('UTR Regulatory Element') +
  theme_bw()

prt_mean <- lapply(names(sets), function(s) data.frame(run_mean_conf(sets[[s]], g, k = 20, 'prt', ), how = s))
prt_mean1 <- do.call(rbind, prt_mean)
p2 <- ggplot(prt_mean1, aes(x=median, xmin=lower, xmax=upper, y = reorder(how, median), fill = how, group = how)) +
  geom_boxplot(color = 'black') +
  xlab('Mean Protein (Log2 abundance)') +
  ylab('') +
  theme_bw()

plot_grid(p1, p2, rel_widths = c(0.41, 0.59))

# 2) comparing standard deviation for expression across all tissues
rna_sd <- lapply(names(sets), function(s) data.frame(run_mean_conf(sets[[s]], f, k = 20, 'rna', ), how = s))
rna_sd1 <- do.call(rbind, rna_sd)

p3 <- ggplot(rna_sd1, aes(x=median, xmin=lower, xmax=upper, y = reorder(how, median), fill = how, group = how)) +
  geom_boxplot(color = 'black', show.legend = F) +
  xlab('Standard Deviation RNA (Log2 TPM)') +
  ylab('UTR Regulatory Element') +
  theme_bw()

prt_sd <- lapply(names(sets), function(s) data.frame(run_mean_conf(sets[[s]], f, k = 20, 'prt', ), how = s))
prt_sd1 <- do.call(rbind, prt_sd)
p4 <- ggplot(prt_sd1, aes(x=median, xmin=lower, xmax=upper, y = reorder(how, median), fill = how, group = how)) +
  geom_boxplot(color = 'black') +
  xlab('Standard Deviation Protein (Log2 abundance)') +
  ylab('') +
  theme_bw()

plot_grid(p3, p4, rel_widths = c(0.41, 0.59))

#

# Get standard deviations across tissue
aggr <- aggregate(rna ~ enstid_version, data = dt, FUN = function(x) sd(x, na.rm = T))
colnames(aggr)[2] <- "rna.tissue.sd"
mrg <- merge(aggr, complexity)
aggr <- aggregate(prt ~ enstid_version, data = dt, FUN = function(x) sd(x, na.rm = T))
colnames(aggr)[2] <- "prt.tissue.sd"
mrg <- merge(aggr, mrg)



# check if U5 ORFs correlate with standard deviation?
summary(lm(prt.tissue.sd ~ u5_ORF + u5_len, data = mrg))
p2a <- ggplot(mrg, aes(x=prt.tissue.sd, y = u5_ORF)) +
  geom_point(alpha = 0.9, color = 'darkblue') +
  xlab('Protein Abundance standard deviation (Across 32 tissues)') +
  ylab('U5 ORFs') +
  geom_vline(xintercept = median(mrg$prt.tissue.sd, na.rm = T), linetype = 'dashed') +
  theme_bw() + 
  ggtitle('U5 ORFs as a function of Standard Deviation across tissues',
          'LM p-value "prt.tissue.sd ~ u5_ORF + u5_len", u5_ORF: 0.498 ')


summary(lm(rna.tissue.sd ~ u5_ORF + u5_len, data = mrg))
p2b <- ggplot(mrg, aes(x=rna.tissue.sd, y = u5_ORF)) +
  geom_point(alpha = 0.9, color = 'red') +
  xlab('RNA Abundance standard deviation (Across 32 tissues)') +
  ylab('U5 ORFs') +
  geom_vline(xintercept = median(mrg$prt.tissue.sd, na.rm = T), linetype = 'dashed') +
  theme_bw() + 
  ggtitle('U5 ORFs as a function of Standard Deviation across tissues',
          'LM p-value "rna.tissue.sd ~ u5_ORF + u5_len", u5_ORF: 0.597 ')

nrow(mrg)
plot_grid(p2b, p2a)










# Having uORF versus not having uORF

dt$u5_ORF_present <- dt$u5_ORF > 0
dt$u5_aug_present <- dt$u5_AUG > 0

fit <- glm(u5_ORF_present ~ u5_aug_present, family = binomial, data = dt)

summary(fit)


summary(lm(rna ~ u5_ORF + log(u5_len) + u5_GC, data = dt))





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













