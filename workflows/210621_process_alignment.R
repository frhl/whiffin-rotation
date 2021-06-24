

# read alignment
d <- fread(file = 'derived/tables/210623_1000_G2P-confirmedprobable_five_prime_utr_alignment_NW.txt', sep = '\t')
d$V1 <- NULL

# remove duplicates genes
genes <- d$gene_symbol

#ensgids_version <- d$ensgid
#ensgids <- unlist(lapply(ensgids_version, function(x) strsplit(x, split = '\\.')[[1]][1]))
indicies <- (1:1000)[duplicated(d$gene_symbol)]
d$ensgid <- NULL

#d$gene_symbol <- NULL
d <- d[-indicies, ]
d <- d[,-(indicies+1), with = F]
rownames(d) <- genes[-indicies]
colnames(d) <- c('gene_symbol', genes[-indicies])
dim(d)
head(d)

# combine data
d <- as.data.table(melt(d, id.vars = 'gene_symbol'))
d <- d[d$value != '',]
colnames(d) <- c('g1','g2','full')
d$g2 <- as.character(d$g2)

# extract variables
d_extracted <- do.call(rbind, strsplit(d$full, split = '\\:'))
d_extracted <- as.data.table(apply(d_extracted, 2, as.numeric))
colnames(d_extracted) <- c('l','m','g')

# combine and use
final <- cbind(d, d_extracted)
final$full <- NULL

# Scoring
score_match <- 5
score_gap <- 2
final$cur_score <- final$m*score_match - (final$g*score_gap)^1.4
final$max_score <- final$l*score_match
final$score <- final$cur_score / final$max_score
#hist(final$score)

# select top 20%
top <- quantile(final$score, probs = 0.99)
final[final$score == max(final$score)]
final[final$score > top]
plot(final[final$score > top]$score, final[final$score > top]$l)

# select the pairs with high individual scores
count <- 0
genes <- unique(final$g1)
pairs <- lapply(genes, function(gt){
  count <<- count + 1
  pair <- final[(final$g1 == gt | final$g2 == gt) & score > top]
  pair$grp <- count
  return(pair)
})
names(pairs) <- genes

# check the group with highest median resemblancel
scores_list <- lapply(pairs, function(d) mean(d$score))
scores <- unlist(scores_list)
top_expression <- head(rev(sort(na.omit(scores))), n = 500)
top_pairs <- names(top_expression)
top_pairs_groups <- lapply(top_pairs, function(gt) unique(c(pairs[[gt]]$g1, pairs[[gt]]$g2)))


# load expression data: is protein expression similar for genes with similar sequence?
expression <- fread('derived/tables/210609_prt_rna_numerical.txt', sep = '\t')
count <- length(top_pairs)
prt <- lapply(top_pairs, function(gt){

  count <<- count -1
  current_genes <- unique(c(pairs[[gt]]$g1, pairs[[gt]]$g2))
  tmp <- expression[expression$gene.id %in% current_genes,]
  tmp$group <- count
  tmp$mean_score <- mean(pairs[[gt]]$score)
  return(tmp)
  
})

# go through all tissue
mat_prt <- do.call(rbind, prt)
tissues <- unique(mat_prt$tissue)

lm.1 <- lapply(tissues, function(ts){
  mat_prt_selected <- mat_prt[mat_prt$tissue == ts]
  fit <- lm(prt ~ mean_score, data =  mat_prt_selected)
  sumfit <- summary(fit)
  out <- as.data.frame(t(as.matrix(sumfit$coefficients[2,])))
  out$tissue <- ts
  return(out)
})
lm.1.res <- do.call(rbind, lm.1)
lm.1.res$FDR <- stats::p.adjust(lm.1.res$`Pr(>|t|)`)
lm.1.res[order(lm.1.res$FDR),]


pdf('derived/plots/210624_500_similarity_score_vs_rna_prt.pdf',width = 8, height = 6)
tissues_selected <- c("Brain_Cerebellum", "Brain_Cortex", 'Liver', 'Nerve_Tibial', 'Heart_Ventricle', 'Muscle_Skeletal')
#tissues_selected <- c("Brain_Cerebellum", "Brain_Cortex", 'Liver','Nerve_Tibial')
ggplot(mat_prt[mat_prt$tissue %in% tissues_selected & !is.na(mat_prt$prt)], aes(x=mean_score, y = prt, group = group)) +
  geom_point(size = 1, alpha = 0.8) +
  geom_smooth(method = "lm", col = "red") +
  facet_wrap(~tissue) +
  xlab('Similarity Score [0-1]') +
  ylab('Protein Expression') + 
  theme_grey()

sd1 <- function(x) sd(x, na.rm = T)
# standard deviation as a function of similarity
df_sd <- do.call(rbind, lapply(tissues, function(ts){
  
  # combine protein / rna data
  mat_prt_selected <- mat_prt[mat_prt$tissue == ts]
  aggr.prt <- aggregate(prt ~ mean_score, data = mat_prt_selected, sd1)
  colnames(aggr.prt)[2] <- 'prt'
  aggr.rna <- aggregate(rna ~ mean_score, data = mat_prt_selected, sd1)
  colnames(aggr.rna)[2] <- 'rna'
  aggr.lens <- aggregate(prt ~ mean_score, data = mat_prt_selected, length)
  colnames(aggr.lens)[2] <- 'obs'
  aggr <- merge(aggr.prt,aggr.rna)
  aggr <- merge(aggr, aggr.lens)
  aggr$tissue <- ts
  
  
  return(aggr)
}))

lm.2 <- lapply(tissues, function(ts){
  df_sd_selected <- df_sd[df_sd$tissue == ts,]
  fit <- lm(rna ~ mean_score, data =  df_sd_selected)
  sumfit <- summary(fit)
  out <- as.data.frame(t(as.matrix(sumfit$coefficients[2,])))
  out$tissue <- ts
  return(out)
})
lm.2.res <- do.call(rbind, lm.2)
lm.2.res$FDR <- stats::p.adjust(lm.2.res$`Pr(>|t|)`)
lm.2.res[order(lm.2.res$FDR),]

#pdf('derived/plots/210624_500_similarity_score_vs_rna_prt.pdf',width = 8, height = 6)
tissues_selected <- c("Brain_Cerebellum", "Brain_Cortex", 'Liver', 'Nerve_Tibial', 'Heart_Ventricle', 'Muscle_Skeletal')
ggplot(df_sd[df_sd$tissue %in% tissues_selected,], aes(x=mean_score, y = prt, color = obs)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", col = "red") +
  facet_wrap(~tissue) +
  xlab('Similarity Score [0-1]') +
  ylab('Standard Deviation of Log2 Protein Expression') + 
  theme_grey() +
  ggtitle('Mean Similarity Score versus SD(Protein Expression)')

ggplot(df_sd[df_sd$tissue %in% tissues_selected,], aes(x=mean_score, y = rna, color = obs)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", col = "lightblue") +
  facet_wrap(~tissue) +
  xlab('Similarity Score [0-1]') +
  ylab('Standard Deviation of Log2 RNA Expression') + 
  theme_grey() +
  ggtitle('Mean Similarity Score versus SD(RNA Expression)')


df_ratio_sd <- do.call(rbind, lapply(tissues, function(ts){
  mat_prt_selected <- mat_prt[mat_prt$tissue == ts]
  mat_prt_selected$prt_rna_ratio <- mat_prt_selected$prt - mat_prt_selected$rna
  aggr.ratio <- aggregate(prt_rna_ratio ~ mean_score, data = mat_prt_selected, sd1)
  aggr.lens <- aggregate(prt_rna_ratio ~ mean_score, data = mat_prt_selected, length)
  colnames(aggr.ratio)[2] <- 'prt_rna_ratio'
  colnames(aggr.lens)[2] <- 'obs'
  aggr <- merge(aggr.ratio, aggr.lens)
  aggr$tissue <- ts
  return(aggr)
}))

tissues_selected <- c("Brain_Cerebellum", "Brain_Cortex", 'Liver', 'Nerve_Tibial')
ggplot(df_ratio_sd[df_ratio_sd$tissue %in% tissues_selected,], aes(x=mean_score, y = prt_rna_ratio, color = obs)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", col = "red") +
  facet_wrap(~tissue) +
  xlab('Similarity Score [0-1]') +
  ylab('Standard Deviation of log(PRT)-log(RNA) expression') + 
  theme_grey()


# do sets of similar UTRs tend to more often be HI / TS
collins2020 <- setDT(read_xlsx('~/Projects/08_genesets/genesets/data/dosage/Collins2021medrxvic_supplementary_data.xlsx', 13))
complexity <- fread('derived/tables/210615_MANE.v0.93.UTR_features.txt', sep = '\t')
#constraints <- fread('~/Projects/08_genesets/genesets/data/gnomad/karczewski2020/supplementary_dataset_11_full_constraint_metrics.tsv') #transcript-level 

lapply(top_pairs, function(gt){
  
  current_genes <- unique(c(pairs[[gt]]$g1, pairs[[gt]]$g2))
  complexity[complexity$ensgid %in% current_genes]
  
})


#count <- length(top_pairs)
#cons <- lapply(top_pairs, function(gt){
#  out <- NULL
#  count <<- count -1
#  
#  current_genes_hgnc <- complexity[complexity$ensgid %in% current_genes,]$gene_symbol
#  if (any(collins2020$gene %in% current_genes_hgnc)) out <- data.frame(collins2020[collins2020$gene %in% current_genes_hgnc,], mean_score = mean(pairs[[gt]]$score, na.rm = T)) 
#  return(out)
#})
#
#mat_cons <- do.call(rbind, cons)
#
#plot(mat_cons$mean_score, mat_cons$pHI)
graphics.off()







