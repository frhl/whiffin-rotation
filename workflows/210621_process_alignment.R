

# read alignment
d <- fread(file = 'derived/tables/210623_500_G2P-confirmedprobable_five_prime_utr_alignment_NW.txt', sep = '\t')
d <- as.data.table(melt(d, id.vars = 'V1'))
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
score_match <- 3
score_gap <- 1
final$cur_score <- final$m*score_match - final$g*score_gap
final$max_score <- final$l*score_match
final$score <- final$cur_score / final$max_score

# select top 20%
top <- quantile(final$score, probs = 0.99)
final[final$score == max(final$score)]
final[final$score > top]

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
top_expression <- head(rev(sort(na.omit(scores))), n = 200)
top_pairs <- names(top_expression)

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
  
  fit <- lm(prt ~ group, data =  mat_prt_selected)
  sumfit <- summary(fit)
  out <- as.data.frame(t(as.matrix(sumfit$coefficients[2,])))
  out$tissue <- ts
  out$mean_score <- unique(mat_prt_selected)
  return(out)
})
lm.1.res <- do.call(rbind, lm.1)
lm.1.res$FDR <- stats::p.adjust(lm.1.res$`Pr(>|t|)`)


ggplot(mat_prt[mat_prt$tissue %in% 'Adrenal_Gland'], aes(x=group, y = prt, group = group)) +
  geom_boxplot() +
  stat_smooth(method = "lm", col = "red") 





