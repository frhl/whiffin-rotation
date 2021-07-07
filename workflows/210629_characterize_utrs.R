#

# get protein / RNA 
#d <- fread('derived/tables/210629_MANE.v0.95.UTR_features.txt', sep = '\t')
d <- fread('derived/tables/210706_MANE.v0.95.UTR_features.txt', sep = '\t')


# non parametric bootstrap
bootstrap_mean_ci <- function(x, k = 100, func = function(x) sd(x, na.rm = T), probs = c(0.025,0.5,0.975)){
  simsamples <- replicate(k, sample(x, replace = T))
  simmedians <- apply(simsamples, 2, func)
  outdf <- t(as.data.frame(quantile(simmedians, probs = probs)))
  rownames(outdf) <- NULL
  return(outdf)
}

pasteci <- function(x, digits = 2){
  x <- round(x, digits = digits)
  paste0(x[2],' [',x[1],';',x[3],']')
} 

eval_u5_features <- function(k = 100, f){
  data.frame(
    transcripts = length(unique(d$ensgid_version)),
    lens = pasteci(bootstrap_mean_ci(d$u5_len, k, f)),
    augs = pasteci(bootstrap_mean_ci(d$u5_AUG, k, f)),
    start_to_cap = pasteci(bootstrap_mean_ci(d$u5_ORF_cap_to_start, k, f)),
    uorf = pasteci(bootstrap_mean_ci(d$u5_ORF, k, f)),
    oorf =  pasteci(bootstrap_mean_ci(d$u5_oORF_altered_cds, k, f)),
    oorf =  pasteci(bootstrap_mean_ci(d$u5_oORF_inframe, k, f)),
    gc =  pasteci(bootstrap_mean_ci(d$u5_GC, k, f))
  )
}

count_u5_features <- function(){
  data.frame(
    transcripts = length(unique(d$ensgid_version)),
    augs = sum(d$u5_AUG > 0),
    uorf = sum(d$u5_ORF > 0),
    oorf =  sum(d$u5_oORF_altered_cds > 0),
    oorf_inframe =  sum(d$u5_oORF_inframe > 0),
    kozak_strong = sum(d$u5_oORF_kozak == 3),
    kozak_moderate = sum(d$u5_oORF_kozak == 2),
    kozak_weak = sum(d$u5_oORF_kozak == 1)
  )
}

count_to_pct <- function(x) { paste0(x, ' (',round(100*(x / x[[1]]),2), '%)')}

sds <- eval_u5_features(k = 1000, f = function(x) sd(x, na.rm = T))
mus <- eval_u5_features(k = 1000, f = function(x) mean(x, na.rm = T))
outmat <- rbind(mus, sds)
rownames(outmat) <- c('Means','SDs')

count_to_pct(count_u5_features())


#





