# constraint modelling in 5' UTR

library(data.table)



gather_expected <- function(path = 'download/210711_hpc_derived/'){
  paths <- list.files(path, pattern = "210711_MANE.GRCh38.v0.95_five_prime_utr_codons_expt_rep", full.names = T)
  dfs <- as.data.table(do.call(cbind, lapply(paths, function(p){ 
    rep <- unlist(strsplit(tools::file_path_sans_ext(basename(p)), split = '_'))[8]
    d <- fread(p)[,1:64]
    colnames(d) <- paste0(rep,'.',colnames(d))
    return(d)
  })))
  return(dfs)
}

get_oe <- function(expt, obs, fun = mean, in_reps = NULL){

  obs <- obs[,get('obs',obs),with = F]
  reps <- unique(gsub('rep','',unique(gsub('\\.expt\\.[A-Z]+','',colnames(expt)))))
  reps <- sort(as.numeric(reps))
  if (!is.null(in_reps)) reps <- reps[reps %in% in_reps]
  matrices <- lapply(reps, function(i){
    name <- paste0('rep',i,'.expt')
    expt <- colSums(expt[, get(name, expt), with = F])
    found <- colSums(obs)
    return(found/expt)
  })
  
  return(apply(simplify2array(matrices), 1, fun))
  
}




d_obs <- fread('derived/210701_MANE.GRCh38.v0.95_codons_obs.csv', sep = ',')
#d_expt <- gather_expected()

sem <- function(x) return(sd(x)/sqrt(length(x)))

get_upper <- function(x) quantile(x, probs = 0.975)
get_lower <- function(x) quantile(x, probs = 0.025)

res_upper <- get_oe(d_expt, d_obs, fun = get_upper, in_reps = 1:35)
res_lower <- get_oe(d_expt, d_obs, fun = get_lower, in_reps = 1:35)
res_mean <- get_oe(d_expt, d_obs, fun = mean, in_reps = 1:2)
codons <- gsub('obs\\.','',names(res_sem))
d <- data.frame(codons = codons, mean = res_mean, upper = res_upper, lower = res_lower)
#d$sem <- apply(expt, 2, sem)

ggplot(d, aes(x = reorder(codons, mean), y = mean, ymax = upper, ymin = lower)) +
  geom_point() +
  geom_errorbar() +
  xlab('Codons') +
  ylab('O/E (95% confidence interval of mean)') +
  ggtitle('Confidence interval for OE 7000 simulations per codon') + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 





