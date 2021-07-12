
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

get_oe <- function(dfs, obs, reps = 1:35, fun = mean){
  matrices <- lapply(reps, function(i){
    name <- paste0('rep',i,'.expt')
    expt <- colSums(dfs[, get(name, dfs), with = F])
    found <- colSums(obs)
    return(found/expt)
  })
  
  return(apply(simplify2array(matrices), 1, fun))

}
