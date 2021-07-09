# Do a GO enrichment analysis
devtools::load_all()

# parallelize the workflow
library(data.table)
library(parallel)
library(doParallel)
cores <- detectCores()
registerDoParallel(cores)

d <- fread('derived/tables/210708_MANE.v0.95.UTR_features.txt', sep = '\t')
outdir = 'derived'

# setup analysis
dbs <- list(
  go_mf = goa_mf_table[c(1,3)],
  go_bp = goa_bp_table[c(1,3)],
  go_cc = goa_cc_table[c(1,3)]
)
dbs <- lapply(dbs, function(x){colnames(x) <- c('gene','pathway'); return(x)})


# without uORF and with uORF
for (db_name in names(dbs)){ 
  
  write(paste0('# running ',db_name,' (doParallel)'), stdout())
  result_go <- (foreach (i=0:4, .combine=rbind) %dopar% {
  
    write(paste0('checking uORF = ',i), stdout())
    db<- dbs[[db_name]]
    db$significant <- TRUE
    
    d_analysis <- data.frame(gene=d$gene_symbol, significant = d$u5_ORF == i)
    d_result <- lapply_calc_hyper(d_analysis, db, col.by = 'pathway', intersectN = F)
    d_result$u5_ORF <- i
    d_result$dataset <- db_name
    return(d_result)
    
  })
  
  outname <- paste0('210709_hypergeom_',db_name,'_uORF_analysis.txt')
  result_go <- as.data.frame(result_go)
  result_go$successInSampleGenes <- NULL
  outpath <- file.path(outdir,outname)
  write(paste0('writing ',outpath,'...'), stdout())
  write.table(result_go, outpath, sep = '\t',row.names = F)
  
}

write("Done.", stdout())







