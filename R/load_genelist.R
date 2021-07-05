# temprory scripts for plotting and retriving genesets

load_mcarthur_list <- function(name, newname, compare){

  path <- paste0('~/Projects/10_mcarthur_genelists/gene_lists/lists/',name)
  if (!file.exists(path)) stop(paste0(path, 'does not exist!'))
  compare$geneset = compare$gene %in% fread(path, header = F)$V1
  d <- as.data.frame(table(compare$geneset, compare$decile))
  colnames(d) <- c('geneset','decile','count')
  d$type <- newname
  d <- d[d$decile != 100,]
  return(d)
}

load_mcarthur_list_loeuf <- function(name, newname, compare){

  d_test0 <- load_mcarthur_list(name, newname, compare)
  d_test0$how <- 'depletion score'
  
  path <- paste0('~/Projects/10_mcarthur_genelists/gene_lists/lists/',name)
  compare[[newname]] = compare$gene %in% fread(path, header = F)$V1
  d_test1 <- as.data.frame(table(compare[[newname]], compare$decile_loeuf))
  colnames(d_test1) <- c('geneset','decile','count')
  d_test1$type <- newname
  d_test1 <- d_test1[d_test1$decile != 100,]
  d_test1$how = 'LOEUF'
  
  com <- rbind(d_test0, d_test1)
  com <- com[com$geneset == TRUE,]
  
  return(com)
  
}

plotit <- function(d){
  
  p <- ggplot(d, aes(y=count, x = decile, fill = how)) +
          geom_bar(stat='identity', position = 'dodge', color = 'black') +
          xlab('Decile') + 
          ylab('Count') +
          ggtitle(unique(d$type))
  
  print(p)
}