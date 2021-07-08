# Do a GO enrichment analysis
devtools::load_all()
library(genoppi)
library(cowplot)
d <- fread('derived/tables/210708_MANE.v0.95.UTR_features.txt', sep = '\t')

gtex_cats <- fread('~/Projects/08_genesets/genesets/data/gtex/GTEX.tstat.categories.genoppi.csv')


### helper starts

format_tissue_enrichment <- function(res){
  
  # pre data
  mat_res <- cbind(as.character(res[[1]][,1]), do.call(cbind, lapply(res, function(x) x[,2])))
  mat_res <- cbind(as.data.frame(mat_res[,1]), as.data.frame(apply(mat_res[,-1], 2, as.numeric)))
  colnames(mat_res) <- c('Tissue.genoppi',paste0(1:length(res)))
  mat_res <- melt(mat_res)
  
  # prettify
  mat_res$significant <- mat_res$value < 0.05
  mat_res$bonfsig <- mat_res$value < 0.05/53
  mat_res$label <- ifelse(mat_res$bonfsig, '**', ifelse(mat_res$significant, '*',''))
  mat_res <- merge(mat_res, gtex_cats)
  
  # set order
  mat_res$order <- order(mat_res$Tissue.category.for.display)
  mat_res$Tissue <- factor(mat_res$Tissue, levels = unique(mat_res$Tissue[mat_res$order]))
  
  #mat_res$Tissue.category.for.display <- factor(mat_res$Tissue.category.for.display)
  tmp <- mat_res
  tmp <- tmp[,c('Tissue','Tissue.category.for.display')]
  tmp <- tmp[!duplicated(tmp$Tissue),]
  tmp$order <- order(tmp$Tissue.category.for.display)
  hlines <- cumsum(unlist(lapply(unique(tmp$Tissue.category.for.display), function(x){
    nrow(tmp[tmp$Tissue.category.for.display %in% x,])
  }))) 
  
  # manually add tissue category
  mat_res$test <- paste0(mat_res$Tissue, ' (',mat_res$Tissue.category.for.display,')')
  mat_res$test <- factor(mat_res$test, levels =  unique(mat_res$test[mat_res$order]))
  return(mat_res)
}





### helper end



#########
# uORFS #
########

# calc tissue enrichment 
res <- lapply(0:5, function(i){
  print(i)
  uorfs <- data.frame(gene = d$gene_symbol, significant = d$u5_ORF == i)
  print(sum(uorfs$significant))
  result <- genoppi::lapply_calc_hyper(uorfs, gtex_rna, col.by = 'tissue', intersectN = F)
  return(result[,c('list_name','pvalue')])
})


# What tissue is enriched?
mat_wo_uorf <- format_tissue_enrichment(res[1:2])
mat_wo_uorf <- mat_wo_uorf[mat_wo_uorf$variable == 1,]
mat_wo_uorf$variable <- factor(0)
mat_uorf <- format_tissue_enrichment(res[2:9])
mat_uorf <- format_tissue_enrichment(res[2:3])


p1 <- ggplot(mat_uorf, aes(x=variable, y = test,fill = log10(value), label = label)) +
  geom_tile(show.legend = F) +
  geom_text() +
  scale_fill_gradient(low='red',high='white') +
  xlab('N uORFs in gene') +
  ylab('GTEx (RNA expression)') +
  geom_hline(yintercept = c(2.5, 10.5, 23.5, 30.5, 34.5), linetype = 'dashed', color = 'black') +
  theme_bw()

p2 <- ggplot(mat_wo_uorf, aes(x=variable, y = test,fill = log10(value), label = label)) +
  geom_tile(show.legend = F) +
  geom_text() +
  scale_fill_gradient(low='blue',high='white') +
  xlab('') +
  ylab('GTEx (RNA expression)') +
  geom_hline(yintercept = c(2.5, 10.5, 23.5, 30.5, 34.5), linetype = 'dashed', color = 'black') +
  theme_bw() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
  

cowplot::plot_grid(plotlist = list(p1, p2),rel_widths = c(0.9,0.1))



# calc tissue enrichment 
res <- lapply(0:2, function(i){
  print(i)
  uorfs <- data.frame(gene = d$gene_symbol, significant = d$u5_oORF == i)
  print(sum(uorfs$significant))
  result <- genoppi::lapply_calc_hyper(uorfs, gtex_rna, col.by = 'tissue', intersectN = F)
  return(result[,c('list_name','pvalue')])
})


# What tissue is enriched?
mat_wo_oorfs <- format_tissue_enrichment(res[1:2])
mat_wo_oorfs <- mat_wo_oorfs[mat_wo_oorfs$variable == 1,]
mat_wo_oorfs$variable <- factor(0)
mat_oorf <- format_tissue_enrichment(res[2:3])

p1 <- ggplot(mat_oorf, aes(x=variable, y = test,fill = log10(value), label = label)) +
  geom_tile(show.legend = F) +
  geom_text() +
  scale_fill_gradient(low='red',high='white') +
  xlab('oORFs (n)') +
  ylab('GTEx (RNA expression)') +
  geom_hline(yintercept = c(2.5, 10.5, 23.5, 30.5, 34.5), linetype = 'dashed', color = 'black') +
  theme_bw()

p2 <- ggplot(mat_wo_oorf, aes(x=variable, y = test,fill = log10(value), label = label)) +
  geom_tile(show.legend = F) +
  geom_text() +
  scale_fill_gradient(low='blue',high='white') +
  xlab('') +
  ylab('GTEx (RNA expression)') +
  geom_hline(yintercept = c(2.5, 10.5, 23.5, 30.5, 34.5), linetype = 'dashed', color = 'black') +
  theme_bw() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())


cowplot::plot_grid(plotlist = list(p1, p2),rel_widths = c(0.9,0.1))



# calc tissue enrichment 
res <- lapply(0:5, function(i){
  print(i)
  uorfs <- data.frame(gene = d$gene_symbol, significant = d$u5_NTE == i)
  print(sum(uorfs$significant))
  result <- genoppi::lapply_calc_hyper(uorfs, gtex_rna, col.by = 'tissue', intersectN = F)
  return(result[,c('list_name','pvalue')])
})


# What tissue is enriched?
mat_wo_nte <- format_tissue_enrichment(res[1:2])
mat_wo_nte <- mat_wo_nte[mat_wo_nte$variable == 1,]
mat_wo_nte$variable <- factor(0)
mat_nte <- format_tissue_enrichment(res[2:6])

p1 <- ggplot(mat_nte, aes(x=variable, y = test,fill = log10(value), label = label)) +
  geom_tile(show.legend = F) +
  geom_text() +
  scale_fill_gradient(low='red',high='white') +
  xlab('NTEs (n)') +
  ylab('GTEx (RNA expression)') +
  geom_hline(yintercept = c(2.5, 10.5, 23.5, 30.5, 34.5), linetype = 'dashed', color = 'black') +
  theme_bw()

p2 <- ggplot(mat_wo_nte, aes(x=variable, y = test,fill = log10(value), label = label)) +
  geom_tile(show.legend = F) +
  geom_text() +
  scale_fill_gradient(low='blue',high='white') +
  xlab('NTEs (n) ') +
  ylab('GTEx (RNA expression)') +
  geom_hline(yintercept = c(2.5, 10.5, 23.5, 30.5, 34.5), linetype = 'dashed', color = 'black') +
  theme_bw() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())


cowplot::plot_grid(plotlist = list(p1, p2),rel_widths = c(0.9,0.1))




