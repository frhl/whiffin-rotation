# constraint modelling in 5' UTR

library(data.table)




d_obs <- fread('derived/210701_MANE.GRCh38.v0.95_codons_obs.csv', sep = ',')
d_expt <- fread('derived/210701_MANE.GRCh38.v0.95_codons_expt.csv', sep = ',')
#UTR <- "5' UTR"

#d_obs_ci <- fread('derived/210705_MANE.GRCh38.v0.95_three_prime_utr_codons_expt_ci.csv', sep = ',')
#d_obs <- fread('derived/210705_MANE.GRCh38.v0.95_three_prime_utr_codons_obs.csv', sep = ',')
#d_expt <- fread('derived/210705_MANE.GRCh38.v0.95_three_prime_utr_codons_expt.csv', sep = ',')
#UTR <- "3' UTR"

    
    
  
pdf('derived/plots/210712_5UTR_features.pdf', width = 5, height = 4)


# setup colors
library(RColorBrewer)
color = brewer.pal(6, 'Set2') 
names(color) <- c('ATG','CGA','TAG','TGA','TAA')
color_scale <- scale_colour_manual(name = "codon",values = color)


# merge observed / expected
mrg <- merge(d_expt, d_obs)
mrg$ensgid_version <- NULL
mrg$ensgid <- NULL

# observed versus expected
index_obs <- 3:66
index_expt <- 67:(67+63)

obs <- mrg[,index_obs, with = F] 
expt <- mrg[,index_expt, with = F] 




features <- fread('derived/tables/210709_MANE.v0.95.UTR_features.txt', sep = '\t')

# finally merge data and check
#combined <- merge(features, mrg, by = c('enstid_version','ensgid'))
#combined$u5_proximity <- as.numeric(unlist(lapply(combined$u5_ORF_cds_proximity, function(x) min(unlist(strsplit(x, split = ';'))))))
#combined <- combined[!is.na(combined$u5_proximity)]
#sum(!is.na(combined$u5_proximity))

# create decile for proximity
#combined$u5_proximity_pct <- combined$u5_proximity / combined$u5_len
#decile_seq <- seq(0,1,by=0.1)
#deciles <- quantile(combined$u5_proximity_pct, probs = decile_seq)
#combined$u5_proximity_decile <- cut(combined$u5_proximity_pct, deciles)
#levels(combined$u5_proximity_decile) <- decile_seq*10
#combined <- combined[!is.na(combined$u5_proximity_decile)]
#plot(combined$u5_proximity_decile, combined$u5_proximity)

#list_proximity <- lapply(unique(combined$u5_proximity_decile), function(decile){
#  row_bool <- combined$u5_proximity_decile == decile
#  obs_get <- colSums(combined[row_bool, get('obs',combined), with = F])
#  expt_get <- colSums(combined[row_bool, get('expt',combined), with = F])
#  oe <- obs_get/expt_get
#  oe_names <- gsub('obs\\.','',names(oe))
#  df <- data.table(decile = decile, codons = oe_names, oe = oe)
#  return(df)
#})

#proximity <- do.call(rbind, list_proximity)
#proximity <- proximity[proximity$codons %in% c('ATG')]# ,'CGA','TAG','TAA', 'CTG'),]
#ggplot(proximity, aes(x=decile, y = oe, fill = codons, group = codons)) +
#  geom_bar(stat='identity')
#  #geom_point() + 
#  #geom_smooth(method = 'lm',)


###############
# gene length #
###############

combined <- merge(features, mrg, by = c('enstid_version'))
combined <- combined[combined$u5_len > 30,]
enstids <- combined$enstid_version

# create decile for proximity
decile_seq <- seq(0,1,by=0.1)
deciles <- quantile(combined$u5_len, probs = decile_seq)
combined$u5_len_decile <- cut(combined$u5_len, deciles)
levels(combined$u5_len_decile) <- decile_seq*10
combined <- combined[!is.na(combined$u5_len_decile)]
plot(combined$u5_len_decile, combined$u5_len)

list_len <- lapply(unique(combined$u5_len_decile), function(decile){
  row_bool <- combined$u5_len_decile == decile
  obs_get <- colSums(combined[row_bool, get('obs',combined), with = F])
  expt_get <- colSums(combined[row_bool, get('expt',combined), with = F])
  oe_sem <- apply(combined[row_bool, get('expt',combined), with = F], 2, sem)
  oe <- obs_get/expt_get
  oe_names <- gsub('obs\\.','',names(oe))
  df <- data.table(decile = decile, codons = oe_names, oe = oe, sem = oe_sem)
  return(df)
})


lens <- do.call(rbind, list_len)
lens$codons <- as.factor(lens$codons)
lens <- lens[lens$codons %in% c('ATG','CGA','TAG','TGA','TAA'),]
ggplot(lens, aes(x=decile, y = oe, color = codons, group = codons, ymax = oe + sem, ymin = oe - sem)) +
  geom_point() + 
  geom_errorbar(width = 0.05, position = 'dodge') +
  geom_smooth(method = 'lm') +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  xlab('UTR Length (Decile)') +
  ylab('Observed / Expected') +
  theme_bw() +
  color_scale +
  ggtitle(paste0(UTR," UTR length"))

##############
# GC content #
##############

combined <- merge(features, mrg, by = c('enstid_version'))
combined <- combined[combined$u5_len > 30,]

# create decile for proximity
decile_seq <- seq(0,1,by=0.1)
deciles <- quantile(combined$u5_GC, probs = decile_seq)
combined$u5_gc_decile <- cut(combined$u5_GC, deciles)
levels(combined$u5_gc_decile) <- decile_seq*10
combined <- combined[!is.na(combined$u5_gc_decile)]
plot(combined$u5_gc_decile, combined$u5_gc)

list_gc <- lapply(unique(combined$u5_gc_decile), function(decile){
  row_bool <- combined$u5_gc_decile == decile
  obs_get <- colSums(combined[row_bool, get('obs',combined), with = F])
  expt_get <- colSums(combined[row_bool, get('expt',combined), with = F])
  oe_sem <- apply(combined[row_bool, get('expt',combined), with = F], 2, sem)
  oe <- obs_get/expt_get
  oe_names <- gsub('obs\\.','',names(oe))
  df <- data.table(decile = decile, codons = oe_names, oe = oe, sem = oe_sem)
  return(df)
})

gcs <- do.call(rbind, list_gc)
gcs <- gcs[gcs$codons %in% c('ATG','CGA','TAG','TAA','TGA'),]
ggplot(gcs, aes(x=decile, y = oe, color = codons, group = codons, ymin = oe-sem, ymax = oe+sem)) +
  geom_point() + 
  geom_errorbar(width = 0.05, alpha = 0.3) +
  geom_smooth(method = 'lm') +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  xlab('UTR GC content (Decile)') +
  ylab('Observed / Expected') +
  theme_bw() +
  color_scale +
  ggtitle(paste0(UTR," UTR GC content"))

#############
# ATG Kozak #
#############

#combined <- merge(features, mrg, by = c('enstid_version'))
#combined <- combined[combined$u5_len > 30,]

# create decile for proximity
#decile_seq <- seq(0,1,by=0.1)
#deciles <- quantile(combined$u5_len, probs = decile_seq)
#combined$u5_len_decile <- cut(combined$u5_len, deciles)
#levels(combined$u5_len_decile) <- decile_seq*10
#combined <- combined[!is.na(combined$u5_len_decile)]

#list_len_kozak <- lapply(unique(combined$u5_len_decile), function(decile){
#  kozak <- do.call(rbind, lapply(1:3, function(k){
#    row_bool <- combined$u5_len_decile == decile & combined$u5_ORF_kozak == k
#    obs_get <- colSums(combined[row_bool, get('obs',combined), with = F])
#    expt_get <- colSums(combined[row_bool, get('expt',combined), with = F])
#    oe <- obs_get/expt_get
#    oe_names <- gsub('obs\\.','',names(oe))
#    df <- data.table(decile = decile, codons = oe_names, oe = oe, kozak = k)
#  }))
#  return(kozak)
#})

#kozak <- do.call(rbind, list_len_kozak)
#bool <- kozak$codons %in% 'ATG'
#kozak <- kozak[bool,]
#kozak$kozak <- as.factor(kozak$kozak)
#ggplot(kozak, aes(x=decile, y = oe, color = kozak, group = kozak)) +
#  geom_point() + 
#  geom_smooth(method = 'lm') +
#  geom_hline(yintercept = 1, linetype = 'dashed') +
#  xlab('Length (Decile)') +
#  ylab('Observed / Expected') +
#  theme_bw() +
#  ggtitle(paste0(UTR," UTR length - ATG depletion"))

############
# genesets #
############

get_genelist <- function(){
  paths <- list.files('~/Projects/10_mcarthur_genelists/gene_lists/lists', full.names = T)
  lst <- lapply(paths, function(x) fread(x, header = F)$V1)
  names(lst) <- tools::file_path_sans_ext(basename(paths))
  lst[['Transcript_annotation_Haploinsufficient_genes_61.tsv']] <- NULL
  df <- as.data.frame(stack(lst))
  colnames(df) <- c('gene_symbol','geneset')
  return(df)
}


#genelists <- list(
#  olfactory = data.table(genelist = 'olfactory receptors', gene_symbol = fread('~/Projects/10_mcarthur_genelists/gene_lists/lists/olfactory_receptors.tsv', header = F)$V1),
#  mgi = data.table(genelist = 'olfactory receptors', gene_symbol = fread('~/Projects/10_mcarthur_genelists/gene_lists/lists/mgi_essential.tsv', header = F)$V1),
#  ar1 = data.table(genelist = 'olfactory receptors', gene_symbol = fread('~/Projects/10_mcarthur_genelists/gene_lists/lists/all_ar.tsv', header = F)$V1)
#)

mcarthur<- get_genelist()

for (geneset in unique(genelist$geneset)){
  
  combined <- merge(features, mrg, by = c('enstid_version'))
  combined$genelist <- combined$gene_symbol %in% mcarthur$gene_symbol[mcarthur$geneset %in% geneset]
  list_genelist <- lapply(c(TRUE,FALSE), function(bool){
    row_bool <- combined$genelist == bool
    obs_get <- colSums(combined[row_bool, get('obs',combined), with = F])
    expt_get <- colSums(combined[row_bool, get('expt',combined), with = F])
    oe_sem <- apply(combined[row_bool, get('expt',combined), with = F], 2, sem)
    oe <- obs_get/expt_get
    oe_names <- gsub('obs\\.','',names(oe))
    df <- data.table(genelist = bool, codons = oe_names, oe = oe, sem = oe_sem)
  })
  
  glist <- do.call(rbind, list_genelist)
  glist <- glist[glist$codons %in% c('ATG','CGA','TAG','TAA','TGA'),]
  p1 <- ggplot(glist, aes(x=genelist, y = oe, color = codons, group = codons, ymin = oe - sem, ymax = oe + sem)) +
    geom_point() + 
    geom_errorbar(width = 0.05, alpha = 0.7) +
    geom_line(linetype = 'dashed') +
    geom_hline(yintercept = 1, linetype = 'dashed') +
    xlab('In genelist') +
    ylab('Observed / Expected') +
    ggtitle(paste0(UTR," Genelist: ", geneset)) +
    color_scale +
    theme_bw()
  print(p1)

}

graphics.off()






