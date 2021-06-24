# author: frederik heymann (21-06-09)
# description: compare complexity of UTR across different tissue/protein expression tables

devtools::load_all()

# load expression and UTR data
expression <- fread('derived/tables/210609_prt_rna_numerical.txt', sep = '\t')
complexity <- fread('derived/tables/210615_MANE.v0.93.UTR_features.txt', sep = '\t')
dt <- merge(complexity, expression, by.x = 'ensgid',by.y = 'gene.id')
dt$rna_std <- (dt$rna - mean(dt$rna, na.rm = T))/sd(dt$rna, na.rm = T)
dt$prt_std <- (dt$prt - mean(dt$prt, na.rm = T))/sd(dt$prt, na.rm = T)
dt <- dt[dt$prt != 0, ]
dt$prt_rna_ratio <-  dt$prt-dt$rna
dt$rna_prt_ratio <-  dt$rna-dt$prt

# function for generating empircal CDF
emp_cdf <- function(x){
  x <-  sort(x)
  y <- seq(1, length(x)) / length(x)
  data.frame(x=x, y=y)
}

# count occurences
count_how <- function(df){
  counts <- table(df$how)
  paste0(names(counts),' N=',counts, collapse = '\n')
}


# kozak strengths
eval_kozak_ecdf <- function(d, col = 'prt'){
  
  cdf0 <- data.frame(emp_cdf(d[[col]][d$u5_ORF_kozak == 0]), how = 'uORF Kozak - 0')
  cdf1 <- data.frame(emp_cdf(d[[col]][d$u5_ORF_kozak == 1]), how = 'uORF Kozak - 1')
  cdf2 <- data.frame(emp_cdf(d[[col]][d$u5_ORF_kozak == 2]), how = 'uORF Kozak - 2')
  cdf3 <- data.frame(emp_cdf(d[[col]][d$u5_ORF_kozak == 3]), how = 'uORF Kozak - 3')
  df <- rbind(cdf0, cdf1,  cdf2, cdf3)
  df$how <- factor(df$how)
  return(df)
  
}

# oorfs
eval_oorf_ecdf <- function(d, col = 'prt'){
  
  cdf0 <- data.frame(emp_cdf(d[[col]][d$u5_oORF_altered_cds == 0]), how = 'no oORF')
  cdf1 <- data.frame(emp_cdf(d[[col]][d$u5_oORF_altered_cds != 0]), how = 'oORF')
  df <- rbind(cdf0, cdf1)
  df$how <- factor(df$how)
  return(df)
  
}

# no uORF versus uORF
eval_uorf_ecdf <- function(d, col = 'prt'){
  
  cdf1 <- data.frame(emp_cdf(d[[col]][d$u5_ORF == 0]), how = 'no uORF')
  cdf2 <- data.frame(emp_cdf(d[[col]][d$u5_ORF != 0]), how = 'uORF')
  df <- rbind(cdf1, cdf2)
  df$how <- factor(df$how)
  return(df)
  
}

# uORF lengths
eval_n_uorf_ecdf <- function(d, col = 'prt'){
  
  tabl <- table(d$u5_ORF)
  tabl <- as.data.frame(tabl[tabl > 50])
  df <- lapply(1:nrow(tabl), function(i)  data.frame(emp_cdf(d[[col]][d$u5_ORF == tabl$Var1[i]]), how = paste(tabl$Var1[i], 'uORFs')))
  df <- do.call(rbind, df)
  df$how <- factor(df$how)
  return(df)
  
}

ggplot_ecdf_protein <- function(df){
  
  ggplot(df, aes(x = x, y = y, color = how)) +
    geom_point(size = 0.5) +
    geom_line(size = 0.3, alpha = 0.8) + 
    annotate(geom = 'text',x = 2, y = 0.1, label = count_how(df), hjust = 0) +
    theme_linedraw() +
    xlab('Protein Expression (Log)') +
    ylab('Fraction of Genes')
  
}

compare_rna <- function(df, col = 'u5_ORF'){
  bool <- df[[col]] == 0  
  t.test(df$rna_std[bool], df$rna_std[!bool])
  # x: no uORF
  # y: with uORF
}


# mRNA versus protein
tissues <- c('Liver','Brain_Cerebellum','Brain_Cortex','Lung')


plot(density(dt$prt), col = 'red')
lines(density(na.omit(dt$rna)), col = 'blue')
plot(density(na.omit(dt$rna_prt_ratio)))


#d_liver_prt <- eval_uorf_ecdf(dt[dt$tissue == tissues[1],],'prt')
#d_liver_prt$what <- 'prt'
#d_liver_rna <- eval_uorf_ecdf(dt[dt$tissue == tissues[1],],'rna')
#d_liver_rna$what <- 'rna'
#df <- rbind(d_liver_prt, d_liver_rna)



## check uORFs

# cumulative expression of genes containing uORFs versus no uORFS
tissues <- c('Liver','Brain_Cerebellum','Brain_Cortex','Lung')
d_liver <- eval_uorf_ecdf(dt[dt$tissue == tissues[1],])
d_liver_rna <- eval_uorf_ecdf(dt[dt$tissue == tissues[1],],'rna')
d_cerebellum <- eval_uorf_ecdf(dt[dt$tissue == tissues[2],])
d_cerebellum_rna <- eval_uorf_ecdf(dt[dt$tissue == tissues[2],],'rna')
d_cortex <- eval_uorf_ecdf(dt[dt$tissue == tissues[3],])
d_cortex_rna <- eval_uorf_ecdf(dt[dt$tissue == tissues[3],],'rna')
d_lung <- eval_uorf_ecdf(dt[dt$tissue == tissues[4],])
d_lung_rna <- eval_uorf_ecdf(dt[dt$tissue == tissues[4],],'rna')


lapply(list(d_liver, d_cerebellum, d_cortex, d_lung), function(j) ks.test(j$x, j$y ))

#compare_rna(dt[dt$tissue == tissues[1]])
#compare_rna(dt[dt$tissue == tissues[2]])
#compare_rna(dt[dt$tissue == tissues[3]])

#  colors
myColors <- c('grey','red')
names(myColors) <- levels(df$how)
colScale <- scale_colour_manual(name = "how",values = myColors)

myColors <- c('grey','lightblue')
names(myColors) <- levels(df$how)
colScale_blue <- scale_colour_manual(name = "how",values = myColors)

# uorf protein expression
pdf('derived/plots/220624_ecdf_uorf_rna_prt_tissues.pdf', width = 8, height = 3)
ggplot_ecdf_protein(d_liver) + colScale + ggtitle('Liver (Protein)','') + xlab('Median Protein expression (Log2)')
ggplot_ecdf_protein(d_liver_rna) + colScale_blue + ggtitle('Liver (RNA)','') + xlab('Median RNA expression (Log2)')
ggplot_ecdf_protein(d_cerebellum) + colScale + ggtitle('Cerebellum (Protein)','') + xlab('Median Protein expression (Log2)')
ggplot_ecdf_protein(d_cerebellum_rna) + colScale_blue + ggtitle('Cerebellum (RNA)','') + xlab('Median RNA expression (Log2)')
ggplot_ecdf_protein(d_cortex) + colScale + ggtitle('Cortex (Protein)','') +xlab('Median Protein expression (Log2)')
ggplot_ecdf_protein(d_cortex_rna) + colScale_blue + ggtitle('Cortex (RNA)','') + xlab('Median RNA expression (Log2)')
graphics.off()

#pdf('derived/plots/220623_ecdf_uorf_all_tissues.pdf', width = 8, height = 3)
#pdf('derived/plots/220623_delta_test.pdf', width = 8, height = 3)
#for (tis in unique(dt$tissue)){
#  p <- ggplot_ecdf_protein( eval_uorf_ecdf(dt[dt$tissue == tis,], 'prt_rna_ratio')) + colScale + xlab('log(PRT/RNA)') + ggtitle(tis)
#  print(p)
#}
#graphics.off()

## check if more uORFs correspond to less protein

# cumulative expression of genes containing uORFs versus no uORFS
tissues <- c('Liver','Brain_Cerebellum','Brain_Cortex','Lung')
d_liver <- eval_n_uorf_ecdf(dt[dt$tissue == tissues[1],])
d_liver_ratio <- eval_n_uorf_ecdf(dt[dt$tissue == tissues[1],], 'rna_prt_ratio')
d_cerebellum <- eval_n_uorf_ecdf(dt[dt$tissue == tissues[2],])
d_cortex <- eval_n_uorf_ecdf(dt[dt$tissue == tissues[3],])
d_lung <- eval_n_uorf_ecdf(dt[dt$tissue == tissues[4],])
lapply(list(d_liver, d_cerebellum, d_cortex, d_lung), function(j) ks.test(j$x, j$y ))

# uorf protein expression
pdf('derived/plots/220623_ecdf_uorf_n_tissues.pdf', width = 8, height = 3)
ggplot_ecdf_protein(d_lung)  + ggtitle('Cortex')
ggplot_ecdf_protein(d_liver_ratio) +  ggtitle('Liver')
ggplot_ecdf_protein(d_cerebellum)  + ggtitle('Cerebellum')
ggplot_ecdf_protein(d_cortex) + ggtitle('Cortex')
graphics.off()


## check oORFs

# cumulative expression of genes containing uORFs versus no uORFS
tissues <- c('Liver','Brain_Cerebellum','Brain_Cortex')
d_liver <- eval_oorf_ecdf(dt[dt$tissue == tissues[1],])
d_cerebellum <- eval_oorf_ecdf(dt[dt$tissue == tissues[2],])
d_cortex <- eval_oorf_ecdf(dt[dt$tissue == tissues[3],])
lapply(list(d_liver, d_cerebellum, d_cortex), function(j) ks.test(j$x, j$y ))

#  colors
colfunc <- colorRampPalette(c("cyan", "red"))
myColors <- colfunc(4)
names(myColors) <- levels(df$how)
colScale <- scale_colour_manual(name = "how",values = myColors)

# uorf protein expression
pdf('derived/plots/220623_ecdf_oorf_tissues.pdf', width = 8, height = 4)
ggplot_ecdf_protein(d_liver) + colScale + ggtitle('Liver')
ggplot_ecdf_protein(d_cerebellum) + colScale + ggtitle('Cerebellum')
ggplot_ecdf_protein(d_cortex) + colScale + ggtitle('Cortex')
graphics.off()


## check Kozak

# cumulative expression of genes containing uORFs versus no uORFS
tissues <- c('Liver','Brain_Cerebellum','Brain_Cortex')
d_liver <- eval_kozak_ecdf(dt[dt$tissue == tissues[1],])
d_cerebellum <- eval_kozak_ecdf(dt[dt$tissue == tissues[2],])
d_cortex <- eval_kozak_ecdf(dt[dt$tissue == tissues[3],])
lapply(list(d_liver, d_cerebellum, d_cortex), function(j) ks.test(j$x, j$y ))

#  colors
colfunc <- colorRampPalette(c("cyan", "red"))
myColors <- colfunc(4)
names(myColors) <- levels(df$how)
colScale <- scale_colour_manual(name = "how",values = myColors)

# uorf protein expression
pdf('derived/plots/220623_ecdf_kozak_tissues.pdf', width = 8, height = 4)
ggplot_ecdf_protein(d_liver) + colScale + ggtitle('Liver','')
ggplot_ecdf_protein(d_cerebellum) + colScale + ggtitle('Cerebellum','')
ggplot_ecdf_protein(d_cortex) + colScale + ggtitle('Cortex','')
graphics.off()



