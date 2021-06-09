# author: Frederik Heymann Lassen, 03-June-2021
# description: Investigate whether the correlation between protein and RNA expression 
# are a result of driving 

library(readxl)
library(data.table)
library(ggplot2)
library(ggrepel)
library(ggsignif)
library(biomaRt)

# get UTR lengths
utr <- fread('extdata/biomart-5prime-utr-canonical-protein_coding.txt')
utr <- utr[!is.na(utr$`5' UTR start`) & !is.na(utr$`5' UTR end`),]
utr$UTR_length <- utr$`5' UTR end`- utr$`5' UTR start`
utr <- utr[,c(1,8)]
colnames(utr)[1] <- 'ensembl_id'

# load data
table_rna <- read_excel('extdata/prt_rna/table_s3.xlsx', sheet = 'B RNA tissue median', skip = 1) # C RNA standard z-score
table_prt <- read_excel('extdata/prt_rna/table_s2.xlsx', sheet = 'E protein tissue median', skip = 1) # F protein standard z-score
table_cor <- read_excel('extdata/prt_rna/table_s4.xlsx', sheet = 'B concordance comparison', skip = 1) # F protein standard z-score

# check table cor instead
head(table_cor)
dis_rna <- table_cor == 'discordance_rna_high'; sum(dis_rna, na.rm = T)
dis_prt <- table_cor == 'discordance_prt_high'; sum(dis_prt, na.rm = T)
concordance <- table_cor == 'concordance'; sum(concordance, na.rm = T)

# extract genes x categories that are disconcordant
df <- table_cor[,c(1,4,6:37)]
long <- melt(setDT(df), id.vars = c(1,2))
long <- long[!is.na(long$value),]
colnames(long)[3] <- c('tissue')
levels(long$tissue) <- unlist(lapply(strsplit(levels(long$tissue), split = '\\_'), function(x) paste(unlist(strsplit(x[3], split = ' ')), collapse = '_')))
long <- merge(long, utr, all.x = T)
long$value <- factor(long$value)
str(long)
long$UTR_length <- NULL
#fwrite(long, 'derived/tables/210609_prt_rna_binary_concordance.txt', sep = '\t')

# merge with constraints
constr <- constraints[,c('gene','loeuf','pTS','pHI')]
colnames(constr)[1] <- 'hgnc_symbol'
long <- merge(long,constr, by = 'hgnc_symbol', all.x = T)


# compare 5' UTR lengths bween concordant / discordant across genes stratified by tissue
#max(log10(na.omit(long$UTR_length)))
#min(log10(na.omit(long$UTR_length)))
top = 2.1 # 3.4
comparisons <- list(c("concordance","discordance_prt_high"), c("concordance","discordance_rna_high"), c("discordance_prt_high","discordance_rna_high"))
p1 <- ggplot(long, aes(x=value,y=pHI, fill = value)) +
  geom_jitter(size = 0.5, alpha = 0.5) +
  geom_violin(alpha = 0.8) +
  xlab('Concordance group') +
  #ylab("LOEUF") +
  #ggsignif::geom_signif(comparisons = comparisons ,map_signif_level=TRUE, y_position = c(top, top+0.4, top+0.8)) +
  facet_wrap(~tissue) +
  theme_bw() +
  #ylim(c(0,top+1.2)) +
  theme(axis.text=element_text(size=12), axis.text.x=element_blank(), axis.title=element_text(size=14,face="bold")) +
  ggtitle('Comparison of LOEUFs using Human Proteome Map Categories (wilcox.test)')

p1

ggsave('derived/plots/210604_condordance_vs_loeuf_all.pdf',width = 14, height = 16)  


# check constraints




