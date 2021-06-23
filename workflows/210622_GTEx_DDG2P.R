# author: frederik heymann 21-06-16
# description Are genes with other organ involvement less likely to be brain specific in 
# expression? Also, are they more or less likely to have discordant RNA/protein 
# levels?

devtools::load_all()

# load DDG2P data
g2p <- fread('extdata/DDG2P_4_6_2021.csv')
g2p <- g2p[g2p$`DDD category` %in% c('confirmed','probable')] # 1809 genes
colnames(g2p) <- gsub(' ','_',colnames(g2p))

# load organ data
organ <- fread('extdata/210616_DDG2P_organ_matrix.csv')
organ <- data.frame(gene_symbol = organ$gene_symbol, 
                    organ_count = apply(organ[,-1],1 ,sum))

# combine with hillary finucance GTEx RNA data
library(genoppi)
gtex <- gtex_rna

# are DDG2P genes more often tissue-specific?
g2p <- g2p[grepl('brain', tolower(g2p$organ_specificity_list)),]
g2p <- g2p[g2p$mutation_consequence %in% 'loss of function']

# check with full background?
#fread('~/Projects/08_genesets/genesets/data/biomart/protein_coding_genes.tsv')

g2p$gene <- g2p$gene_symbol
g2p$significant <- TRUE
g2p_df <- g2p[,c('gene','significant')]
g2p_df <- g2p_df[!duplicated(g2p_df),]


g2p_tissue <- lapply_calc_hyper(g2p_df, gtex, intersectN = T, verbose = T)
gg_bar(g2p_tissue, bonf = 0.05/53)


# tissue-specifcity vs count
#tabl <- as.matrix(table(mrg$significant, mrg$organ_count))[,-1]
#pct <- tabl[2,] / colSums(tabl)
#plot(x=c(1:10,12), y=pct, main = 'Percent tissue-specific genes as a function of\n Organ Involvement count', xlab = 'Count') # no real trend
lapply_calc_hyper(organ)

# Aggregate data and test if having 1 gene in organ involment is likely to be tissue-specific
p.values <- lapply(1:10, function(i){
  organ$significant <- organ$organ_count %in% i
  organ_selected <- organ
  #organ_selected <- organ[organ$gene_symbol %in% g2p$gene_symbol,]
  statistics <- calc_hyper(organ_selected, gtex, intersectDf = data.frame(intersectN = F))
  statistics$statistics$pvalue
})

# plot the p-value for hypergeometric overlap test
names(p.values) <- 1:10
d <- stack(p.values)
colnames(d) <- c('pvalue', 'list_name')
d$pvalue <- as.numeric(d$pvalue)
gg_bar(d, bonf = 0.05 / 10)






