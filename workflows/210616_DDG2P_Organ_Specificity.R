
# load DDG2P data
g2p <- fread('extdata/DDG2P_4_6_2021.csv')
g2p <- g2p[g2p$`DDD category` %in% c('confirmed','probable')]
colnames(g2p) <- gsub(' ','_',colnames(g2p))
dt <- g2p

# convert organ development boolean matrix into matrix stratified by gene.
organ <- strsplit(dt$organ_specificity_list, split = ';')
names(organ) <- dt$gene_symbol
organ <- stack(organ)
organ <- organ[!duplicated(organ),c(2,1)]
colnames(organ) <- c('gene_symbol','organ')
organ$involvement <- TRUE
organ <- data.table::dcast(setDT(organ), formula = gene_symbol ~ organ + involvement)
colnames(organ) <- gsub('_TRUE','',colnames(organ))
colnames(organ)[-1] <- paste0('organ ',colnames(organ)[-1])
colnames(organ) <- gsub(' ','_',colnames(organ))
organ[is.na(organ)] <- F

# check organ eye..
sum(organ$organ_Eye)
x <- organ[,grepl('Eye',colnames(organ)), with = F]
bool_eye <- apply(x, 1, any)

# organ eye
organ$organ_Eye <- NULL
organ$`organ_Eye:_Connective_tissue` <- NULL
organ$`organ_Eye:_Cornea` <- NULL
organ$`organ_Eye:_Lens` <- NULL
organ$`organ_Eye:_Optic_nerve` <- NULL
organ$`organ_Eye:_Physiology` <- NULL
organ$`organ_Eye:_Retina`<- NULL
organ$`organ_Eye:_Structure` <- NULL
organ$`organ_Eye:_Vasculature` <- NULL
organ$`organ_Eye:_Vitreous` <- NULL
organ$organ_Eye <- bool_eye

# remove cancer pre_disposition
organ$organ_Cancer_predisposition <- NULL

# some stats
hist(apply(organ[,-1], 1, sum))

fwrite(organ,'extdata/210616_DDG2P_organ_matrix.csv')



