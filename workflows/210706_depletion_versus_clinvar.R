
library(data.table)
d_obs_ci <- fread('derived/210701_MANE.GRCh38.v0.95_codons_expt_ci.csv', sep = ',')
d_expt <- fread('derived/210701_MANE.GRCh38.v0.95_codons_expt.csv', sep = ',')
d_obs <- fread('derived/210701_MANE.GRCh38.v0.95_codons_obs.csv', sep = ',')

# merge observed / expected
mrg <- merge(d_expt, d_obs)
mrg$ensgid_version <- NULL
mrg$ensgid <- NULL
mrg

# observed versus expected
index_obs <- 3:66
index_expt <- 67:(67+63)
bool <- rowSums(obs)

obs <- mrg[,index_expt, with = F] 
expt <- mrg[,index_obs, with = F] 



# plot data
d <- as.data.frame(colSums(obs) / colSums(expt))
colnames(d) <- 'oe'
d$codon <- unlist(lapply(strsplit(rownames(d), split = '\\.'), function(x) x[2]))
# only codons that can be mutated into AUG
newcodons <- data.frame(codon = c('ATA','ATC','ATT', 'AAG','ACG','AGG', 'CTG','TTG','GTG'),
             mut = c('A>G','C>G','T>G', 'A>T','C>T','G>T', 'C>A','T>A','G>A'))




d <- merge(d, newcodons)
d$label <- paste0(d$codon, '\n(',d$mut,')')

depletion_order <- d$codon[order(d$oe)]
ggplot(d, aes(x=reorder(label, oe), y = oe)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  ggtitle('Observed versus Expected 5" UTR pre-AUG codons',
          '1000 simulations preservering di-nt frequency for each sequnce') +
  ylab('Observed / Expected') +
  xlab('Codon') +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1))

# are uAUG creating variants enriched?
features <- fread('derived/tables/210629_MANE.v0.95.UTR_features.txt', sep = '\t')
features <- features[,c('gene_symbol','ensgid_version','enstid_version')] # for mapping
utr <- fread('../../08_genesets/genesets/data/MANE/210705_MANE.GRCh38.v0.95.combined-table.txt')
aug <- fread('extdata/uAUG-creating_all_possible_annotated.txt')
aug$found <- !is.na(aug$gnomAD_AC) & aug$gnomAD_AC > 0


#aug <- aug[!is.na(aug$gnomAD_AC) & aug$gnomAD_AC != 0,]
#aug <- aug[aug$effect == 'uORF_created',]

# split clinvar codons to see what they were before
table_aug <- as.data.frame(table(aug$ref, aug$alt, aug$found))
colnames(table_aug) <- c('ref','alt','gnomad','freq')
bool <- duplicated(table_aug[,c('ref','alt')])
stopifnot(sum(colSums(table_aug[bool,] == table_aug[!bool,])) == 36)
table_aug <- cbind(table_aug[!bool,], data.frame(gnomad_freq = table_aug[bool,]$freq))
table_aug$ref_alt <- paste0(table_aug$ref,'>',table_aug$alt)
table_aug <- table_aug[table_aug$freq != 0,]
table_aug$fraction <- table_aug$gnomad_freq/table_aug$freq
ggplot(table_aug, aes(x=reorder(ref_alt, fraction), y=fraction)) +
  geom_bar(stat='identity', color = 'black', fill = 'firebrick3') + 
  ggtitle('Frequency of uAUG creating variants','uAUG-creating_all_possible_annotated.txt') + 
  xlab('Reference & Alternate allele') +
  ylab('Variant frequency')


res <- merge(d, table_aug, by.x = 'mut', by.y = 'ref_alt')
ggplot(res, aes(x=oe, y=fraction, label = label)) +
  geom_point(size=3) +
  #geom_bar(stat='identity', color = 'black', fill = 'firebrick3') + 
  ggtitle('5UTR pre-AUG variants versus O/E ratio','uAUG-creating_all_possible_annotated.txt') + 
  xlab('O/E') +
  ylab('Variant frequency') +
  ggrepel::geom_label_repel(size = 4, box.padding = 0.4)
