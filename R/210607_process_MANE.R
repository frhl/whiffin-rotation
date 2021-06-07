# this is a messy workflow to create UTR and sequence files. Load them into a seperate script.


library(data.table)


# read data and convert to matrix
x <- readLines('extdata/MANE.GRCh38.v0.93.select_ensembl_rna.fna')
#x <- x[1:200]

# add charcaters for future splitting
x[grepl('>',x)] <- paste0(x[grepl('>',x)] , '<')
x <- paste0(x, collapse = '')
sequences <- data.frame(unlist(strsplit(x, split = '>')))
colnames(sequences) <- 'full'

# conver into 2xn data.frame
mat <- as.data.table(do.call(rbind, lapply(sequences$full, function(x) unlist(strsplit(x, split = '<')))))
colnames(mat) <- c('id','seq')

# format header
enstid_version <- unlist(lapply(strsplit(mat$id, split = ' '), function(x) x[1]))
enstid <- unlist(lapply(strsplit(enstid_version, split = '\\.'), function(x) x[1]))
ensgid_version <- unlist(lapply(strsplit(gsub('gene:','',mat$id), split = ' '), function(x) x[4]))
ensgid <- unlist(lapply(strsplit(ensgid_version, split = '\\.'), function(x) x[1]))
biotype <- unlist(lapply(strsplit(gsub('gene_biotype:','',mat$id), split = ' '), function(x) x[5]))
biotype_transcript <- unlist(lapply(strsplit(gsub('transcript_biotype:','',mat$id), split = ' '), function(x) x[6]))
gene_symbol <- unlist(lapply(strsplit(gsub('gene_symbol:','',mat$id), split = ' '), function(x) x[7]))
loc <- unlist(lapply(strsplit(gsub('chromosome:GRCh38:','',mat$id), split = ' '), function(x) x[3]))
loc_chr <- unlist(lapply(strsplit(loc, split = '\\:'), function(x) x[1]))
loc_start <- as.numeric(unlist(lapply(strsplit(loc, split = '\\:'), function(x) x[2])))
loc_end <- as.numeric(unlist(lapply(strsplit(loc, split = '\\:'), function(x) x[3])))
loc_direction <- as.numeric(unlist(lapply(strsplit(loc, split = '\\:'), function(x) x[4])))

# generate new matrix with these details
dt <- data.table(gene_symbol, ensgid, ensgid_version, enstid, enstid_version, biotype, biotype_transcript, chr = loc_chr, bp_start = loc_start, bp_end = loc_end, direction = loc_direction, seq = mat$seq) # ~ 17k genes
fwrite(dt, 'derived/210607_MANE.GRCh38.v0.93.select_ensembl_rna_matrix.txt', sep = '\t')
fwrite(dt, '~/Projects/08_genesets/genesets/data/MANE/210607_MANE.GRCh38.v0.93.select_ensembl_rna_matrix.txt', sep = '\t')

# find uAUG, uORFs
#manually search: uAUG = ATG, uORF = ATG + a stop codon (TAG, TAA, TGA) in frame

x <- read.table('extdata/MANE.GRCh38.v0.93.select_ensembl_genomic.gff', sep = '\t')
y <- x
table(y$V3)


cds <- y[ y$V3 %in% c('CDS','exon','start_codon' ,'stop_codon'),]
utrs <- y[ y$V3 %in% c('five_prime_UTR','three_prime_UTR'),]


## deal with CDS first
mat_cds <- do.call(rbind, strsplit(cds$V9,split =';'))

# standard
chr <- gsub('chr','',cds$V1)
origin <- cds$V2
type <- cds$V3
bp_start <- cds$V4
bp_end <- cds$V5
direction <- cds$V7
extra <- cds$V8

# different
id <-  gsub('ID=','',mat_cds[,1])
parent <- gsub('Parent=','',mat_cds[,2])
ensgid_version <-  gsub('gene_id=','',mat_cds[,3])
enstid_version <-  gsub('transcript_id=','',mat_cds[,4])
biotype <-  gsub('gene_type=','',mat_cds[,5])
biotype_transcript <-  gsub('transcript_type=','',mat_cds[,7])
gene_symbol <-  gsub('gene_name=','',mat_cds[,6])
exon <-  gsub('exon_number=','',mat_cds[,9])
exon_id <-  gsub('exon_id=','',mat_cds[,10])

# keep in CDS
cds_dt <- data.table(id, type, origin, chr, bp_start, bp_end, direction, extra, gene_symbol, parent, ensgid_version, enstid_version, biotype, biotype_transcript, exon, exon_id)
fwrite(cds_dt, '~/Projects/08_genesets/genesets/data/MANE/210607_MANE.GRCh38.v0.93.select_ensembl_genomic_matrix_cds.txt', sep = '\t')


## deal with CDS first
mat_utr <- do.call(rbind, strsplit(utrs$V9,split =';'))

# standard
chr <- gsub('chr','',utrs$V1)
origin <- utrs$V2
type <- utrs$V3
bp_start <- utrs$V4
bp_end <- utrs$V5
direction <- utrs$V7
extra <- utrs$V8

# different
id <-  gsub('ID=','',mat_utr[,1])
parent <- gsub('Parent=','',mat_utr[,2])
ensgid_version <-  gsub('gene_id=','',mat_utr[,3])
enstid_version <-  gsub('transcript_id=','',mat_utr[,4])
biotype <-  gsub('gene_type=','',mat_utr[,5])
biotype_transcript <-  gsub('transcript_type=','',mat_utr[,7])
gene_symbol <-  gsub('gene_name=','',mat_utr[,6])
exon <-  gsub('exon_number=','',mat_utr[,9])
exon_id <-  gsub('exon_id=','',mat_utr[,10])

# keep in CDS
utr_dt <- data.table(id, type, origin, chr, bp_start, bp_end, direction, extra, gene_symbol, parent, ensgid_version, enstid_version, biotype, biotype_transcript, exon, exon_id)
fwrite(utr_dt, '~/Projects/08_genesets/genesets/data/MANE/210607_MANE.GRCh38.v0.93.select_ensembl_genomic_matrix_UTRs.txt', sep = '\t')

final_dt <- rbind(utr_dt, cds_dt)
fwrite(final_dt, 'derived/210607_MANE.GRCh38.v0.93.select_ensembl_genomic_matrix.txt', sep = '\t')


fwrite(final_dt, '~/Projects/08_genesets/genesets/data/MANE/210607_MANE.GRCh38.v0.93.select_ensembl_genomic_matrix.txt', sep = '\t')




final_dt
