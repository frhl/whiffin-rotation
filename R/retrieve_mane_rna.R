

retrieve_mane_rna <- function(file = 'extdata/MANE.GRCh38.v0.95.select_ensembl_rna.fna'){
  
  # read data and convert to matrix
  x <- readLines(file)
  
  # define boundaries for identifiers
  x[grepl('>',x)] <- paste0(x[grepl('>',x)] , '<')
  x <- paste0(x, collapse = '')
  
  # all data is merged, now split them by our boundary elements
  sequences <- data.frame(unlist(strsplit(x, split = '>')))
  colnames(sequences) <- 'full'
  
  # And convert into 2xN data.frame with columns ID and sequence
  mat <- as.data.table(do.call(rbind, lapply(sequences$full, function(x) unlist(strsplit(x, split = '<')))))
  colnames(mat) <- c('id','seq')
  
  # format the header (extract the relevant columns)
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
  return(dt)
  
}