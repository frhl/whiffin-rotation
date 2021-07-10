
# setup compute environment
library(parallel)
library(doParallel)
cores <- detectCores()
registerDoParallel(cores)

devtools::load_all()

# get data
d <- fread('~/Projects/08_genesets/genesets/data/MANE/210709_MANE.GRCh38.v0.95.combined-table.txt')
d <- d[d$type == 'five_prime_UTR']
#d <- head(d)

# enumerate codons that can be 
codon_target <- 'CGA'
pre_codon <- get_pre_codons(codon_target)

#result_list <- (foreach (i=1:nrow(d)) %dopar% {
result_list <- lapply(head(1:nrow(d)), function(i){
  
  row <- d[i,]
  sequence <- row$seq
  splitted_sequence <- split_seq(sequence)
  len <- nchar(sequence)
  
  mapping <- make_mapping(sequence, row$bp_cdna, row$bp)
  
  # go over every codon possible
  mat_codons <- do.call(rbind, lapply(1:nrow(pre_codon), function(j){
    
    cur_codon <- pre_codon[j,]
    positions <- find_codon(sequence, cur_codon$codon)
    
    # For all matching codons get positons
    mat <- as.data.frame(do.call(rbind, lapply(positions, function(p){
      position <- p
      range <- max(0, position-5):min(len, position+5)
      context <- paste0(splitted_sequence[range], collapse = '')
      
      data.table(
        gene_symbol = row$gene_symbol,
        enstid_version = row$ensgid_version,
        ensgid_version = row$ensgid_version,
        type = row$type,
        chr = row$chr,
        bp = row$bp,
        bp_cdna = row$bp_cdna,
        context = context, 
        position = position,
        codon_before = cur_codon$codon,
        codon_target = codon_target,
        before = cur_codon$before,
        after = cur_codon$after,
        cdna_pos = p - 2
      )

    })))
    

    return(mat)
    
  }))
  
  # map cdna postion to global position
  if (nrow(mat_codons) > 0){
    mat$grch38_bp <- unlist(lapply(mat_codons$cdna_pos, function(pos){mapping$bp[mapping$cdna == pos]}))
  }
  
  
  
  
})

result <- as.data.table(do.call(rbind, result_list))
result$type <- factor(result$type)
result$chr <- factor(result$chr)

head(result)

fwrite(result, 'extdata/210910_CGA_creating_variants.txt',sep = '\n')


