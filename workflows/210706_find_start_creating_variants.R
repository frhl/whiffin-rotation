
# setup compute environment
library(parallel)
library(doParallel)
cores <- detectCores()
registerDoParallel(cores)

devtools::load_all()

# get data
d <- fread('~/Projects/08_genesets/genesets/data/MANE/210710_MANE.GRCh38.v0.95.combined-table.txt')
d <- d[d$type == 'five_prime_UTR']
d <- d[d$seq != '',]
#d <- head(d)

# enumerate codons that can be 
codon_target <- 'ATG'
pre_codon <- get_pre_codons(codon_target)

result_list <- (foreach (i=1:nrow(d)) %dopar% {
#result_list <- lapply(head(1:nrow(d)), function(i){
  
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
    mat_codons$grch38_bp = unlist(lapply(mat_codons$cdna_pos, function(pos){mapping$bp[mapping$cdna == pos]}))
  } 
  
  return(mat_codons)
  
  
  
})

result <- as.data.table(do.call(rbind, result_list))


#fwrite(result, 'extdata/210910_uATG_creating_variants.txt',sep = '\t')
result <- fread('extdata/210910_uATG_creating_variants.txt', sep = '\t')
#result <- fread('extdata/210910_uCGA_creating_variants.txt', sep = '\t')
result$ref <- result$before
result$alt <- result$after


variants <- fread('extdata/clinvar/clinvar_unpacked.txt')
colnames(variants) <- c('chr','pos','ref','alt', 'info')
variants$info <- unlist(lapply(variants$info, function(x){
  splitted <- unlist(strsplit(x, split = '\\;'))
  extract <- splitted[grepl('CLNSIG',splitted) & !grepl('CLNSIGINCL',splitted)]
  res <- gsub('CLNSIG\\=','',extract[1])
  res <- gsub(',_other','',res)
  res <- gsub(',_risk_factor','',res)
  res <- gsub(',_protective','',res)
  res <- gsub(',_drug_response','',res)
  res <- gsub(',_confers_sensitivity','',res)
  res <- gsub(',_Affects','',res)
  res <- gsub(',_association','',res)
  res <- gsub(',_confers_sensitivity','',res)
  return(res)
}))
variants$info <- factor(variants$info, levels = unique(variants$info))
colnames(variants)[2:3] <- 'grch38_bp'
#fwrite(variants, 'extdata/clinvar_unpacked_clnsig.txt',sep = '\t')

#variants <- fread('extdata/clinvar/clinvar_20210626_chr1_22.txt')
#colnames(variants) <- c('chr','pos','ref','alt')





#result <- fread('extdata/uAUG-creating_all_possible_annotated.txt')
#variants <- fread('extdata/clinvar/clinvar_20210626_chr1_22.txt')
#colnames(variants) <- c('chr','pos','ref','alt')
categories <- c('Uncertain_significance', 'Likely_benign','Benign/Likely_benign', 'Benign', 'Likely_pathogenic','Pathogenic')
mrg <- merge(result, variants, by = c('chr','grch38_bp','ref','alt'))
mrg <- mrg[mrg$info %in% categories]
mrg$info <- factor(mrg$info, levels = categories)

convert <- mrg[,c('codon_before','ref','alt')]
convert$change <- apply(convert[,c(2,3)], 1, paste0, collapse = '>')
convert$ref <- NULL
convert$alt <- NULL

d <- as.data.frame(table(mrg$info, mrg$codon_before))
p1 <- ggplot(d, aes(x=reorder(Var2, Freq), y=Var1, fill = log10(Freq), label = Freq)) +
  scale_fill_gradient(low = 'white',high ='red') +
  geom_tile() +
  geom_text() +
  theme_bw() +
  xlab('pre-ATG codon') +
  ylab('Clinvar significance') +
  geom_hline(yintercept = 1.5, linetype = 'dashed') +
  geom_hline(yintercept = 4.5, linetype = 'dashed')

p2 <- ggplot(d, aes(x=reorder(Var2, Freq), y=Freq)) +
  geom_bar(stat='identity') +
  theme_minimal()
  



mrg[mrg$ref.x == mrg$ref.y,]

#sort(table(mrg[mrg$before == mrg$ref,]$codon_before))
#x <- sort(table(mrg[mrg$before == mrg$ref,]$codon_before))

as.data.frame(x/sum(x))
