
devtools::load_all()


d1 <- fread('~/Projects/08_genesets/genesets/data/MANE/210705_MANE.GRCh38.v0.95.combined-table.txt')
d1 <- d1[d1$type == 'five_prime_UTR']
d2 <- fread('derived/tables/210708_MANE.v0.95.UTR_features.txt', sep = '\t')
mrg <- merge(d1, d2)
mrg1 <- mrg[mrg$u5_len > 50 & mrg$u5_AUG > 3 & mrg$enstid_version %in% plus,]
s <- mrg1[22,]

# bcftools query -f '%CHROM %POS %REF %ALT\n' clinvar_20210626.vcf.gz | awk '$1==16 && $2>2340573 && $2<2340728' 
s




# this matches!!
create_rna_seq <- function(x, bpin){
  len <- nchar(x$seq)
  data.frame(index = (len):1,
             bp = bpin:(bpin-len+1),
             ref = rev(unlist(strsplit(x$seq, split = '')))
             )
  
}

q <- create_rna_seq(s, 109621175)
q[q$bp == 109621025,]
q[q$bp == 109621064,]
q[q$bp == 109621101,]
q[q$bp == 109621109,]
q[q$bp == 109621124,]



# helpers
create_seq <- function(w){
  
  sequence <- seq(1,length(w)-4, by = 4)
  m <- as.data.frame(do.call(rbind, lapply(sequence, function(i)w[i:(i+3)])))
  
  colnames(m) <- c('chr','bp','ref','alt')
  m$bp <- as.numeric(m$bp)
  bpseq <- min(m$bp):max(m$bp)
  outseq <- unlist(lapply(bpseq, function(i){
    cur <- m$ref[m$bp == i]
    if (length(unique(cur)) > 1) stop('!!')
    if (length(cur) > 0) return(cur[1]) else return('.')
  }))
  
  #if (rev) return(paste0(rev(outseq), collapse = ''))
  return(paste0(outseq, collapse = ''))
}



#seq <- s$seq
seq <- s$seq
seqrev <- paste0(rev(unlist(strsplit(s$seq, split = ''))), collapse = '')


99782756 - 99782640
find_codon(s$seq, )



start <-  99782640 
99782756 - start + 1
unlist(strsplit(s$seq, split = ''))[117]

99782822 - start + 1
unlist(strsplit(s$seq, split = ''))[183]

99782805 - start + 1
unlist(strsplit(s$seq, split = ''))[166]

99782821 - start + 1
unlist(strsplit(s$seq, split = ''))[182]

99782822 - start + 1
unlist(strsplit(s$seq, split = ''))[183]






s_AMPD2 <- s
## stop this time???
109621064 - 109621175 - 1
unlist(strsplit(s_AMPD2$seq, split = ''))[112] # match

109621124 - 109621175 - 1
unlist(strsplit(s_AMPD2$seq, split = ''))[52] # match

109621025 - 109621175 - 1
unlist(strsplit(s_AMPD2$seq, split = ''))[151] # match

109621109 - 109621175 - 1
unlist(strsplit(s_AMPD2$seq, split = ''))[67] # match alternate

109621101 - 109621175 - 1
unlist(strsplit(s_AMPD2$seq, split = ''))[75] # match alternate

unlist(strsplit(s_AMPD2$seq, split = ''))[1945]

str_locate(s$seq, create_seq(w2))




99782802 - start
99782805 - start




# veryify lengths



extra = 427 - (80654728-80654490)

start <- 80654490
  

find_codon(s$seq, create_seq(w2))

find_codon(s$seq, create_seq(w3))

find_codon(s$seq, create_seq(w5))

unlist(strsplit(s$seq, split = ''))[80654974 - start + extra]

80654939 - start


create_seq()

(80654815 - 80654490) - extra


unlist(strsplit(s$seq, split = ''))[137]

regex <- 'AT..C.CG.C'
find_codon(s$seq, regex) # looks like the position is shifted
regex <- 'AT..C'
find_codon(s$seq, regex) # looks like the position is shifted


regex <- 'GGTTC.G'
find_codon(s$seq, regex) # looks like the position is shifted


extra <- 429 - (80654728-80654490) - nchar(regex)


22 == 2340595-2340573

30 == 2340603 - 2340573







