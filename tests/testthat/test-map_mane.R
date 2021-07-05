


test_that("basic mapping", {
  
  d <- data.frame(enstid = c('X'), bp = c('5-10;15-20'), seq = c(paste0(letters[c(5:10,15:20)], collapse = '')))
  expect_equal(map_mane(d, 5, adjust_pos = -2), 7)
  
})

test_that("Find codons in original sequence", {
  
  # generate fake sequence
  dna <- 'xxxxATGGTAxxxxxCGATAGCxxxxxTAGGAT'
  spacer <- strsplit(dna, split = '[^x]+')[[1]]
  rna <- strsplit(dna, split = 'x')[[1]]
  rna <- rna[rna != '']
  seq <- paste0(rna, collapse = '')
  
  # get lens
  spacer_lens <- unlist(lapply(spacer, nchar))
  rna_lens <- unlist(lapply(rna, nchar))
  
  # bp start
  sums_spacer_rna <- colSums(rbind(spacer_lens, rna_lens))
  position_matrix <- as.data.frame(cbind((cumsum(sums)-rna_lens),cumsum(sums)))
  position_matrix$V1 <- position_matrix$V1 + 1
  bps <- apply(position_matrix, 1, paste, collapse = '-')
  strsplit(dna, split = '')[[1]][5:10]
  strsplit(dna, split = '')[[1]][16:22]
  
  d <- data.frame(bp = paste0(bps, collapse = ';'), seq = seq)
  
  # first test ATG
  codon_atg <- find_codon(d$seq,'ATG')
  res_atg <- map_mane(d, codon_atg, adjust_pos = -2)
  expect_equal(res_atg, 5)
  
  # then test multiple codons
  nchar(d$seq)
  codon_at <- find_codon(d$seq,'AT')
  res_at <- map_mane(d, codon_at, adjust_pos = -2)
  expect_equal(strsplit(dna, split = '')[[1]][5:6], 'AT')
  expect_equal(strsplit(dna, split = '')[[1]][18:19], 'AT')
  expect_equal(strsplit(dna, split = '')[[1]][32:33], 'AT')
  
})


test_that("test with clinvar", {
  
  d <- fread('~/Projects/08_genesets/genesets/data/MANE/210705_MANE.GRCh38.v0.95.combined-table.txt')
  d <- d[d$gene_symbol %in% 'IRF6' &d$type == 'five_prime_UTR',]
  
  #d <- d[8,]
  
  codons_pos <- find_codon(d$seq, 'AT[^G]')
  map_mane(d, codons_pos, adjust_pos = -2)
  
  

  
  

})

