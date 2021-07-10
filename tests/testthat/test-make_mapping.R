d <- fread('~/Projects/08_genesets/genesets/data/MANE/210710_MANE.GRCh38.v0.95.combined-table.txt')
d <- d[d$type == 'five_prime_UTR']
d <- d[d$bp != '-']

test_that("basic functionality", {
  
  
  row <- d[1, ]
  result <- make_mapping(row$seq, row$bp_cdna, row$bp)
  expect_equal(paste0(result$ref, collapse = ''), row$seq)
  
  row <- d[2, ]
  result <- make_mapping(row$seq, row$bp_cdna, row$bp)
  expect_equal(paste0(result$ref, collapse = ''), row$seq)
  
  row <- d[3, ]
  result <- make_mapping(row$seq, row$bp_cdna, row$bp)
  expect_equal(paste0(result$ref, collapse = ''), row$seq)
  
  row <- d[4, ]
  result <- make_mapping(row$seq, row$bp_cdna, row$bp)
  expect_equal(paste0(result$ref, collapse = ''), row$seq)
  
  row <- d[5, ]
  result <- make_mapping(row$seq, row$bp_cdna, row$bp)
  expect_equal(paste0(result$ref, collapse = ''), row$seq)
  
})


test_that("clinvar matches (forward strand)", {
  
  # load clinvar
  variants <- fread('extdata/clinvar/clinvar_20210626_chr1_22.txt')
  colnames(variants) <- c('chr','bp','ref','alt')
  variants <- variants[nchar(variants$ref) == 1]
  
  # first
  row <- d[d$enstid_version == 'ENST00000361915.8', ]
  chrom <- unique(row$chr)
  map <- make_mapping(row$seq, row$bp_cdna, row$bp)
  res  <- map[map$bp %in% variants$bp[variants$chr %in% chrom],]
  ref <- variants[variants$bp %in% map$bp & variants$chr %in% chrom]
  expect_equal(res$ref, ref$ref)
  expect_equal(paste0(map$ref, collapse = ''), row$seq)
  
  # ITGB4
  row <- d[d$enstid_version == 'ENST00000200181.8']
  chrom <- unique(row$chr)
  map <- make_mapping(row$seq, row$bp_cdna, row$bp)
  res  <- map[map$bp %in% variants$bp[variants$chr %in% chrom],]
  ref <- variants[variants$bp %in% map$bp & variants$chr %in% chrom]
  expect_equal(res$ref, ref$ref)
  expect_equal(paste0(map$ref, collapse = ''), row$seq)
  
  # RIN 3
  row <- d[d$enstid_version == 'ENST00000216487.12']
  chrom <- unique(row$chr)
  map <- make_mapping(row$seq, row$bp_cdna, row$bp)
  res  <- map[map$bp %in% variants$bp[variants$chr %in% chrom],]
  ref <- variants[variants$bp %in% map$bp & variants$chr %in% chrom]
  expect_equal(res$ref, ref$ref)
  expect_equal(paste0(map$ref, collapse = ''), row$seq)
  
  # RPN
  row <- d[d$enstid_version == 'ENST00000220676.2']
  chrom <- unique(row$chr)
  map <- make_mapping(row$seq, row$bp_cdna, row$bp)
  res  <- map[map$bp %in% variants$bp[variants$chr %in% chrom],]
  ref <- variants[variants$bp %in% map$bp & variants$chr %in% chrom]
  expect_equal(res$ref, ref$ref)
  expect_equal(paste0(map$ref, collapse = ''), row$seq)
  
  # EIF1B
  row <- d[d$enstid_version == 'ENST00000232905.4']
  chrom <- unique(row$chr)
  map <- make_mapping(row$seq, row$bp_cdna, row$bp)
  res  <- map[map$bp %in% variants$bp[variants$chr %in% chrom],]
  ref <- variants[variants$bp %in% map$bp & variants$chr %in% chrom]
  expect_equal(res$ref, ref$ref)
  expect_equal(paste0(map$ref, collapse = ''), row$seq)
    
})

test_that('clinvar matches (reverse strand)', {
  
  # FGF20 (reverse strand)
  row <- d[d$enstid_version == 'ENST00000180166.6']
  chrom <- unique(row$chr)
  map <- make_mapping(row$seq, row$bp_cdna, row$bp)
  res  <- map[map$bp %in% variants$bp[variants$chr %in% chrom],]
  ref <- variants[variants$bp %in% map$bp & variants$chr %in% chrom]
  expect_equal(res$ref, ref$ref)
  expect_equal(paste0(map$ref, collapse = ''), row$seq)
  
  # RPA3 (reverse strand)
  row <- d[d$enstid_version == 'ENST00000223129.8']
  chrom <- unique(row$chr)
  map <- make_mapping(row$seq, row$bp_cdna, row$bp)
  res  <- map[map$bp %in% variants$bp[variants$chr %in% chrom],]
  ref <- variants[variants$bp %in% map$bp & variants$chr %in% chrom]
  expect_equal(res$ref, ref$ref)
  expect_equal(paste0(map$ref, collapse = ''), row$seq)
  
  # SLC19A2 (reverse strand)
  row <- d[d$enstid_version == 'ENST00000236137.10']
  chrom <- unique(row$chr)
  map <- make_mapping(row$seq, row$bp_cdna, row$bp)
  res  <- map[map$bp %in% variants$bp[variants$chr %in% chrom],]
  ref <- variants[variants$bp %in% map$bp & variants$chr %in% chrom]
  #expect_equal(res$ref, ref$ref)
  #expect_equal(paste0(map$ref, collapse = ''), row$seq)
  
  utr.ens <- 'AGAGGCGTCTGTAGGGTAAAGCTGGGGGTTCTGTAGCCGGAGGCGGCGGCGAGTCCAGAACGTCCTGGCCTTACAGGGAGAAGGCGTCACTCGCGGTTACAAGTGCCTGACCCTCACTCCAGTTGGCGGAGGAGGAGAAGGAAGGGGCCGGGCCGGGTCCCCTCCCCTCGCGCCCCGG'
  
  # check first 100
  #comp <- lapply(1:1000, function(i){
  #for (i in 1:1000){
  #  row <- d[i, ]
  #  print(i)
  #  map <- make_mapping(row$seq, row$bp_cdna, row$bp)
  #  res  <- map[map$bp %in% variants$bp,]
  #  ref <- variants[variants$bp %in% map$bp]
  #  return(res$ref == ref$ref)
  #  
  #})
  
})


