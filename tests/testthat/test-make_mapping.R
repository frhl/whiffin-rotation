d <- fread('~/Projects/08_genesets/genesets/data/MANE/210710_MANE.GRCh38.v0.95.combined-table.txt')
d <- d[d$type == 'five_prime_UTR']

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
