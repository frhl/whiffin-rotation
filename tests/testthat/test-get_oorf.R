context('get_oorf')

test_that("basic usage",{
  
  # get out of frame overlapping (oORF)
  utr <- 'xxxATGxxxXXXxx'
  cds <- 'ATGXxxxTAGxxTAG'
  result1 <- get_oorf(utr,cds,inframe = F)
  expect_equal(extract_starts(result1), 6)
  result2 <- get_oorf(utr,cds,inframe = NULL) 
  expect_equal(extract_starts(result2), 6)
  
  # get N-terminal extension (in-frame)
  utr <- 'xxxATGxxxXXXxxY'
  cds <- 'ATGXxxxTAGxxTAG'
  result3 <- get_oorf(utr,cds,inframe = T, share_stops = F)
  expect_equal(extract_starts(result1), 6)
  
})


test_that("easy example", {

  # in frame
  x1 <- 'ATGAAA'
  x2 <- 'ATGAAATAG'
  expect_equal(get_oorf(x1, x2)[[1]], "ATGAAAATGAAATAG")
  
  # add one base to 5' end
  x1 <- 'XATGAAA'
  x2 <- 'ATGAAATAG'
  get_oorf(x1,x2)
  expect_equal(get_oorf(x1, x2)[[1]], "ATGAAAATGAAATAG")
  
  # only in seperate but not overlapping
  x1 <- 'ATGXXXTAG'
  x2 <- 'ATGYYYTAG'
  expect_equal(length(get_oorf(x1,x2)), 0)
  
})


test_that("in/out of frame with CDS", {
  
  # in frame
  x1 <- 'ATGAAA'
  x2 <- 'ATGAAATAG'
  expect_equal(get_oorf(x1, x2, inframe = T)[[1]], "ATGAAAATGAAATAG")
  expect_equal(length(get_oorf(x1,x2, inframe = F)), 0)
  
  # example 2
  x1 <- 'XXXATGAAA'
  x2 <- 'ATGAAATAG'
  expect_equal(get_oorf(x1, x2, inframe = T)[[1]], "ATGAAAATGAAATAG")
  expect_equal(length(get_oorf(x1,x2, inframe = F)), 0)
  
  # get out of cds frame oORF
  x1 <- 'ATG..'
  x2 <- 'ATG.TAG..TAG'
  expect_equal(get_oorf(x1, x2, inframe = F)[[1]], "ATG..ATG.TAG")
  expect_equal(length(get_oorf(x1,x2, inframe = T)), 0)
  
  # two out of frame
  x1 <- 'ATGxxxXATGxxXX'
  x2 <- 'ATGxxxTAGxxTAGxxTAG'
  expect_equal(extract_stops(get_oorf(x1, x2, inframe = F, share_stops = F)),c(28,33))
  expect_true(length(get_oorf(x1, x2, inframe = T, share_stops = F)) == 0)
  
  
})

test_that("N Terminal Extension (NTE)", {
  
  # basic example
  x1 <- 'xxxATGxxx'
  x2 <- 'ATGxxxTAG'
  expect_equal(extract_starts(get_oorf(x1,x2)), 6)
  
  # Two overlapping with shared stops
  x1 <- 'xxxATGxxxXXXGxxATGxxx'
  x2 <- 'ATGxxxTAG'
  expect_equal(extract_starts(get_oorf(x1,x2, share_stops = T)), c(6,18))
  expect_equal(extract_stops(get_oorf(x1,x2, share_stops = T)), c(30,30))
  
  # no shared stops (only strongest koza)
  expect_equal(extract_starts(get_oorf(x1,x2, share_stops = F)), c(6))
  expect_equal(extract_stops(get_oorf(x1,x2, share_stops = F)), c(30))  
  
  
  # two out of frame
  x1 <- 'ATGxxxXATGxxXX'
  x2 <- 'ATGxxxTAGxxTAGxxTAG'
  #expect_equal(extract_stops(get_oorf(x1, x2, inframe = F)),c(28,33))
  #expect_true(length(get_oorf(x1, x2, inframe = T)) == 0)
  
  
})

test_that("expexted to not find anything",{
  
  # not overlapping (overlapping)
  x1 <- 'ATGXX'
  x2 <- 'ATGYYYTAG'
  expect_equal(length(get_oorf(x1,x2)), 0)
  
  # in UTR and CDS but not overlapping
  x1 <- 'ATGXXXTAG'
  x2 <- 'ATGYYYTAG'
  expect_equal(length(get_oorf(x1,x2)), 0)
  
})


test_that("Different kozak contexts in 5' UTR",{
  
  # not overlapping (overlapping)
  x1 <- 'GxxATGGxxXXXxxxATGxxx'
  x2 <- 'ATGYYYTAG'
  get_kozak_strength(x1)
  
  # both codons returned
  expect_equal(extract_starts(get_oorf(x1,x2, share_stops = T)), c(6,18))

  # Non-shared stops returns only the one with strongest kozak (first codon)
  expect_equal(extract_starts(get_oorf(x1,x2, share_stops = F)), c(6))
  
})

