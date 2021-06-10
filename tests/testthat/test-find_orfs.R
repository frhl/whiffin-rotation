context('get_orf')

test_that("basic usage single ORF", {

  # basic usage 1)
  x = 'ATGTAG'
  result <- get_orf(x)
  expect_equal(names(result),'3_6')
  expect_equal(result[[1]], 'ATGTAG')
  
  # basic usage 2)
  x = 'ATGXXXTAG'
  result <- get_orf(x)
  expect_equal(names(result),'3_9')
  expect_equal(result[[1]], 'ATGXXXTAG')
  
  # basic usage 3)
  x = 'XATGXXXTAGX'
  result <- get_orf(x)
  expect_equal(names(result),'4_10')
  expect_equal(result[[1]], 'ATGXXXTAG')
  
})

test_that("basic usage multiple ORFs", {
  
  # Two start codons
  x = 'XATGATGTAG'
  result <- get_orf(x)
  expect_equal(names(result),c('4_10','7_10'))
  expect_equal(result[[1]], 'ATGATGTAG')
  expect_equal(result[[2]], 'ATGTAG')
  
  # multiple ORFs 
  x = 'ATGTAGATGTAG'
  result <- get_orf(x)
  expect_equal(names(result),c('3_6','9_12'))
  expect_equal(result[[1]], "ATGTAG")
  #expect_equal(result[[2]], "ATGTAGATGTAG")
  expect_equal(result[[2]], 'ATGTAG')
  
  # multiple ORFs (overlapping)
  x = 'ATGxATGxATGxTAGxTAGxTAG'
  result <- get_orf(x)
  expect_equal(names(result),c('3_15','7_19','11_23'))
  expect_equal(result[[1]], "ATGxATGxATGxTAG")
  expect_equal(result[[2]], "ATGxATGxTAGxTAG")
  expect_equal(result[[3]], "ATGxTAGxTAGxTAG")
  

})

test_that("ORF stops after first stop codon", {
  
  x = 'ATGXXXTAGXXXTAG'
  result <- get_orf(x)
  expect_equal(length(result), 1)
  
})

test_that("expected no ORF", {
  
  # only start
  x = 'ATGAAA'
  result <- get_orf(x)
  expect_equal(length(result), 0)
  
  # only end
  x = 'AAATAG'
  result <- get_orf(x)
  expect_equal(length(result), 0)
  
  # stop then start
  x = 'TAGATG'
  result <- get_orf(x)
  expect_equal(length(result), 0)
  
  # stop then start
  x = 'ATGXXTAG'
  result <- get_orf(x)
  expect_equal(length(result), 0)
  
  # two starts
  x = 'ATGATG'
  result <- get_orf(x)
  expect_equal(length(result), 0)
  
})




