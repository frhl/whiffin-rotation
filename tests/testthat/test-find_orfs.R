context('find_orfs')

test_that("basic usage single ORF", {

  # basic usage 1)
  x = 'ATGTAG'
  result <- find_orfs(x)
  expect_equal(names(result),'3_6')
  expect_equal(result[[1]], 'ATGTAG')
  
  # basic usage 2)
  x = 'ATGXXXTAG'
  result <- find_orfs(x)
  expect_equal(names(result),'3_9')
  expect_equal(result[[1]], 'ATGXXXTAG')
  
  # basic usage 3)
  x = 'XATGXXXTAGX'
  result <- find_orfs(x)
  expect_equal(names(result),'4_10')
  expect_equal(result[[1]], 'ATGXXXTAG')
  
})

test_that("basic usage multiple ORFs", {
  
  # Two start codons
  x = 'XATGATGTAG'
  result <- find_orfs(x)
  expect_equal(names(result),c('4_10','7_10'))
  expect_equal(result[[1]], 'ATGATGTAG')
  expect_equal(result[[2]], 'ATGTAG')
  
  # multiple ORFs 
  x = 'ATGTAGATGTAG'
  result <- find_orfs(x)
  expect_equal(names(result),c('3_6','9_12'))
  expect_equal(result[[1]], "ATGTAG")
  #expect_equal(result[[2]], "ATGTAGATGTAG")
  expect_equal(result[[2]], 'ATGTAG')
  
  # multiple ORFs (overlapping)
  x = 'ATGxATGxATGxTAGxTAGxTAG'
  result <- find_orfs(x)
  expect_equal(names(result),c('3_15','7_19','11_23'))
  expect_equal(result[[1]], "ATGxATGxATGxTAG")
  expect_equal(result[[2]], "ATGxATGxTAGxTAG")
  expect_equal(result[[3]], "ATGxTAGxTAGxTAG")
  

})

test_that("ORF stops after first stop codon", {
  
  x = 'ATGXXXTAGXXXTAG'
  result <- find_orfs(x)
  expect_equal(length(result), 1)
  
})

test_that("expected no ORF", {
  
  # only start
  x = 'ATGAAA'
  result <- find_orfs(x)
  expect_equal(length(result), 0)
  
  # only end
  x = 'AAATAG'
  result <- find_orfs(x)
  expect_equal(length(result), 0)
  
  # stop then start
  x = 'TAGATG'
  result <- find_orfs(x)
  expect_equal(length(result), 0)
  
  # stop then start
  x = 'ATGXXTAG'
  result <- find_orfs(x)
  expect_equal(length(result), 0)
  
  # two starts
  x = 'ATGATG'
  result <- find_orfs(x)
  expect_equal(length(result), 0)
  
})




