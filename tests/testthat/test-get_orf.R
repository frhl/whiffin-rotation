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
  result <- get_orf(x, share_stops = T)
  expect_equal(names(result),c('4_10','7_10'))
  expect_equal(result[[1]], 'ATGATGTAG')
  expect_equal(result[[2]], 'ATGTAG')
  
  # multiple ORFs 
  x = 'ATGTAGATGTAG'
  result <- get_orf(x, share_stops = T)
  expect_equal(names(result),c('3_6','9_12'))
  expect_equal(result[[1]], "ATGTAG")
  #expect_equal(result[[2]], "ATGTAGATGTAG")
  expect_equal(result[[2]], 'ATGTAG')
  
  # multiple ORFs (overlapping)
  x = 'ATGxATGxATGxTAGxTAGxTAG'
  result <- get_orf(x, share_stops = T)
  expect_equal(names(result),c('3_15','7_19','11_23'))
  expect_equal(result[[1]], "ATGxATGxATGxTAG")
  expect_equal(result[[2]], "ATGxATGxTAGxTAG")
  expect_equal(result[[3]], "ATGxTAGxTAGxTAG")
  
  # Multiple in-frame ORFs
  x = 'ATGxxxATGxxxTAG'
  result <- get_orf(x, share_stops = T)
  expect_equal(names(result),c('3_15','9_15'))
  

})

test_that("ORF stops after first stop codon", {
  
  x = 'ATGXXXTAGXXXTAG'
  result <- get_orf(x)
  expect_equal(length(result), 1)
  
})

test_that("shared stops (different kozak contexts)", {
  
  # two shared/overlapping ORFs
  x = 'ATGXXXATGXXXTAG'
  result <- get_orf(x, share_stops = T)
  expect_equal(length(result), 2)
  result <- get_orf(x, share_stops = F)
  expect_equal(length(result), 1)
  
  # four shared/overlapping ORFs
  x = 'ATGATGATGATGTAG'
  result <- get_orf(x, share_stops = T)
  expect_equal(length(result), 4)
  result <- get_orf(x, share_stops = F)
  expect_equal(length(result), 1)
  
  # One shared stuff weak and strong kozak
  x = 'AGGATGGXXxxxATGxxxTAG'
  expect_equal(as.numeric(unlist(get_kozak_strength(x))),c(3,1))
  result <- get_orf(x, share_stops = F)
  expect_equal(length(result), 1)
  result <- get_orf(x, share_stops = T)
  expect_equal(length(result), 2)
  
  # Two moderate kozaks result in the first being returned
  x = 'xxxATGGXXxxxATGGxxTAG'
  expect_equal(as.numeric(unlist(get_kozak_strength(x))),c(2,2))
  result <- get_orf(x, share_stops = F)
  expect_equal(length(result), 1)
  expect_equal(result[[1]], 'ATGGXXxxxATGGxxTAG')
  result <- get_orf(x, share_stops = T)
  expect_equal(length(result), 2)
  
})

test_that("expected no ORF", {
  
  # nothing
  x = 'AAAAAA'
  result <- get_orf(x)
  expect_equal(length(result), 0)
  
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




