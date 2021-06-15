context("get_orf_kozak")

test_that("basic functionality", {
  
  # weak
  x1 = 'XXXATGXXXXXXTAG'
  expect_equal(get_orf_kozak(x1)[[1]], 1)
  
  # moderate
  x1 = 'XXXATGGXXXXXTAG'
  expect_equal(get_orf_kozak(x1)[[1]], 2)
  
  # weak followed by strong 
  x = 'ATGxxATGGxxTAG'
  expect_equal(get_orf_kozak(x)[[1]], 3)
  
  # weak and strong in-frame
  x = 'ATGxxxAxxATGGxxTAG'
  expect_equal(get_orf_kozak(x)[[1]], 1)
  expect_equal(get_orf_kozak(x)[[2]], 3)
  
})

test_that("expect nothing to be returned", {
  
  # no stop codon
  x = 'ATGxxATGGxxXXX'
  expect_equal(length(get_orf_kozak(x)), 0)
  
  # only stop codon
  x = 'XXXxxXXXGxxTAG'
  expect_equal(length(get_orf_kozak(x)), 0)
  
})


