context("get_oorf_kozak")

test_that("basic functionality", {
  
  # weak
  x1 = 'XXXATGXXX'
  x2 = 'ATGTAG'
  expect_equal(get_oorf_kozak(x1, x2)[[1]], 1)
  
  # moderate
  x1 = 'XXXATGGXX'
  x2 = 'ATGTAG'
  expect_equal(get_oorf_kozak(x1, x2)[[1]], 2)
  
  # One truncation and one cds alteration
  x1 <- 'ATGXXATGXXX'
  x2 <- 'ATGNNNTAG'
  expect_equal(get_oorf_kozak(x1, x2)[[1]], 1)
  expect_equal(get_oorf_kozak(x1, x2)[[2]], 2)
  
})

test_that("expect nothing to be returned", {
  
  # nothing
  x1 = 'ATGTAG'
  x2 = 'ATGTAG'
  expect_equal(length(get_oorf_kozak(x1, x2)), 0)
  
  # nothing
  x1 = 'AAAAAA'
  x2 = 'AAAAAA'
  expect_equal(length(get_oorf_kozak(x1, x2)), 0)

  
})


