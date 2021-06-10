context('find_truncating_augs')

test_that("basic usage", {
  
  x1 <- 'ATGXXATGXXXTGA'
  expect_equal(get_truncating_augs(x1)[[1]], "ATGXXATGXXXTGA")
  expect_equal(get_orf(x1)[[1]], "ATGXXXTGA")
  
  x1 <- 'XATGXXATGXXXTGA'
  expect_equal(get_truncating_augs(x1)[[1]], "ATGXXATGXXXTGA")
  expect_equal(get_orf(x1)[[1]], "ATGXXXTGA")
  
})
