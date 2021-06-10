context('find_truncating_augs')

test_that("basic usage", {
  
  x1 <- 'ATGXXATGXXXTGA'
  expect_equal(find_truncating_augs(x1)[[1]], "ATGXXATGXXXTGA")
  expect_equal(find_orfs(x1)[[1]], "ATGXXXTGA")
  
})
