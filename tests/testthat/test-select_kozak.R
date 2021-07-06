
test_that("basic usage", {
  
  # strong kozak
  x <- 'CAAATGGxxxx'
  expect_equal(select_kozak(x), 6)
  
  # strong and weak
  x <- 'CAAATGGxxXXXxxxATGxxx'
  expect_equal(select_kozak(x), 6)
  
  # weak and moderate
  x <- 'ATGxxxXXXAxxATGxxx'
  expect_equal(select_kozak(x), 15)
  
})
