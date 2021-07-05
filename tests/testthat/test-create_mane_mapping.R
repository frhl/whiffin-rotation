


test_that("basic mapping", {
  
  d <- data.frame(bp = c('5-10;15-20'))
  res <- create_mane_mapping(d)
  expect_equal(as.numeric(res[1,]), c(5, 10, 5, 1, 6, 5))
  expect_equal(as.numeric(res[2,]), c(15, 20, 5, 7, 12, 5))
  
})
