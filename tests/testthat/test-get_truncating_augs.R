

test_that("basic usage", {
  
  x <- 'ATGxxxXXATGxXXXxxxTAG'
  expect_equal(extract_starts(get_orf(x)),3)
  expect_equal(extract_starts(get_truncating_augs(x)),11)
  
})


test_that("don't expect anything", {
  
  x <- 'ATGxxxXXATGxTAGxxTAG'
  get_kozak_strength(x)
  expect_equal(extract_starts(get_orf(x, share_stops = T)), c(3,11))
  expect_true(length(get_truncating_augs(x)) == 0)
  
})
