
test_that("basic usage", {
  
  seq <- 'XXATGxxxXXXxxxXXXxxxTAGxxxXXXxxx'
  expect_equal(get_cds_proximity(seq), 9)
  expect_equal(get_leader_proximity(seq), 2)
  
  seq <- 'ATGxxxXXXxxxXXXxxxTAG'
  expect_equal(get_cds_proximity(seq), 0)
  expect_equal(get_leader_proximity(seq), 0)
  
})
