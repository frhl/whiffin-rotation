context('get_altered_cds')

test_that("basic usage", {
  
  # truncation
  x1 <- 'XXXATGXX'
  x2 <- 'ATGNNNTAG'
  expect_equal(get_altered_cds(x1, x2)[[1]], "ATGXXATGNNNTAG")
  
  # elongated CDS
  x1 <- 'XXXATGXXX'
  x2 <- 'ATGNNNTAG'
  expect_equal(get_altered_cds(x1, x2)[[1]], "ATGXXXATGNNNTAG")
  
  # In frame and truncation
  x1 <- 'ATGXATGXXTAG'
  x2 <- 'ATGNNNTAG'
  expect_equal(get_altered_cds(x1, x2)[[1]], "ATGXXTAGATGNNNTAG")
  expect_equal(length(get_oorf(x1, x2)), 0)
  
  # One truncation and one cds alteration
  x1 <- 'ATGXXATGXXX'
  x2 <- 'ATGNNNTAG'
  res <- get_altered_cds(x1, x2)
  expect_equal(names(res), c("8_20", "3_0"))
  
  # Two truncations
  x1 <- 'ATGXXATGXXX'
  x2 <- 'ATGNNTAG'
  res <- get_altered_cds(x1, x2)
  expect_equal(names(res), c("3_0", "8_0"))

})


test_that("expext nothing to be returned", {
  
  # nothing
  x1 <- 'AAAAAAA'
  x2 <- 'TTTTTT'
  expect_equal(length(get_altered_cds(x1, x2)), 0)
  expect_true(is.list(get_altered_cds(x1, x2)))
  
  # ORFs in UTR and CDS but not overlapping
  x1 <- 'ATGTAG'
  x2 <- 'ATGTAG'
  expect_equal(length(get_altered_cds(x1, x2)), 0)
  expect_true(is.list(get_altered_cds(x1, x2)))
  
  # ORFs in UTR and CDS but not overlapping
  x1 <- 'ATGTAG'
  x2 <- 'ATGTAG'
  expect_equal(length(get_altered_cds(x1, x2)), 0)
  expect_true(is.list(get_altered_cds(x1, x2)))
  
})


