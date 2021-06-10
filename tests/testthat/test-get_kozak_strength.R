context('get_kozak_strength')

test_that("basic usage", {
  
  # weak
  x1 <- 'NNNATGNNNTAG'
  expect_equal(get_kozak_strength(x1)[[1]], 1)
  
  # moderate
  x2 <- 'ANNATGNNNTAG'
  expect_equal(get_kozak_strength(x2)[[1]], 2)
  x2 <- 'NNNATGGNNTAG'
  expect_equal(get_kozak_strength(x2)[[1]], 2)
  
  # strong
  x3 <- 'ANNATGGNNTAG'
  expect_equal(get_kozak_strength(x3)[[1]], 3)
  x3 <- 'GNNATGGNNTAG'
  expect_equal(get_kozak_strength(x3)[[1]], 3)
  
})

test_that("more than one kozak", {
  
  # moderate and weak
  x1 <- 'NNNATGGNNTAGNNNNATGNNNTAG'
  expect_equal(as.numeric(unlist(get_kozak_strength(x1))), c(2,1))
  
})

test_that("seperate between augs and just open reading frames", {
  
  x1 <- 'NNNATGNNNNTAG'
  expect_equal(get_kozak_strength(x1, only_orf = F)[[1]], 1)
  expect_equal(names(get_kozak_strength(x1, only_orf = F)), '6')
  expect_null(get_kozak_strength(x1, only_orf = T))

})

test_that("Kozak in start/end of sequence", {
  
  x1 <- 'ATGNNNNTAG'
  expect_equal(get_kozak_strength(x1, only_orf = F)[[1]], 1)

  x2 <- 'TAGATG'
  expect_equal(get_kozak_strength(x1, only_orf = F)[[1]], 1)
  
})

