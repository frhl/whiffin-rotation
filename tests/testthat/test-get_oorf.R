
source('R/null_omit.R')
source('R/find_orfs.R')

test_that("easy example", {

  # in frame
  x1 <- 'ATGAAA'
  x2 <- 'ATGAAATAG'
  
  
  expect_equal(get_oorf(x1, x2)[[1]], "ATGAAAATGAAATAG")
  
  # add one base to 5' end
  x1 <- 'XATGAAA'
  x2 <- 'ATGAAATAG'
  get_oorf(x1,x2)
  expect_equal(get_oorf(x1, x2)[[1]], "ATGAAAATGAAATAG")
  
  # only in seperate but not overlapping
  x1 <- 'ATGXXXTAG'
  x2 <- 'ATGYYYTAG'
  expect_equal(length(get_oorf(x1,x2)), 0)
  
})


test_that("in/out of frame with CDS", {
  
  # in frame
  x1 <- 'ATGAAA'
  x2 <- 'ATGAAATAG'
  expect_equal(get_oorf(x1, x2, inframe = T)[[1]], "ATGAAAATGAAATAG")
  expect_equal(length(get_oorf(x1,x2, inframe = F)), 0)
  
  # example 2
  x1 <- 'XXXATGAAA'
  x2 <- 'ATGAAATAG'
  expect_equal(get_oorf(x1, x2, inframe = T)[[1]], "ATGAAAATGAAATAG")
  expect_equal(length(get_oorf(x1,x2, inframe = F)), 0)
  
  # get out of cds frame oORF
  x1 <- 'ATG..'
  x2 <- 'ATG.TAG..TAG'
  expect_equal(get_oorf(x1, x2, inframe = F)[[1]], "ATG..ATG.TAG")
  expect_equal(length(get_oorf(x1,x2, inframe = T)), 0)
  
})

test_that("",{
  
  # not overlapping (overlapping)
  x1 <- 'ATGXX'
  x2 <- 'ATGYYYTAG'
  get_oorf(x1,x2)
  
  
  
  
})

