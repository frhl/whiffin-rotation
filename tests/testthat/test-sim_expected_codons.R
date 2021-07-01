context('simulate expected codons')

test_that("Basic cases", {
  
  #seq <- 'AAAAA'
  #res <- sim_expected_codons(seq, codons = c('XXX','AAA'), k = 2, iter = 1000)
  #expect_equal(res[[1]][1], 0)
  #expect_equal(res[[1]][2], 1)
  
  #seq <- c('AAAAAA','XXXXXX')
  #res <- sim_expected_codons(seq, codons = c('XXX','AAA'), k = 2, iter = 1000)
  #expect_equal(as.numeric(res[[1]]), c(0,1))
  #expect_equal(as.numeric(res[[2]]), c(1,0))
  
})


test_that("Simulation", {
  
  # 33 % chance of occurrence
  #seq <- 'ATGTAG'
  #res <- sim_expected_codons(seq, codons = c('ATG','TAG'), k = 2, iter = 1000)
  #expect_equal(round(res[[1]][1],1), round(res[[1]][2],1))
  
})
