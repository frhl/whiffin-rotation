context('find_codon')

test_that("basic usage", {
  
  expect_equal(find_codon('ATG','ATG'), 3)
  expect_equal(find_codon('XATG','ATG'), 4)
  expect_equal(find_codon('ATGXXATG','ATG'), c(3,8))
  expect_null(find_codon('XXX','ATG'))
  expect_equal(find_codon('AXATG','A.ATG'), 3)
})
