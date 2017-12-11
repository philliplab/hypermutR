library(hypermutR)

context("Data sets")

test_that("The supplied datasets still match the orignals", {
  expect_equal(length(ld_seqs), 872)
  expect_equal(length(hd_seqs), 691)

  expect_equal(nchar(paste(attr(ld_seqs, 'name'), collapse='')), 15692)
  expect_equal(nchar(paste(attr(hd_seqs, 'name'), collapse='')), 12013)

  expect_equal(attr(ld_seqs, 'name')[1], "CGGTATTACATCACA_43")
  expect_equal(attr(hd_seqs, 'name')[1], "ACCCGGGGAAGGCGA_5")
})
