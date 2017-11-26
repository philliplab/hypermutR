library(hypermutR)

context("Data sets")

test_that("The supplied datasets still match the orignals", {
  expect_equal(length(ld_seqs), 872)
  expect_equal(length(hd_seqs), 691)

  expect_equal(nchar(paste(names(ld_seqs), collapse='')), 15692)
  expect_equal(nchar(paste(names(hd_seqs), collapse='')), 12013)

  expect_equal(names(ld_seqs)[1], "CGGTATTACATCACA_43")
  expect_equal(names(hd_seqs)[1], "ACCCGGGGAAGGCGA_5")
})
