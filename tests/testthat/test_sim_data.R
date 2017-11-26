library(hypermutR)

context("Simulation")

test_that("placeholder", {
  expect_equal(1, 1)

  ld_mat <- convert_alignment_to_matrix(ld_seqs)
  class(ld_mat)
  expect_equal(class(ld_mat), "matrix")
  expect_equal(nrow(ld_mat), 872)
  expect_equal(ncol(ld_mat), 471)
  row_samp <- c(656L, 167L, 775L, 544L, 278L, 508L, 326L, 441L, 265L, 586L,
                90L, 780L, 760L, 216L, 466L, 222L, 282L, 225L, 94L, 302L)
  col_samp <- c(263L, 246L, 454L, 383L, 288L, 293L, 193L, 116L, 229L, 240L,
                337L, 95L, 8L, 143L, 108L, 310L, 66L, 46L, 267L, 361L)
  let_samp <- mapply(function(row, col) {return(ld_mat[row, col])},
                     row = row_samp, col = col_samp)
  expect_equal(paste(let_samp, collapse = ''), "gcgtacacaacgccggtaat")
})
