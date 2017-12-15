library(hypermutR)

context("consensusMatrix")

test_that("consensusMatrix works", {
  cm <- consensusMatrix_seqinr(ld_seqs)
  the_rownames <- c("A", "C", "G", "T", "M", "R", "W", "S", "Y", "K", "V", "H",
                    "D", "B", "N", "-", "+", ".")
  expect_equal(attr(cm, "dimnames")[[1]], the_rownames)
  expect_equal(unique(apply(cm, 2, sum)), 872)
  expect_equal(apply(cm, 1, sum)['A'], structure(173525, names = 'A'))
  expect_equal(dim(cm), c(18, 471))

  expect_error(consensusMatrix_seqinr(ld_seqs[[1]]), ' == "list"')

  x <- list(ld_seqs[[1]])
  names(x) <- 'k'
  cm <- consensusMatrix_seqinr(x)
  expect_equal(attr(cm, "dimnames")[[1]], the_rownames)
  expect_equal(unique(apply(cm, 2, sum)), 1)
  expect_equal(apply(cm, 1, sum)['A'], structure(198, names = 'A'))
  expect_equal(dim(cm), c(18, 471))
})



