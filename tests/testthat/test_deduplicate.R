library(hypermutR)

context("Deduplicate")

test_that("The deduplicate function works", {
  seqs <- c('aaa', 'aaa', 'aab', 'bbb', 'aaa', 'bbb')
  names(seqs) <- paste('s', 1:length(seqs), sep = '')
  dseqs <- deduplicate_seqs(seqs)
  expect_equal(length(dseqs), 3)
  expect_equal(paste(dseqs[[1]]$dup_names, collapse=''), 's1s2s5')
  expect_equal(paste(dseqs[[2]]$dup_names, collapse=''), 's3')
  expect_equal(paste(dseqs[[3]]$dup_names, collapse=''), 's4s6')
})
