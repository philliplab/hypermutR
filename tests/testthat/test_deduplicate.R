library(hypermutR)

context("Deduplicate")

test_that("deduplicate_seqs works on basic character strings", {
  seqs <- c('aaa', 'aaa', 'aab', 'bbb', 'aaa', 'bbb')
  names(seqs) <- paste('s', 1:length(seqs), sep = '')
  dseqs <- deduplicate_seqs(seqs)
  expect_equal(length(dseqs), 3)
  expect_equal(paste(dseqs[[1]]$dup_names, collapse=''), 's1s2s5')
  expect_equal(paste(dseqs[[2]]$dup_names, collapse=''), 's3')
  expect_equal(paste(dseqs[[3]]$dup_names, collapse=''), 's4s6')
})

test_that("deduplicate_seqs works on the list based seqinr format", {
  dseqs <- deduplicate_seqs(ld_seqs)
  ld_seqs_char <- sapply(ld_seqs, function(x){ paste(x, collapse = '') })

  expect_equal(length(dseqs), 20)
  expect_equal(length(unlist(sapply(dseqs, `[`, 'dup_names'))), length(ld_seqs))
  for (i in 1:length(dseqs)){
    expect_true(dseqs[[i]]$the_seq %in% ld_seqs_char)
  }
})

