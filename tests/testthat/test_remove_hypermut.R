library(hypermutR)

context("remove_hypermut")

test_that("remove_hypermut works", {
  base_seq <- deduplicate_seqs(ld_seqs)[[1]]$the_seq
  hyp_seq_one <- sim_hyper(base_seq, 1, 'all', 0, seed = 1)
  con_seq_one <- sim_hyper(base_seq, 1, 0, 'all', seed = 2)

  # 1 mut 19 nothing
  dat_20_1_0 <- c(hyp_seq_one,
                  rep(base_seq, 19))
  names(dat_20_1_0) <- paste("s", 1:20, sep = '_')
  dat_20_1_0_l <- make_list_SeqFastadna(dat_20_1_0)
  
  expect_output(remove_hypermut(dat_20_1_0_l), " because p value of ")
  x <- remove_hypermut(dat_20_1_0_l, verbose = FALSE)
  expect_equal(length(x$seq_results), 19)
  expect_equal(length(x$seq_hypermutants), 1)
  expect_equal(names(x$seq_hypermutants), 's_1')
  for (i in 1:length(x$seq_results)){
    expect_equal(x$seq_results[[i]], base_seq)
  }
  expect_equal(x$seq_hypermutants$s_1, as.character(hyp_seq_one))

  # 1 control 19 nothing
  dat_20_0_1 <- c(con_seq_one,
                  rep(base_seq, 19))
  names(dat_20_0_1) <- paste("c", 1:20, sep = '_')
  dat_20_0_1_l <- make_list_SeqFastadna(dat_20_0_1)

  x <- remove_hypermut(dat_20_0_1_l, verbose = FALSE)
  expect_equal(length(x$seq_results), 20)
  expect_equal(length(x$seq_hypermutants), 0)
  expect_null(x$seq_hypermutants)
  expect_equal(x$seq_results[[1]], as.character(con_seq_one))
  for (i in 2:length(x$seq_results)){
    expect_equal(x$seq_results[[i]], base_seq)
  }

  # 20 nothing
  dat_20_0_0 <- c(rep(base_seq, 20))
  names(dat_20_0_0) <- paste("n", 1:20, sep = '_')
  dat_20_0_0_l <- make_list_SeqFastadna(dat_20_0_0)

  x <- remove_hypermut(dat_20_0_0_l, verbose = FALSE)
  expect_equal(length(x$seq_results), 20)
  expect_equal(length(x$seq_hypermutants), 0)
  expect_null(x$seq_hypermutants)
  for (i in 1:length(x$seq_results)){
    expect_equal(x$seq_results[[i]], base_seq)
  }

  # 20 mut
  dat_20_20_0 <- c(rep(hyp_seq_one, 20))
  names(dat_20_20_0) <- paste("h", 1:20, sep = '_')
  dat_20_20_0_l <- make_list_SeqFastadna(dat_20_20_0)

  x <- remove_hypermut(dat_20_20_0_l, verbose = FALSE)
  expect_equal(length(x$seq_results), 20)
  expect_equal(length(x$seq_hypermutants), 0)
  expect_null(x$seq_hypermutants)
  for (i in 1:length(x$seq_results)){
    expect_equal(x$seq_results[[i]], as.character(hyp_seq_one))
  }

  # ld_seqs plus one hypermutant
  ld_seqs_test <- ld_seqs
  ld_seqs_test$hyp_seq_one <- as.SeqFastadna(strsplit(hyp_seq_one, '')[[1]], name = 'hyp_seq_one')

  expect_output(x <- remove_hypermut(ld_seqs_test), " because p value of ")
  expect_equal(length(x$seq_results), 872)
  expect_equal(length(x$seq_hypermutants), 1)
  expect_equal(names(x$seq_hypermutants), 'hyp_seq_one')
  expect_equal(x$seq_hypermutants$hyp_seq_one, as.character(hyp_seq_one))
})

test_that("remove_hypermut reports mut and control positions", {
  x <- remove_hypermut(hd_seqs)
  x <- hd_seqs
  y <- Biostrings::DNAStringSet(as.character(x))
  names(y) <- attr(x, "name")
  Biostrings::writeXStringSet(y, '/tmp/hd_seqs.fasta')
  z <- read.fasta('/tmp/hd_seqs.fasta')
  hd_seqs <- z
  devtools::use_data(hd_seqs, overwrite = T)
})

