library(hypermutR)

context("remove_hypermut")

test_that("remove_hypermut works", {
  base_seq <- deduplicate_seqs(ld_seqs)[[1]]$the_seq
  hyp_seq_one <- sim_hyper(base_seq, 1, 'all', 0, seed = 1)

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



#  dat_20_0_1
#  dat_20_0_0
#  dat_20_20_1
#  ld_plus_1
})

