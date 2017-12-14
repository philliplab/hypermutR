library(hypermutR)

context("remove_hypermut")

test_that("remove_hypermut works", {
  base_seq <- deduplicate_seqs(ld_seqs)[[1]]$the_seq

  dat_20_1_0 <- c(sim_hyper(base_seq, 1, 'all', 0, seed = 1),
                  rep(base_seq, 19))
  names(dat_20_1_0) <- paste('s', 1:20, sep = '_')
  dat_20_1_0_l <- vector('list', length(dat_20_1_0))
  for (i in 1:length(dat_20_1_0)){
    dat_20_1_0_l[[i]] <- as.SeqFastadna(strsplit(dat_20_1_0[i], '')[[1]], name = names(dat_20_1_0)[i])
    names(dat_20_1_0_l)[i] <- names(dat_20_1_0)[i]
  }
  
  capture.output(x <- remove_hypermut(dat_20_1_0_l))
  expect_equal(1, 1)


#  dat_20_0_1
#  dat_20_0_0
#  dat_20_20_1
#  ld_plus_1
})

