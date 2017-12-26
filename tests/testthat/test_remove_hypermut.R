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
})

#test_that("remove_hypermut works with real data (SLOW)", {
#  # ld_seqs plus one hypermutant
#  ld_seqs_test <- ld_seqs
#  ld_seqs_test$hyp_seq_one <- as.SeqFastadna(strsplit(hyp_seq_one, '')[[1]], name = 'hyp_seq_one')
#
#  expect_output(x <- remove_hypermut(ld_seqs_test), " because p value of ")
#  expect_equal(length(x$seq_results), 872)
#  expect_equal(length(x$seq_hypermutants), 1)
#  expect_equal(names(x$seq_hypermutants), 'hyp_seq_one')
#  expect_equal(x$seq_hypermutants$hyp_seq_one, as.character(hyp_seq_one))
#})

#test_that("remove_hypermut reports mut and control positions (SLOW)", {
#  x <- remove_hypermut(hd_seqs)
#  all_mut_pos <- x$all_mut_pos
#  names(x)
#
#  expect_equal(class(all_mut_pos), 'data.frame')
#  true_names <- c("base.in.query", "full_seq", "muted", "pos", "seq_name", "type")
#  expect_equal(sort(names(all_mut_pos)), true_names)
#  expect_equal(length(x$seq_results), length(unique(all_mut_pos$seq_name)))
#
#  randoms <- c(205L, 242L, 193L, 527L, 251L, 210L, 406L, 570L, 449L, 331L,
#               622L, 307L, 233L, 291L, 472L, 576L, 557L, 214L, 149L, 280L,
#               503L, 192L, 620L, 572L, 646L)
#  samp_indx <- randoms[4]
#  for (samp_indx in randoms){
#    expect_true(all_mut_pos[samp_indx,'seq_name'] %in% names(x$seq_result))
#  }
#})

test_that("remove_hypermut handles gaps correctly", {
  test_seqs <- list('s1' = as.SeqFastadna(c('c','t','g','-','c','c'), name = 's1'),
                    's2' = as.SeqFastadna(c('c','t','g','-','c','c'), name = 's2'),
                    's3' = as.SeqFastadna(c('c','t','g','-','c','c'), name = 's3'))
  x <- remove_hypermut(test_seqs)
  expect_equal(nrow(x$all_mut_pos), 3)
  expect_equal(unique(x$all_mut_pos$pos), 3)
})

test_that("remove_hypermut fixes correctly", {
  base_seq <- deduplicate_seqs(ld_seqs)[[1]]$the_seq
  hyp_seq_one <- sim_hyper(base_seq, 1, 'all', 0, seed = 1)
  con_seq_one <- sim_hyper(base_seq, 1, 0, 'all', seed = 2)

  # 1 mut 19 nothing
  dat_20_1_0 <- c(hyp_seq_one,
                  rep(base_seq, 19))
  names(dat_20_1_0) <- paste("s", 1:20, sep = '_')
  dat_20_1_0_l <- make_list_SeqFastadna(dat_20_1_0)
  
  expect_output(remove_hypermut(dat_20_1_0_l, fix_with = 'r'), " because p value of ")
  x <- remove_hypermut(dat_20_1_0_l, fix_with = 'r', verbose = FALSE)
  expect_equal(length(x$seq_results), 20)
  expect_equal(length(x$seq_hypermutants), 1)
  expect_equal(names(x$seq_hypermutants), 's_1')
  expect_equal(table(strsplit(x$seq_hypermutants$s_1, '')[[1]])['r'], structure(43, names = 'r'))

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

})

test_that('ancestral sequence specification is processed correctly', {
  hd_seqs_cons <- "GATCCAAGTGGACTGGGACTTCTAACCTGTAAAAAAAAGACAGTCGCTGGACCAGGCCAATGCAATACTCTCGGCATATTACAATGGGCCCGTGTGATGATTCCAACGATATCAACTCAACGAATGTTAAATCAGACCCGAGTAGAAGGAGATATAATAATTAGATGTGAAGAATTTCCTTTTAGTGCTAAAACATCAATAGTACAAGTTCATGTATCAATAGACATTGACTGTAGGATTCGGGGCAATCATACAAGAATAAGTGTCGGGATAGGAATGAGACTAGAACACACAGTCTATGCAACACTTGCTACACGAAAAGATGCGATACAAGCACCTTGTAACGTTACTTCGACAACTCTCTACACAACTTTTCACAGAGTAAAAGACGTGTTAAAGGAAGACTGCCATAAGACTGTCGAGTTCTTACCATCGCCTGGTAGTGATCTGGAATTGACAACGCATATGTTT"
  result <- ancestor_processing('consensus', hd_seqs)
  expect_equal(result$cons, hd_seqs_cons)
  expect_equal(result$dat, hd_seqs)

  result <- ancestor_processing('first', hd_seqs)
  expect_equal(result$cons, paste(as.character(hd_seqs[[1]]), collapse = ''))
  expect_equal(result$dat[[1]], hd_seqs[[2]])
  expect_equal(names(result$dat)[1], names(hd_seqs)[2])
  expect_equal(length(result$dat), length(hd_seqs)-1)

  result <- ancestor_processing('ACCGTG', ld_seqs)
  expect_equal(result$cons, 'ACCGTG')
  expect_equal(result$dat, ld_seqs)
})
