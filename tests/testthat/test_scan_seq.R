library(hypermutR)

context("scan_seq")

test_that("scan_seq works", {

  cons_seq <- deduplicate_seqs(ld_seqs)[[1]]$the_seq
  
  # all muts mutatated
  query_seq <- sim_hyper(cons_seq, 1, 'all', 0, seed = 1)
  result <- scan_seq(cons_seq, query_seq)
  expect_equal(result$num.mut, 43)
  expect_equal(result$num_potential_mut, 43)
  expect_equal(result$num.control, 0)
  expect_equal(result$num.potential.control, 44)
  expect_equal(result$p_value, 1.50657928786482e-25)
  expect_equal(fisher.test(matrix(c(0, 44, 43, 0),
                                  nrow = 2,
                                  byrow = T),
                                  alternative = 'less')$p.value,
               result$p_value)

  # all control mutatated
  query_seq <- sim_hyper(cons_seq, 1, 0, 'all', seed = 2)
  result <- scan_seq(cons_seq, query_seq)
  expect_equal(result$num.mut, 0)
  expect_equal(result$num_potential_mut, 43)
  expect_equal(result$num.control, 44)
  expect_equal(result$num.potential.control, 44)
  expect_equal(result$p_value, 1)
  expect_equal(fisher.test(matrix(c(44, 0, 0, 43),
                                  nrow = 2,
                                  byrow = T),
                                  alternative = 'less')$p.value,
               result$p_value)

  # no mutations
  result <- scan_seq(cons_seq, cons_seq)
  expect_equal(result$num.mut, 0)
  expect_equal(result$num_potential_mut, 43)
  expect_equal(result$num.control, 0)
  expect_equal(result$num.potential.control, 44)
  expect_equal(result$p_value, 1)
  expect_equal(fisher.test(matrix(c(0, 44, 0, 43),
                                  nrow = 2,
                                  byrow = T),
                                  alternative = 'less')$p.value,
               result$p_value)

  # all control and mut mutated
  query_seq <- sim_hyper(cons_seq, 1, 'all', 'all', seed = 3)
  result <- scan_seq(cons_seq, query_seq)
  expect_equal(result$num.mut, 43)
  expect_equal(result$num_potential_mut, 43)
  expect_equal(result$num.control, 44)
  expect_equal(result$num.potential.control, 44)
  expect_equal(result$p_value, 1)
  expect_equal(fisher.test(matrix(c(44, 0, 0, 43),
                                  nrow = 2,
                                  byrow = T),
                                  alternative = 'less')$p.value,
               result$p_value)

  # 50% mut and 50% control
  query_seq <- sim_hyper(cons_seq, 1, .5, .5, seed = 3)
  result <- scan_seq(cons_seq, query_seq)
  expect_equal(result$num.mut, 21)
  expect_equal(result$num_potential_mut, 43)
  expect_equal(result$num.control, 22)
  expect_equal(result$num.potential.control, 44)
  expect_equal(result$p_value, 0.626499558301043)
  expect_equal(fisher.test(matrix(c(22, 22, 21, 22),
                                  nrow = 2,
                                  byrow = T),
                                  alternative = 'less')$p.value,
               result$p_value)

  # 50% mut and 0% control
  query_seq <- sim_hyper(cons_seq, 1, .5, 0, seed = 3)
  result <- scan_seq(cons_seq, query_seq)
  expect_equal(result$num.mut, 21)
  expect_equal(result$num_potential_mut, 43)
  expect_equal(result$num.control, 0)
  expect_equal(result$num.potential.control, 44)
  expect_equal(result$p_value, 1.38814127000831e-08)
  expect_equal(fisher.test(matrix(c(0, 44, 21, 22),
                                  nrow = 2,
                                  byrow = T),
                                  alternative = 'less')$p.value,
               result$p_value)

  # 0% mut and 50% control
  query_seq <- sim_hyper(cons_seq, 1, 0, .5, seed = 3)
  result <- scan_seq(cons_seq, query_seq)
  expect_equal(result$num.mut, 0)
  expect_equal(result$num_potential_mut, 43)
  expect_equal(result$num.control, 22)
  expect_equal(result$num.potential.control, 44)
  expect_equal(result$p_value, 1)
  expect_equal(fisher.test(matrix(c(22, 44, 0, 43),
                                  nrow = 2,
                                  byrow = T),
                                  alternative = 'less')$p.value,
               result$p_value)
})

test_that("scan_seq detects the correct positions", {
  cons_seq <-  "cccgggcccgat"
  query_seq <- "cccaggcccgat"
  result <- scan_seq(cons_seq, query_seq)
  expect_equal(nrow(result$all_mut_pos), 4)
  expect_equal(class(result$all_mut_pos), 'data.frame')
  expect_equal(sort(result$all_mut_pos$pos), c(4, 5, 6, 10))
  expect_equal(sort(result$all_mut_pos$base.in.query), c('A', 'G', 'G', 'G'))
  expect_equal(sort(result$all_mut_pos$full_seq), c('AGG', 'GAT', 'GCC', 'GGC'))
  expect_equal(sort(result$all_mut_pos$type), c('mut', 'mut', 'pot', 'pot'))
  expect_equal(sort(result$all_mut_pos$muted), c(F, F, F, T))
})

test_that("scan_seq handles gaps correctly", {
  the_seq <- "--GAA"
  result <- scan_seq(the_seq, the_seq)
  expect_equal(result$num_potential_mut + result$num.potential.control, 1)
  expect_equal(result$all_mut_pos$pos, 3)

  the_seq <- "CCG--AA"
  result <- scan_seq(the_seq, the_seq)
  expect_equal(result$num_potential_mut + result$num.potential.control, 1)
  expect_equal(result$all_mut_pos$pos, 3)

  the_seq <- "CGA-T"
  result <- scan_seq(the_seq, the_seq)
  expect_equal(result$num_potential_mut + result$num.potential.control, 1)
  expect_equal(result$all_mut_pos$pos, 2)

  the_seq <- "--GCA"
  result <- scan_seq(the_seq, the_seq)
  expect_equal(result$num_potential_mut + result$num.potential.control, 1)
  expect_equal(result$all_mut_pos$pos, 3)

  the_seq <- "CCG--AC"
  result <- scan_seq(the_seq, the_seq)
  expect_equal(result$num_potential_mut + result$num.potential.control, 1)
  expect_equal(result$all_mut_pos$pos, 3)

  the_seq <- "CGC-C"
  result <- scan_seq(the_seq, the_seq)
  expect_equal(result$num_potential_mut + result$num.potential.control, 1)
  expect_equal(result$all_mut_pos$pos, 2)
})

test_that("scan_seq fixes sequences correctly", {
  the_seq <- "CCAAA"
  the_con <- "CCGAA"
  expect_error(result <- scan_seq(the_con, the_seq, fix_with = 'hello'))
  result <- scan_seq(the_con, the_seq, fix_with = 'r')
  expect_true('the_seq' %in% names(result))
  expect_true(all(tolower(result$the_seq) %in% c(letters, '-')))
  expect_true('r' %in% result$the_seq)
})

test_that("scan_seq requires the consensus and the sequence to be of equal length", {
  the_seq <- "CCAA"
  the_con <- "CCGAA"
  expect_error(result <- scan_seq(the_con, the_seq), 'length\\(cons\\) == length\\(the_seq\\) is not TRUE')
})
