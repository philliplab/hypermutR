library(hypermutR)

context("scan_seq")

test_that("scan_seq works", {

  cons_seq <- deduplicate_seqs(ld_seqs)[[1]]$the_seq
  
  # all muts mutatated
  query_seq <- sim_hyper(cons_seq, 1, 'all', 0, seed = 1)
  result <- scan_seq(cons_seq, query_seq)
  expect_equal(result$num.mut, 43)
  expect_equal(result$num.potential.mut, 43)
  expect_equal(result$num.control, 0)
  expect_equal(result$num.potential.control, 44)
  expect_equal(result$p.value, 1.50657928786482e-25)
  expect_equal(fisher.test(matrix(c(0, 44, 43, 0),
                                  nrow = 2,
                                  byrow = T),
                                  alternative = 'less')$p.value,
               result$p.value)

  # all control mutatated
  query_seq <- sim_hyper(cons_seq, 1, 0, 'all', seed = 2)
  result <- scan_seq(cons_seq, query_seq)
  expect_equal(result$num.mut, 0)
  expect_equal(result$num.potential.mut, 43)
  expect_equal(result$num.control, 44)
  expect_equal(result$num.potential.control, 44)
  expect_equal(result$p.value, 1)
  expect_equal(fisher.test(matrix(c(44, 0, 0, 43),
                                  nrow = 2,
                                  byrow = T),
                                  alternative = 'less')$p.value,
               result$p.value)

  # no mutations
  result <- scan_seq(cons_seq, cons_seq)
  expect_equal(result$num.mut, 0)
  expect_equal(result$num.potential.mut, 43)
  expect_equal(result$num.control, 0)
  expect_equal(result$num.potential.control, 44)
  expect_equal(result$p.value, 1)
  expect_equal(fisher.test(matrix(c(0, 44, 0, 43),
                                  nrow = 2,
                                  byrow = T),
                                  alternative = 'less')$p.value,
               result$p.value)

  # all control and mut mutated
  query_seq <- sim_hyper(cons_seq, 1, 'all', 'all', seed = 3)
  result <- scan_seq(cons_seq, query_seq)
  expect_equal(result$num.mut, 43)
  expect_equal(result$num.potential.mut, 43)
  expect_equal(result$num.control, 44)
  expect_equal(result$num.potential.control, 44)
  expect_equal(result$p.value, 1)
  expect_equal(fisher.test(matrix(c(44, 0, 0, 43),
                                  nrow = 2,
                                  byrow = T),
                                  alternative = 'less')$p.value,
               result$p.value)

  # 50% mut and 50% control
  query_seq <- sim_hyper(cons_seq, 1, .5, .5, seed = 3)
  result <- scan_seq(cons_seq, query_seq)
  expect_equal(result$num.mut, 21)
  expect_equal(result$num.potential.mut, 43)
  expect_equal(result$num.control, 22)
  expect_equal(result$num.potential.control, 44)
  expect_equal(result$p.value, 0.626499558301043)
  expect_equal(fisher.test(matrix(c(22, 22, 21, 22),
                                  nrow = 2,
                                  byrow = T),
                                  alternative = 'less')$p.value,
               result$p.value)

  # 50% mut and 0% control
  query_seq <- sim_hyper(cons_seq, 1, .5, 0, seed = 3)
  result <- scan_seq(cons_seq, query_seq)
  expect_equal(result$num.mut, 21)
  expect_equal(result$num.potential.mut, 43)
  expect_equal(result$num.control, 0)
  expect_equal(result$num.potential.control, 44)
  expect_equal(result$p.value, 1.38814127000831e-08)
  expect_equal(fisher.test(matrix(c(0, 44, 21, 22),
                                  nrow = 2,
                                  byrow = T),
                                  alternative = 'less')$p.value,
               result$p.value)

  # 0% mut and 50% control
  query_seq <- sim_hyper(cons_seq, 1, 0, .5, seed = 3)
  result <- scan_seq(cons_seq, query_seq)
  expect_equal(result$num.mut, 0)
  expect_equal(result$num.potential.mut, 43)
  expect_equal(result$num.control, 22)
  expect_equal(result$num.potential.control, 44)
  expect_equal(result$p.value, 1)
  expect_equal(fisher.test(matrix(c(22, 44, 0, 43),
                                  nrow = 2,
                                  byrow = T),
                                  alternative = 'less')$p.value,
               result$p.value)
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
  expect_equal(sort(result$all_mut_pos$hyper), c(F, F, F, T))
})

