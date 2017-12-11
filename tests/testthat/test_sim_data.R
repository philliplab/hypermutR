library(hypermutR)

context("test_sim_data")

test_that("convert_alignment_to_matrix works", {
  expect_equal(1, 1)

  ld_mat <- convert_alignment_to_matrix(ld_seqs)
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

#sim_hyper <- function(dat, n1, n2, n3, seed = NULL, verbose = FALSE){

test_that("sim_hyper works", {
  test_seq <- ld_seqs[1]
  test_seq <- as.character(test_seq)

  n_g <- table(strsplit(test_seq, "")[[1]])['G']

  m10h <- sim_hyper(test_seq, 1, 10, 0, 1)
  g_count <- table(strsplit(as.character(m10h), "")[[1]])['g']
  g_gone <- n_g - g_count
  names(g_gone) <- NULL
  expect_equal(g_gone, 10)

  m10p <- sim_hyper(test_seq, 1, 0, 20, 1)
  g_count <- table(strsplit(as.character(m10p), "")[[1]])['g']
  g_gone <- n_g - g_count
  names(g_gone) <- NULL
  expect_equal(g_gone, 20)

  no_g <- sim_hyper(test_seq, 1, 'all', 'all', 1)
  no_g <- as.character(no_g)
  expect_false(grepl("g", no_g))
})

