library(hypermutR)

context("test_sim_data")

test_that("convert_alignment_to_matrix works", {
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
  test_seq <- ld_seqs[[1]]
  test_seq <- paste(as.character(test_seq), collapse = '')

  n_g <- table(strsplit(test_seq, "")[[1]])['g']

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

test_that("make_list_SeqFastadna works", {

  dduped <- deduplicate_seqs(ld_seqs)
  seqs <- rep("", 20)
  seqs[1] <- dduped[[1]]$the_seq
  seqs[2] <- dduped[[2]]$the_seq
  seqs[3:20] <- rep(dduped[[3]]$the_seq, length(3:20))
  names(seqs) <- paste("seq", 1:20, sep = "")

  x <- make_list_SeqFastadna(seqs[1])
  expect_equal(length(x), 1)  
  expect_equal(names(x), "seq1")
  expect_equal(length(x$seq1), 471)
  expect_equal(class(x), 'list')
  expect_equal(class(x[[1]]), 'SeqFastadna')
  expect_equal(paste(x[[1]][c(1,5,10,34,100,400)], collapse = ''), 'gcgaag')

  x <- make_list_SeqFastadna(seqs)
  expect_equal(length(x), 20)  
  expect_equal(names(x)[10], "seq10")
  expect_equal(length(x$seq15), 471)
  expect_equal(class(x), 'list')
  expect_equal(class(x[[19]]), 'SeqFastadna')
  expect_equal(paste(x[[1]][c(1,5,10,34,100,400)], collapse = ''), 'gcgaag')
})
