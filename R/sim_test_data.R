#TODO: Remove:
library(Biostrings)
library(ape)

#' Removes ambiguous letters from a sequence by replacing them with a randomly
#' selected letter they represent
#' Move character swaps to strings in base R because it is faster
#' @param seq_dat The sequence in which to remove the ambiguity characters
#' @export
randomize_ambig <- function(seq_dat){
  if (class(seq_dat) == 'DNAString'){
    seq_dat <- DNAStringSet(seq_dat)
  }
  ambig_char <- IUPAC_CODE_MAP[!(names(IUPAC_CODE_MAP) %in% c('A', 'C', 'G', 'T'))]
  for (j in 1:length(seq_dat)){
    curr_seq <- strsplit(as.character(seq_dat[j]), '')[[1]]
    new_seq <- curr_seq
    for (i in 1:length(curr_seq)){
      curr_char <- curr_seq[i]
      if (!(curr_char %in% c('A', 'C', 'G', 'T', '-'))){
        n_options <- nchar(ambig_char[curr_char])
        sample_indx <- sample(1:n_options, 1)
        new_let <- substr(ambig_char[curr_char], sample_indx, sample_indx)
        new_seq[i] <- new_let
      }
    }
    seq_dat[j] <- DNAStringSet(paste(new_seq,sep='',collapse=''))
  }
  return(seq_dat)
}

simTestSeq <- function(n_unmut_hypermut, n_mut_hypermut,
                       n_unmut_control, n_mut_control,
                       shuffle = FALSE){
  unmut_hypermut <- '--GRD--'
  mut_hypermut <- '--ARD--'

  unmut_control1 <- '--GYN--'
  unmut_control2 <- '--GRC--'
  mut_control1 <- '--AYN--'
  mut_control2 <- '--ARC--'

  unmut_control <- unmut_control1
  mut_control <- mut_control1
  all_frags <- c(
    rep(unmut_hypermut, n_unmut_hypermut),
    rep(mut_hypermut, n_mut_hypermut),
    rep(unmut_control, n_unmut_control),
    rep(mut_control, n_mut_control))
  if (shuffle) {all_frags <- all_frags[sample(1:length(all_frags), length(all_frags))]}
  sim_seq <- DNAString(paste(all_frags, sep = '', collapse = ''))
  sim_seq <- randomize_ambig(sim_seq)
  names(sim_seq) <- paste('ref_seq', n_unmut_hypermut, 'uhm',
                          n_mut_hypermut, 'mhm',
                          n_unmut_control, 'ucm',
                          n_mut_control, 'mcm', sep = '_')
  return(sim_seq)
}

copyAndHyperMutate <- function(seq_dat, n_new_hypermut, n_new_controlmut){
  stopifnot(length(seq_dat) == 1)
  seq_dat <- as.character(seq_dat)
  if (n_new_hypermut > 0){
    for (i in 1:n_new_hypermut){
      seq_dat <- sub("--G([AG][AGT])", "--A\\1", seq_dat)
    }
  }
  if (n_new_controlmut > 0){
    for (i in 1:n_new_controlmut){
      seq_dat <- sub("--G(([CT][ACGT])|([AG]C))", "--A\\1", seq_dat)
    }
  }
  seq_dat <- DNAStringSet(seq_dat)
  names(seq_dat) <- paste('seq_+', n_new_hypermut, 'hm_+', n_new_controlmut, 'cm', sep = '')
  return(seq_dat)
}

# n_unmut_hypermut <- 10
# n_mut_hypermut <- 10
# n_unmut_control <- 10
# n_mut_control <- 10
# 
# ref_seq <- simTestSeq(10,10,10,10)
# test_file1 <- DNAStringSet(c(
# ref_seq,
# copyAndHyperMutate(ref_seq, 5, 5),
# copyAndHyperMutate(ref_seq, 1, 1),
# copyAndHyperMutate(ref_seq, 8, 8),
# copyAndHyperMutate(ref_seq, 10, 10),
# copyAndHyperMutate(ref_seq, 5, 0),
# copyAndHyperMutate(ref_seq, 1, 0),
# copyAndHyperMutate(ref_seq, 8, 0),
# copyAndHyperMutate(ref_seq, 10, 0)
# ))
# 
# 
# writeXStringSet(test_file1, '/tmp/test_file1.fasta', width=20000)
# 
# ref_seq <- simTestSeq(1,1,1,1)
# test_file_1111 <- DNAStringSet(c(
# ref_seq,
# copyAndHyperMutate(ref_seq, 0, 0),
# copyAndHyperMutate(ref_seq, 1, 0),
# copyAndHyperMutate(ref_seq, 0, 1),
# copyAndHyperMutate(ref_seq, 1, 1),
# copyAndHyperMutate(ref_seq, 2, 2)
# ))
# 
# writeXStringSet(test_file_1111, '/tmp/test_file_1111.fasta', width=20000)
# 
# 
# ref_seq <- simTestSeq(2,2,2,2)
# test_file_4x2 <- DNAStringSet(c(
# ref_seq,
# copyAndHyperMutate(ref_seq, 0, 0),
# copyAndHyperMutate(ref_seq, 1, 0),
# copyAndHyperMutate(ref_seq, 0, 1),
# copyAndHyperMutate(ref_seq, 1, 1),
# copyAndHyperMutate(ref_seq, 2, 0),
# copyAndHyperMutate(ref_seq, 2, 2)
# ))
# 
# writeXStringSet(test_file_4x2, '/tmp/test_file_4x2.fasta', width=20000)



