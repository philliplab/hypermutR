#' Removes hypermutation from sequence data
#'
#' @param dat The sequence data
#' @export

remove_hypermut <- function(dat){
  cons_ints <- apply(consensusMatrix(dat), 2, function(x){which(x == max(x))[1]})
  cons <- paste(row.names(consensusMatrix(dat))[cons_ints], collapse = '')
  rm(cons_ints)

  results <- NULL
  hypermutants <- NULL
  
  ddat <- deduplicate_seqs(dat)
  for (i in 1:length(ddat)){
    result_scan <- scan_seq(cons, ddat[[i]]$the_seq, 'hyper')

    if (result_scan$p.value < 0.1){
      print(paste("Removing ", ddat[[i]]$dup_names, " because p value of ", result_scan$p.value, sep = ''))
      hypermutants <- rbind(hypermutants,
        data.frame(seq_name = ddat[[i]]$dup_names,
                   the_seq = ddat[[i]]$the_seq,
                   stringsAsFactors = FALSE)
        )
    } else {
      results <- rbind(results,
        data.frame(seq_name = ddat[[i]]$dup_names,
                   the_seq = ddat[[i]]$the_seq,
                   stringsAsFactors = FALSE)
      )
    }
  }

  seq_results <- DNAStringSet(results$the_seq)
  names(seq_results) <- results$seq_name

  seq_hypermutants <- DNAStringSet(hypermutants$the_seq)
  names(seq_hypermutants) <- hypermutants$seq_name

  return(
    list(seq_results = seq_results,
         seq_hypermutants = seq_hypermutants)
    )
  if (FALSE){
    # TEMP store for testing/debugging code
    # to be moved to unit tests sooner
    dat <- read.fasta('/home/phillipl/projects/hack_hypermut/KID141.fasta')
    dat <- read.fasta('/fridge/data/methods_paper/profiles/profile.fasta')
    remove_hypermut(ld_seqs)
    dat <- sim_hyper(ld_seqs, 20, 0.8, 0, verbose = TRUE)
    dat <- ld_seqs
    i <- 1
  }
}
