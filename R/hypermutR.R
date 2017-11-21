#' Removes hypermutation from sequence data
#'
#' @param dat The sequence data
#' @export

remove_hypermut <- function(dat){
  cons_ints <- apply(consensusMatrix(dat), 2, function(x){which(x == max(x))[1]})
  cons <- paste(row.names(consensusMatrix(dat))[cons_ints], collapse = '')
  rm(cons_ints)
  
  ddat <- deduplicate_seqs(dat)
  for (i in length(ddat)){
    # call scan_hyper
    result_hyper <- scan_seq(cons, ddat[[i]]$the_seq, 'hyper')
    # call scan_control
    # comp p.value
    # if significant
      # if fix, then fix and include
    # if not significant, then include
  }

  return(NULL)
  if (FALSE){
    # TEMP store for testing/debugging code
    # to be moved to unit tests sooner
    remove_hypermut(ld_seqs)
    dat <- ld_seqs
  }
}
