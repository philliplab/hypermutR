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
  }
  return(NULL)
  if (FALSE){
    # TEMP store for testing/debugging code
    # to be moved to unit tests sooner
    remove_hypermut(ld_seqs)
    dat <- ld_seqs
  }
}
