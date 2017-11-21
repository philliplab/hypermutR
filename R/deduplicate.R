#' Deduplicates sequence data
#'
#' @param dat The sequence data (DNAStringSet)
#' @export

deduplicate_seqs <- function(dat){

  udat <- unique(dat)
  
  ddat_lookup <- character(length(udat))
  ddat <- list()
  for (i in 1:length(udat)){
    cseq <- as.character(udat[i])
    names(cseq) <- NULL
    ddat[[i]] <- list( first_name = names(udat)[i],
                       the_seq = cseq,
                       dup_names = character(0) )
    ddat_lookup[i] <- cseq
  }

  for (i in 1:length(dat)){
    matched_to <- which(as.character(dat[i]) == ddat_lookup)
    ddat[[matched_to]]$dup_names <- c(ddat[[matched_to]]$dup_names, names(dat[i]))
  }

  for (i in 1:length(ddat)){
    stopifnot(ddat[[i]]$first_name %in% ddat[[i]]$dup_names)
  }

  return(ddat)

  if (FALSE){
    # This is a TEMP store from some testing code - move to unit tests
    dat <- ld_seqs
    deduplicate_seqs(ld_seqs)
    i <- 1
    x <- as.character(dat[i])
    names(x) <- NULL
    x
  }
}
