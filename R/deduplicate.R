#' Deduplicates sequence data
#'
#' @param dat The sequence data (DNAStringSet)
#' @export

deduplicate_seqs <- function(dat){

  ndat <- as.character(dat)
  names(ndat) <- names(dat)
  dat <- ndat
  udat <- sort(unique(dat))
  
  ddat_lookup <- character(length(udat))
  ddat <- list()
  for (i in 1:length(udat)){
    cseq <- udat[i]
    names(cseq) <- NULL
    ddat[[i]] <- list(the_seq = cseq,
                      dup_names = character(0) )
    ddat_lookup[i] <- cseq
  }

  for (i in 1:length(dat)){
    matched_to <- which(dat[i] == ddat_lookup)
    ddat[[matched_to]]$dup_names <- c(ddat[[matched_to]]$dup_names, names(dat[i]))
  }

  return(ddat)
}
