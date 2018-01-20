#' Deduplicates sequence data
#'
#' @description Deduplicates sequence data
#'
#' @details
#' The dataset is converted to a vector of character strings and the unique sequences are selected with the unique function. Looping over the unique sequences, a list is constructed in which each element corresponds to a unique sequence. Each element is also a list with the elements the_seq containing the actual sequences and dup_names, a vector of character strings listing the names of all sequences that matches the unique sequences stored in the_seq.
#' @param dat The sequence data (SeqFastadna)
#'
#' @return A list in which each unique sequence has an entry consisting of:
#' 1. the_seq containing the sequence as a character string.
#' 2. dup_names which is a character vector of the headers of all the
#'    sequences that had that exact sequence.
#'
#' @examples
#' deduplicate_seqs(ld_seqs)
#' deduplicate_seqs(c('aaa', 'aaa', 'aab'))
#'
#' @export

deduplicate_seqs <- function(dat){

  if (class(dat) == 'list'){
    tdat <- character(length(dat))
    for (i in 1:length(dat)){
      tdat[i] <- paste(dat[[i]], collapse = '')
      names(tdat)[i] <- attr(dat[[i]], 'name')
    }
    dat <- tdat
  }
  ndat <- as.character(dat)
  names(ndat) <- names(dat)
  if (is.null(names(ndat))){stop('Vector must be named, use unique function instead')}
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
