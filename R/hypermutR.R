#' Removes hypermutation from sequence data
#'
#' @param dat The sequence data
#' @export

remove_hypermut <- function(dat){
  ddat <- deduplicate_seqs(dat)
  return(NULL)
}
