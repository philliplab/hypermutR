#' @useDynLib hypermutR
#' @importFrom Rcpp sourceCpp
NULL

#' Resolve the ancestor argument of remove_hypermut
#'
#' @description
#' A simple function that resolves the ancestor argument passed to the
#' \code{remove_hypermut} function.
#'
#' @details
#' Three options exists for specifying the ancestral sequence to compare the
#' query sequences in the dataset to. 
#' \describe{
#'   \item{consensus}{If the value 'consensus' is specified via the ancestor
#'   parameter, a consensus sequence will be computed from the sequences in the
#'   input file. The letter that most frequently occurs is placed in the
#'   consensus sequence. In the case of ties, the first letter, when arranged
#'   alphabetically, is used.}
#'   \item{first}{Include the ancestral sequence as the first sequence in the
#'   input file and to set the value of the ancestor parameter to 'first'. In
#'   this case, the first sequence will be removed from the dataset before
#'   proceeding.}
#'   \item{Ancestral sequence}{The ancestor parameter can be assigned the
#'   ancestral sequence itself. The only validation that is performed on the
#'   last of the three options is to check that the sequence assigned to
#'   ancestor has the same length as the sequences in the input file.}
#' }
#' 
#' @param ancestor Either 'consensus' to indicate that the consensus sequences
#' must be computed, or 'first' to indicate that the first sequence in the
#' dataset should be considered to be the ancestral sequence, or the ancestral
#' sequence itself.
#' @param dat The sequence data
#'
#' @return A list with two elements:
#' \describe{
#'   \item{cons}{The ancestral sequence. Originally, only the consensus
#'   sequence option was supported, hence it is called 'cons' for historical
#'   reasons.}
#'   \item{dat}{The sequence dataset. The dataset that was passed into this
#'   function via the dat parameter with the first sequence removed in the case
#'   of the 'first' option and unaltered in case of the other two options.}
#' }
#'
#' @examples
#' ancestors <- ancestor_processing('consensus', ld_seqs)
#' print(ancestors$cons)
#' @export

ancestor_processing <- function(ancestor, dat){
  if (tolower(ancestor == 'consensus')){
    cons_ints <- apply(consensusMatrix_seqinr(dat), 2, function(x){which(x == max(x))[1]})
    cons <- paste(row.names(consensusMatrix_seqinr(dat))[cons_ints], collapse = '')
    rm(cons_ints)
  } else if (tolower(ancestor == 'first')){
    cons <- paste(as.character(dat[[1]]), collapse = '')
    n_dat <- vector('list', length(dat)-1)
    for (i in 1:length(n_dat)){
      n_dat[[i]] <- dat[[i+1]]
      names(n_dat)[i] <- names(dat)[i+1]
    }
    dat <- n_dat
  } else {
    cons <- ancestor
  }
  return(list('cons' = cons,
              'dat' = dat))
}
