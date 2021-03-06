#' @title Detect and process hypermutated sequences
#'
#' @description 
#' This function is a wrapper that will detect and either remove or 'fix'
#' hypermutated sequences depending of the value of the 'fix_with' argument. 
#'
#' @details
#' It calls \code{ancestor_processing} to obtain the ancestral sequence to
#' compare the query sequences to, then calls \code{deduplicate_seqs} to remove
#' duplicate sequences for performance reasons, next loops over each unique
#' sequence, comparing it to the ancestral sequence with \code{scan_seq} and
#' finally collates the results.
#'
#' @seealso \code{\link{ancestor_processing}}, \code{\link{deduplicate_seqs}},
#' and \code{\link{scan_seq}}
#'
#' @param dat The sequence data. The structure must match the format produced
#' by read.fasta from the seqinr package. This is a list in which each element
#' represents a single sequence. Each element is of class \code{SeqFastadna}
#' and consists of a vector of single letters of class character with optional
#' attributes names and Annot. 
#' @param verbose If TRUE, print the name and p-value of removed/fixed
#' sequences.
#' @param fix_with Either FALSE or a single letter. If not FALSE, then replace
#' the hypermutated base with the letter indicated.
#' @param ancestor Either 'consensus' to indicate that the consensus sequences
#' must be computed, or 'first' to indicate that the first sequence in the
#' dataset should be considered to be the ancestral sequence, or the ancestral
#' sequence itself.
#' @param p_value The p-value used by the one-sided fischer test.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{all_mut_pos}{A data.frame that contains all the positions in all the
#'   sequences that are either a hypermutation or control position.}
#'   \item{seq_results}{A list that stores all the sequences that did not
#'   contain any hypermutation and, in the case that the fix_with parameter was
#'   set, those sequences with hypermutation that was corrected.}
#'   \item{seq_hypermutants}{A list of all the sequences that contains hypermutation.}
#' }
#'
#' @examples
#' result <- remove_hypermut(hd_seqs)
#' names(result)
#' str(result)
#' @export

remove_hypermut <- function(dat, verbose = TRUE, fix_with = FALSE, ancestor = 'consensus', p_value = 0.05){
  cons <- ancestor_processing(ancestor, dat)
  dat <- cons$dat
  cons <- cons$cons

  results <- NULL
  hypermutants <- NULL
  all_mut_pos <- NULL
  
  ddat <- deduplicate_seqs(dat)
  for (i in 1:length(ddat)){
    result_scan <- rcpp_scan_seq(cons, ddat[[i]]$the_seq, fix_with = fix_with)
    c_all_mut_pos <- NULL
    for (c_name in ddat[[i]]$dup_name){
      c_all_mut_pos <- rbind(c_all_mut_pos,
                             cbind(data.frame(seq_name = c_name,
                                              stringsAsFactors = FALSE),
                                   result_scan$all_mut_pos))
    }
    all_mut_pos <- rbind(all_mut_pos, c_all_mut_pos)
    if (result_scan$p_value < p_value){
      if (verbose) {
        if (fix_with != FALSE) {
          print(paste("Fixing ", ddat[[i]]$dup_names, " because p value of ", result_scan$p_value, 
                      ". Replacing hypermutated positions with ", fix_with, 
                      ". Sequences will be left in the seq_results element of the return list and will", 
                      " also appear in the seq_hypermutants element.", sep = ''))
        } else {
          print(paste("Removing ", ddat[[i]]$dup_names, " because p value of ", result_scan$p_value, sep = ''))
        }
      }
      hypermutants <- rbind(hypermutants,
        data.frame(seq_name = ddat[[i]]$dup_names,
                   the_seq = tolower(paste(result_scan$the_seq, collapse = '')),
                   stringsAsFactors = FALSE)
        )
    } else {
      results <- rbind(results,
        data.frame(seq_name = ddat[[i]]$dup_names,
                   the_seq = tolower(paste(result_scan$the_seq, collapse = '')),
                   stringsAsFactors = FALSE)
      )
    }
  }
  if (fix_with != FALSE){
    results <- rbind(results, hypermutants)
  }

  results_list <- vector('list', nrow(results))
  for (i in 1:nrow(results)){
    results_list[[i]] <- results[i,'the_seq']
    names(results_list)[i] <- results[i, 'seq_name']
  }

  if (!is.null(hypermutants)){
    hypermutants_list <- vector('list', nrow(hypermutants))
    for (i in 1:nrow(hypermutants)){
      hypermutants_list[[i]] <- hypermutants[i, 'the_seq']
      names(hypermutants_list)[i] <- hypermutants[i, 'seq_name']
    }
  } else {
    hypermutants_list <- NULL
  }

  return(
    list(seq_results = results_list,
         seq_hypermutants = hypermutants_list,
         all_mut_pos = all_mut_pos)
    )
}
