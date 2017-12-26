#' Internal function to resolve the ancestor argument of remove_hypermut
#'
#' @param ancestor Either 'consensus' to indicate that the consensus sequences must be computed, or 'first' to indicate that the first sequence in the dataset should be considered to be the ancestral sequence, or the ancestral sequence itself.
#' @param dat The sequence data

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


#' Removes hypermutation from sequence data
#'
#' @param dat The sequence data
#' @param verbose Prints name and p-value of removed/fixed sequences
#' @param fix_with Either false or a single letter. If not FALSE, then replace the hypermutated base with the letter indicated.
#' @param ancestor Either 'consensus' to indicate that the consensus sequences must be computed, or 'first' to indicate that the first sequence in the dataset should be considered to be the ancestral sequence, or the ancestral sequence itself.
#' @export

remove_hypermut <- function(dat, verbose = TRUE, fix_with = FALSE, ancestor = 'consensus'){
  cons <- ancestor_processing(ancestor, dat)

  results <- NULL
  hypermutants <- NULL
  all_mut_pos <- NULL
  
  ddat <- deduplicate_seqs(dat)
  for (i in 1:length(ddat)){
    result_scan <- scan_seq(cons, ddat[[i]]$the_seq, 'hyper', fix_with = fix_with)
    c_all_mut_pos <- NULL
    for (c_name in ddat[[i]]$dup_name){
      c_all_mut_pos <- rbind(c_all_mut_pos,
                             cbind(data.frame(seq_name = c_name,
                                              stringsAsFactors = FALSE),
                                   result_scan$all_mut_pos))
    }
    all_mut_pos <- rbind(all_mut_pos, c_all_mut_pos)
    if (result_scan$p.value < 0.1){
      if (verbose) {
        if (fix_with != FALSE) {
          print(paste("Fixing ", ddat[[i]]$dup_names, " because p value of ", result_scan$p.value, ". Replacing hypermutated positions with ", fix_with, ". Sequences will be left in the seq_results element of the return list and will also appear in the seq_hypermutants element.", sep = ''))
        } else {
          print(paste("Removing ", ddat[[i]]$dup_names, " because p value of ", result_scan$p.value, sep = ''))
        }
      }
      hypermutants <- rbind(hypermutants,
        data.frame(seq_name = ddat[[i]]$dup_names,
# fix_with                   the_seq = ddat[[i]]$the_seq,
                   the_seq = tolower(paste(result_scan$the_seq, collapse = '')),
                   stringsAsFactors = FALSE)
        )
    } else {
      results <- rbind(results,
        data.frame(seq_name = ddat[[i]]$dup_names,
# fix_with                   the_seq = ddat[[i]]$the_seq,
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
