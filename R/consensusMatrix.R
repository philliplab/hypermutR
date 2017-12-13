#' Tallies characters at each position
#'
#' Designed to function like consensusMatrix from Biostrings.
#'
#' Note that seqinr uses lower case letters and Biostrings uses upper case letters.
#'
#' This function can handle upper and lower case input, but will only produce upper case output.
#'
#' To prevent a dependency on a bioconductor package, this is not properly set
#' up as a method for the consensusMatrix generic defined in Biostrings, but
#' rather uses a hacky name.
#'
#' @param x The sequence data of type SeqFastadna
#' @export

consensusMatrix_seqinr <- function(x, ...){
  max_length <- max(sapply(x, length))

  cons_matrix <- matrix(0, nrow = 18, ncol = max_length)
  row.names(cons_matrix) <- c("A", "C", "G", "T", "M", "R", "W", "S", "Y",
                              "K", "V", "H", "D", "B", "N", "-", "+", ".")
  
  seq_matrix <- matrix(NA_character_, nrow = length(x), ncol = max_length)
  for (i in 1:max_length){
    seq_matrix[,i] <- sapply(x, function(y, n){y[n]}, n = i)
    col_tab <- table(toupper(seq_matrix[,i]))
    cons_matrix[match(names(col_tab), row.names(cons_matrix)), i] <- col_tab
  }
  cons_matrix
}

