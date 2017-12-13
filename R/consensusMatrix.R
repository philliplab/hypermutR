#' Tallies characters at each position
#'
#' Designed to function like consensusMatrix from Biostrings.
#' In seqinr, a method is added for the data type SeqFastadna.
#'
#' Note that seqinr uses lower case letters and Biostrings uses upper case letters.
#'
#' This method can handle upper and lower case input, but will only produce upper case output.
#'
#' @param x The sequence data of type SeqFastadna 
#' @export

setGeneric('consensusMatrix')

setMethod("consensusMatrix", signature(x = "SeqFastadna"),
  function(x, ...){
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
