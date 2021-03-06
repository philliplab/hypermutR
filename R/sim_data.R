#' Formats character vector as seqinr::read.fasta
#'
#' seqinr's read.fasta returns a list in each element is a SeqFastadna object.
#'
#' @param dat Named vector of character strings
#' @export

make_list_SeqFastadna <- function(dat){
  
  dat_l <- vector('list', length(dat))
  for (i in 1:length(dat)){
    dat_l[[i]] <- as.SeqFastadna(strsplit(dat[i], '')[[1]], 
                                        name = names(dat)[i])
    names(dat_l)[i] <- names(dat)[i]
  }
  dat_l
}

#' Converts an alignment to a matrix
#'
#' @param dat An alignment as produced by seqinr::read.fasta
#' @return A character matrix
#' @export

convert_alignment_to_matrix <- function(dat){
  if (class(dat) == "list"){
    stopifnot(all(sapply(dat, length) == max(sapply(dat, length))))
    max_length <- max(sapply(dat, length))
    mdat <- matrix(NA_character_, nrow = length(dat), ncol = max_length)
    for (i in 1:max_length){
      mdat[,i] <- sapply(dat, function(y, n){y[n]}, n = i)
    }
  } else {
    # There is a super weird format thing
    stopifnot(length(unique(nchar(dat))) == 1) # all sequences of equal length
    mdat <- matrix('-', nrow = length(dat), ncol = nchar(dat)[1])
    for (i in 1:length(dat)){
      mdat[i,] <- strsplit(tolower(as.character(dat[i])), '')[[1]]
    }
  }
  mdat
}

#' Simulates hypermutation in sequence datasets
#'
#' @param dat An alignment as produced by seqinr::read.fasta
#' @param n1 If less than one, the proportion of sequences to mutate, else
#' the number of sequences to mutate.
#' @param n2 If less than one, the proportion of potential hypermutation
#' positions to mutate, else the number of potential hypermutation positions to
#' mutate. The value 'all' is allowed.
#' @param n3 If less than one, the proportion of potential control
#' positions to mutate, else the number of potential control positions to
#' mutate. The value 'all' is allowed.
#' @param seed A seed value to use for the simulations.
#' @export

sim_hyper <- function(dat, n1, n2, n3, seed = NULL, verbose = FALSE){
  mdat <- convert_alignment_to_matrix(dat)
  if (n1 == 'all') {n1 <- length(dat)}
  if (n1 < 1) {n1 <- floor(length(dat)*n1)}
  if (n1 > length(dat)){
    warning('n1 is larger than the number of sequences, setting it to the number of sequences')
    n1 <- length(dat)
  }
  if (length(n1) > 1){
    seqs_to_mut <- min(n1, length(dat))
  } else {
    seqs_to_mut <- sample(1:length(dat), n1)
  }
  for (i in seqs_to_mut){
    hypermut_pos <- numeric(0)
    control_pos <- numeric(0)
    for (j in 1:(ncol(mdat)-2)){
      if (mdat[i,j]   ==   'g' &
          mdat[i,j+1] %in% c('a', 'g') & 
          mdat[i,j+2] %in% c('a', 'g', 't')){
        hypermut_pos <- c(hypermut_pos, j)
      }
      if (mdat[i,j]   ==   'g' &
          mdat[i,j+1] %in% c('c', 't')){
        control_pos <- c(control_pos, j)
      }
      if (mdat[i,j]   ==   'g' &
          mdat[i,j+1] %in% c('a', 'g') & 
          mdat[i,j+2] ==   'c'){
        control_pos <- c(control_pos, j)
      }
    }
    # now do the n2 and n3 mutations
    if (n2 == 'all') {n2 <- length(hypermut_pos)}
    if (n3 == 'all') {n3 <- length(control_pos)}
    if (n2 < 1) {n2 <- floor(length(hypermut_pos)*n2)}
    if (n3 < 1) {n3 <- floor(length(control_pos)*n3)}
    if (length(hypermut_pos) < n2) {stop("not enough spots to hypermutate")}
    if (length(control_pos) < n3) {stop("not enough control spots to hypermutate")}
#    hypermut_pos_dist <- c(hypermut_pos_dist, length(hypermut_pos))
#    control_pos_dist <- c(control_pos_dist, length(control_pos))

    names(dat)[i] <- paste('h', n2, 'c', n3, '_', 
                           names(dat)[i], sep = '')
    if (verbose) {print(names(dat)[i])}

    if (!is.null(seed)) {set.seed(seed)}
    n2_spots <- sample(hypermut_pos, n2)
    n3_spots <- sample(control_pos, n3)
    
    for (c_n2_spot in n2_spots){
      mdat[i, c_n2_spot] <- 'a'
    }
    for (c_n3_spot in n3_spots){
      mdat[i, c_n3_spot] <- 'a'
    }
  }
  new_dat <- as.SeqFastadna(apply(mdat, 1, paste, sep = '', collapse = ''),
                         name = attr(dat, 'name'))
  new_dat
}
