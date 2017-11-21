#' Scans over sequences finding changes around given patterns
#'
#' @param cons The ancestral sequences to compare against
#' @param the_seq The query sequence
#' @param the_pattern The pattern to sequence for. Valid values: 'hyper' and 'control'
#' @export

scan_seq <- function(cons, the_seq, the_pattern){

  cons <- strsplit(toupper(cons), '')[[1]]
  the_seq <- strsplit(toupper(the_seq), '')[[1]]
  num.potential.mut <- 0
  num.mut <- 0
  num.potential.control <- 0
  num.control <- 0

  for( window.start.i in 1:(length(cons) - 2) ) {
      context.indx1 <- 1
      context.indx2 <- 2
      if( ( cons[window.start.i + 0 ] == "G" ) && # Reference must mutate from G
          ( the_seq[window.start.i + context.indx1 ] %in% c( "A", "G" ) ) && # Context position 1 must match R = [AG] in query
          ( the_seq[window.start.i + context.indx2 ] %in% c( "A", "G", "T" ) ) ){ # Context position 2 must match D = [AGT] in query
          num.potential.mut <- num.potential.mut + 1;
#          potential.pos <- rbind(potential.pos,
#            data.frame(seq.name = row.names(the_seq)[seq.i],
#                       pos = window.start.i,
#                       base.in.query = as.character( the_seq[ seq.i, window.start.i + 0 ] ),
#                       stringsAsFactors = F))
          if( ( as.character( the_seq[ window.start.i + 0 ] ) == "A" ) ) { # If G -> A mutation occurred
              num.mut <- num.mut + 1;
              # NOTE: rbind is VERY slow - check here for performance issues
#              if( fix.sequence ) {
#                  the_seq[ seq.i, window.start.i ] <<- as.DNAbin( fix.with );
#              }
          }
      }
      if( ( cons[ window.start.i + 0 ] == "G" ) && # Reference must mutate from G
          ( ( ( the_seq[ window.start.i + context.indx1 ] %in% c( "C", "T" ) ) && # Option 1 Context position 1 must match Y = [CT] in query
              ( the_seq[ window.start.i + context.indx2 ] %in% c( "A", "C", "G", "T" ) )) || # Option 1 Context position 2 must match N = [ACGT] in query
            ( ( the_seq[ window.start.i + context.indx1 ] %in% c( "A", "G" ) ) && # Option 2 Context position 1 must match R = [AG] in query
              ( the_seq[ window.start.i + context.indx2 ] ) == "C" ) ) ){ # Option 2 Context position 2 must match C in query
          num.potential.control <- num.potential.control + 1;
          if( as.character( the_seq[ window.start.i + 0 ] ) == "A" ) { # If G -> A mutation occureed
              #print( window.start.i );
              #print( as.character( the_seq[ seq.i, window.start.i + 0:2 ] ) );
              #print( as.character( cons[ 1, window.start.i + 0:2 ] ) );
              num.control <- num.control + 1;
          }
      }
  } # for window.start.i
#  p.value <- fisher.test( ( matrix( c( num.control, ( num.potential.control - num.control ), num.mut, ( num.potential.mut - num.mut ) ), nrow = 2, byrow = T ) ) )$p.value;
#  ## TODO: REMOVE
#  # print( row.names(the_seq)[seq.i])
#  # print( c( num.mut = num.mut, num.potential.mut = num.potential.mut, num.control = num.control, num.potential.control = num.potential.control, p.value = p.value ) );
#  row.names(potential.pos) <- NULL
#  return( list(p.value = p.value,
#               potential.pos = potential.pos) );
  return(list(num.mut = num.mut,
              num.potential.mut = num.potential.mut,
              num.control = num.control,
              num.potential.control = num.potential.control)
  )
  if (FALSE){
    # This is a TEMP store from some testing code - move to unit tests
    cons <- "GCTCCAGCTGGATTTGCGATTCTAAAGTGTAATAATAAGACATTCAATGGAACAGGACCATGCAATAATGTCAGCACAGTACAATGTACACATGGAATTAAGCCAGTGATATCAACTCAACTACTGTTAAATGGTAGCCTAGCAGAAGGGGAAATAATAATTAGATCTGAAAATATTACAGACAATGCTAAAACAATAATAGTACAATTTAATGAATCTATACTAATTAACTGTACAAGACCCGGCAATAATACAAGAAAAAGTGTAAGGATAGGAATAGGACGAGGACACACATTCTATGCAACAGGTGCTATAGGAAAAGATGCAAGACAAGCACATTGTAACGTTAGTACAACAGCATGGAACAAAACTTTACAAAGTGTAAAAGAGAAATTAAAAGAATACTTCAACAAAACAATAGAGTTCGCACCATCCTCAGGAGGGGATCTAGAAGTTACAACACATATTTTT"
    scan_seq(cons, as.character(ld_seqs[1]))
  }
}


      # if the window has any gaps in either sequence, skip it.
      # print(window.start.i)
      # Do not have to check position 0 for gaps since they will fail the == 'g' or == 'a' tests

      # Move the context forward if gaps are encountered to ensure that the
      # context pattern is matched to the sequence that the APOBEC enzyme
      # would have encountered.
      # Lots of IFs to prevent attempting to access values outside of the
      # valid range
      # I am sure that this can be done more elegantly...
#      while( as.character( the_seq[ seq.i, window.start.i + context.indx1 ] ) == "-" ){
#          context.indx1 <- context.indx1 + 1
#          context.indx2 <- context.indx2 + 1
#
#          if (window.start.i + context.indx2 > ncol(the_seq)){
#            break
#          }
#      }
#      if (window.start.i + context.indx2 > ncol(the_seq)){
#        next
#      }
#      while( as.character( the_seq[ seq.i, window.start.i + context.indx2 ] ) == "-" ){
#          context.indx2 <- context.indx2 + 1
#          if (window.start.i + context.indx2 > ncol(the_seq)){
#            break
#          }
#      }
#      if (window.start.i + context.indx2 > ncol(the_seq)){
#        next
#      }

