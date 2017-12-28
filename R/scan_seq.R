#' Scans over sequences finding changes around given patterns
#'
#' @param cons The ancestral sequences to compare against
#' @param the_seq The query sequence
#' @param the_pattern The pattern to sequence for. Valid values: 'hyper' and 'control'
#' @param fix_with Either false or a single letter. If not FALSE, then replace the hypermutated base with the letter indicated.
#' @export

scan_seq <- function(cons, the_seq, the_pattern, fix_with = FALSE){
  if (fix_with != FALSE){
    fix_with <- tolower(fix_with)
    stopifnot(tolower(fix_with) %in% letters)
  }
  cons <- strsplit(toupper(cons), '')[[1]]
  the_seq <- strsplit(toupper(the_seq), '')[[1]]
  stopifnot(length(cons) == length(the_seq))
  num_potential_mut <- 0
  num_mut <- 0
  num.potential.control <- 0
  num.control <- 0
  all_mut_pos <- NULL

  for( window.start.i in 1:(length(cons) - 2) ) {
      context.indx1 <- 1
      context.indx2 <- 2

      # if the window has any gaps in either sequence, skip it.
      # print(window.start.i)
      # Do not have to check position 0 for gaps since they will fail the == 'g' or == 'a' tests

      # Move the context forward if gaps are encountered to ensure that the
      # context pattern is matched to the sequence that the APOBEC enzyme
      # would have encountered.
      # Lots of IFs to prevent attempting to access values outside of the
      # valid range
      # I am sure that this can be done more elegantly...
      while( as.character( the_seq[ window.start.i + context.indx1 ] ) == "-" ){
          context.indx1 <- context.indx1 + 1
          context.indx2 <- context.indx2 + 1

          if (window.start.i + context.indx2 > length(cons)){
            break
          }
      }
      if (window.start.i + context.indx2 > length(cons)){
        next
      }
      while( as.character( the_seq[ window.start.i + context.indx2 ] ) == "-" ){
          context.indx2 <- context.indx2 + 1
          if (window.start.i + context.indx2 > length(cons)){
            break
          }
      }
      if (window.start.i + context.indx2 > length(cons)){
        next
      }

      # Check for hypermutated spots
      if( ( cons[window.start.i + 0 ] == "G" ) && 
            # Reference must mutate from G
          ( the_seq[window.start.i + context.indx1 ] %in% c( "A", "G" ) ) && 
            # Context position 1 must match R = [AG] in query
          ( the_seq[window.start.i + context.indx2 ] %in% c( "A", "G", "T" ) ) ){ 
            # Context position 2 must match D = [AGT] in query
          num_potential_mut <- num_potential_mut + 1;
          hyper_muted <- as.character( the_seq[ window.start.i + 0 ] ) == "A"
          all_mut_pos <- rbind(all_mut_pos,
            data.frame(pos = window.start.i,
                       base.in.query = as.character( the_seq[ window.start.i + 0 ] ),
                       full_seq = paste(as.character( the_seq[ window.start.i + c(0, context.indx1, context.indx2) ] ), 
                                        sep = '', collapse = ''),
                       type = 'mut',
                       muted = hyper_muted,
                       stringsAsFactors = F)
            )
          if( hyper_muted ) { 
              # If G -> A mutation occurred
              num_mut <- num_mut + 1;
              if ( fix_with != FALSE ) {
                the_seq[ window.start.i ] <- fix_with
              }
          }
      }
      # Check for control spots
      if( ( cons[ window.start.i + 0 ] == "G" ) && 
          # Reference must mutate from G
          ( ( ( the_seq[ window.start.i + context.indx1 ] %in% c( "C", "T" ) ) && 
                # Option 1 Context position 1 must match Y = [CT] in query
              ( the_seq[ window.start.i + context.indx2 ] %in% c( "A", "C", "G", "T" ) )) || 
                # Option 1 Context position 2 must match N = [ACGT] in query
            ( ( the_seq[ window.start.i + context.indx1 ] %in% c( "A", "G" ) ) && 
                # Option 2 Context position 1 must match R = [AG] in query
              ( the_seq[ window.start.i + context.indx2 ] ) == "C" ) ) ){ 
                # Option 2 Context position 2 must match C in query
          num.potential.control <- num.potential.control + 1;
          control_muted <- as.character( the_seq[ window.start.i + 0 ] ) == "A"
          all_mut_pos <- rbind(all_mut_pos,
            data.frame(pos = window.start.i,
                       base.in.query = as.character( the_seq[ window.start.i + 0 ] ),
                       full_seq = paste(as.character( the_seq[ window.start.i + c(0, context.indx1, context.indx2) ] ), 
                                        sep = '', collapse = ''),
                       type = 'pot',
                       muted = control_muted,
                       stringsAsFactors = F)
            )


          if( control_muted ) { 
              # If G -> A mutation occureed
              num.control <- num.control + 1;
          }
      }
  } # for window.start.i
  p_value <- fisher.test( matrix( c( num.control, ( num.potential.control - num.control ), num_mut, ( num_potential_mut - num_mut ) ), nrow = 2, byrow = T ), alternative = 'less' )$p.value;
  return(list(num_mut = num_mut,
              num_potential_mut = num_potential_mut,
              num.control = num.control,
              num.potential.control = num.potential.control,
              p_value = p_value,
              all_mut_pos = all_mut_pos,
              the_seq = the_seq)
  )
}
