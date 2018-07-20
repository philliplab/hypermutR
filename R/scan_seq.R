#' Sequence scanner
#' 
#' @description Scans over sequences finding changes around given patterns
#'
#' @details The scan_seq function simultaneously passes two sliding windows
#' along the ancestral and query sequences. The sliding window is of length 3,
#' corresponding to the potentially hypermutated position and the 2 downstream
#' positions. At each position, the size of the window is increased until it
#' covers 3 non-gap characters in the query sequence. If a G is located at the
#' first position of the window, the position is considered a position of
#' interest and the query sequence is inspected to classify it as either a
#' hypermutation or control position, incrementing either the num_potential_mut
#' variable or the num_potential_control variable. The query sequence is
#' checked next and if the G mutated to an A, then the tally of the number of
#' possible hypermutations (num_mut) or the number of control mutations
#' (num_control) is incremented. 
#'
#' @param cons The ancestral sequences to compare against
#' @param the_seq The query sequence
#' @param the_pattern The pattern to sequence for. Valid values: 'hyper' and 'control'
#' @param fix_with Either false or a single letter. If not FALSE, then replace the hypermutated base with the letter indicated.
#'
#' @return The return value from scan_seq is a list that contains the number of
#' mutated hypermutation and control positions, the total number of potential
#' hypermutation and control positions, the p-value of the one-sided Fischer
#' exact test, the (possibly corrected) query sequence and the data.frame that
#' catalogs each individual position.
#'
#' @examples
#' scan_seq(paste(as.character(hd_seqs[1][[1]]), collapse = ''), paste(as.character(hd_seqs[2][[1]]), collapse = ''))
#'
#' @export

scan_seq <- function(cons, the_seq, the_pattern, fix_with = FALSE){
  if (fix_with != FALSE){
    fix_with <- tolower(fix_with)
    stopifnot(tolower(fix_with) %in% letters)
  }
  cons <- strsplit(toupper(cons), '')[[1]]
  the_seq <- strsplit(toupper(the_seq), '')[[1]]
  stopifnot(length(cons) == length(the_seq))
  if( length(cons) < 3){
    return(list(num_mut = 0,
                num_potential_mut = 0,
                num_control = 0,
                num_potential_control = 0,
                p_value = 1,
                all_mut_pos = NULL,
                the_seq = cons)
    )
  }
  num_potential_mut <- 0
  num_mut <- 0
  num_potential_control <- 0
  num_control <- 0
  all_mut_pos <- NULL

  for( window_start_i in 1:(length(cons) - 2) ) {
      context_indx1 <- 1
      context_indx2 <- 2

      # Move the context forward if gaps are encountered to ensure that the
      # context pattern is matched to the sequence that the APOBEC enzyme
      # would have encountered.
      while( as.character( the_seq[ window_start_i + context_indx1 ] ) == "-" ){
          context_indx1 <- context_indx1 + 1
          context_indx2 <- context_indx2 + 1

          if (window_start_i + context_indx2 > length(cons)){
            break
          }
      }
      if (window_start_i + context_indx2 > length(cons)){
        next
      }
      while( as.character( the_seq[ window_start_i + context_indx2 ] ) == "-" ){
          context_indx2 <- context_indx2 + 1
          if (window_start_i + context_indx2 > length(cons)){
            break
          }
      }
      if (window_start_i + context_indx2 > length(cons)){
        next
      }

      # Check for hypermutated spots
      if( ( cons[window_start_i + 0 ] == "G" ) && 
            # Reference must mutate from G
          ( the_seq[window_start_i + context_indx1 ] %in% c( "A", "G" ) ) && 
            # Context position 1 must match R = [AG] in query
          ( the_seq[window_start_i + context_indx2 ] %in% c( "A", "G", "T" ) ) ){ 
            # Context position 2 must match D = [AGT] in query
          num_potential_mut <- num_potential_mut + 1;
          hyper_muted <- as.character( the_seq[ window_start_i + 0 ] ) == "A"
          all_mut_pos <- rbind(all_mut_pos,
            data.frame(pos = window_start_i,
                       base_in_query = as.character( the_seq[ window_start_i + 0 ] ),
                       full_seq = paste(as.character( the_seq[ window_start_i + c(0, context_indx1, context_indx2) ] ), 
                                        sep = '', collapse = ''),
                       type = 'mut',
                       muted = hyper_muted,
                       stringsAsFactors = F)
            )
          if( hyper_muted ) { 
              # If G -> A mutation occurred
              num_mut <- num_mut + 1;
              if ( fix_with != FALSE ) {
                the_seq[ window_start_i ] <- fix_with
              }
          }
      }
      # Check for control spots
      if( ( cons[ window_start_i + 0 ] == "G" ) && 
          # Reference must mutate from G
          ( ( ( the_seq[ window_start_i + context_indx1 ] %in% c( "C", "T" ) ) && 
                # Option 1 Context position 1 must match Y = [CT] in query
              ( the_seq[ window_start_i + context_indx2 ] %in% c( "A", "C", "G", "T" ) )) || 
                # Option 1 Context position 2 must match N = [ACGT] in query
            ( ( the_seq[ window_start_i + context_indx1 ] %in% c( "A", "G" ) ) && 
                # Option 2 Context position 1 must match R = [AG] in query
              ( the_seq[ window_start_i + context_indx2 ] ) == "C" ) ) ){ 
                # Option 2 Context position 2 must match C in query
          num_potential_control <- num_potential_control + 1;
          control_muted <- as.character( the_seq[ window_start_i + 0 ] ) == "A"
          all_mut_pos <- rbind(all_mut_pos,
            data.frame(pos = window_start_i,
                       base_in_query = as.character( the_seq[ window_start_i + 0 ] ),
                       full_seq = paste(as.character( the_seq[ window_start_i + c(0, context_indx1, context_indx2) ] ), 
                                        sep = '', collapse = ''),
                       type = 'pot',
                       muted = control_muted,
                       stringsAsFactors = F)
            )


          if( control_muted ) { 
              # If G -> A mutation occureed
              num_control <- num_control + 1;
          }
      }
  } # for window_start_i
  p_value <- fisher.test( matrix( c( num_control, ( num_potential_control - num_control ), num_mut, ( num_potential_mut - num_mut ) ), nrow = 2, byrow = T ), alternative = 'less' )$p.value;
  return(list(num_mut = num_mut,
              num_potential_mut = num_potential_mut,
              num_control = num_control,
              num_potential_control = num_potential_control,
              p_value = p_value,
              all_mut_pos = all_mut_pos,
              the_seq = the_seq)
  )
}
