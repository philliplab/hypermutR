
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

Rcpp::List rcpp_scan_seq(CharacterVector r_cons, CharacterVector r_the_seq, 
  CharacterVector r_the_pattern, CharacterVector r_fix_with){

  std::cout << r_fix_with << std::endl;

  int num_potential_mut = 0;
  int num_mut = 0;
  int num_potential_control = 0;
  int num_control = 0;
  std::vector< int > all_mut_pos;
  bool hypermuted = false;
  bool fix = false;
  std::string fix_with_str;
  char fix_with;
  if (r_fix_with[0] != "FALSE" && r_fix_with[0] != "F"){
    fix_with_str = r_fix_with[0];
    if (fix_with_str.size() != 1){
      throw std::exception();
    } 
    fix_with = fix_with_str[0];
    fix = true;
  }
  
  std::string cons;
  for (int i = 0; i < r_cons[0].size(); ++i){
    cons += r_cons[0][i];
  }
  std::string the_seq;
  for (int i = 0; i < r_the_seq[0].size(); ++i){
    the_seq += r_the_seq[0][i];
  }

  int context_indx1;
  int context_indx2;

//  std::cout << r_cons.size() << std::endl;
//  std::cout << r_cons << std::endl;
//  std::cout << cons.size() << std::endl;
//  std::cout << cons << std::endl;
//  std::cout << the_seq.size() << std::endl;
//  std::cout << the_seq << std::endl;

  for (int window_start_i = 0; window_start_i < cons.size(); ++window_start_i){
    context_indx1 = 1;
    context_indx2 = 2;

    // Move the context forward if gaps are encountered to ensure that the
    // context pattern is matched to the sequence that the APOBEC enzyme
    // would have encountered.
    while ( the_seq[ window_start_i + context_indx1 ] == '-' ) {
      ++context_indx1;
      ++context_indx2;
      if (window_start_i + context_indx2 >= cons.size()){
        break;
      }
    }
    if (window_start_i + context_indx2 >= cons.size()){
      continue;
    }
    while ( the_seq[ window_start_i + context_indx2 ] == '-' ) {
      ++context_indx2;
      if (window_start_i + context_indx2 >= cons.size()){
        break;
      }
    }
    if (window_start_i + context_indx2 >= cons.size()){
      continue;
    }

    // Check for hypermutated spots
    if( ( cons[window_start_i + 0 ] == 'G' ) &&
          // Reference must mutate from G
        ( (the_seq[window_start_i + context_indx1 ] == 'A') ||
          (the_seq[window_start_i + context_indx1 ] == 'G') ) &&
          // Context position 1 must match R = [AG] in query
        ( ( the_seq[window_start_i + context_indx2 ] == 'A' ) ||
          ( the_seq[window_start_i + context_indx2 ] == 'G' ) ||
          ( the_seq[window_start_i + context_indx2 ] == 'T' ) ) ) {

      num_potential_mut++;

      if ( the_seq[ window_start_i + 0 ] == 'A'){
        hypermuted = true;
        num_mut++;
        all_mut_pos.push_back(window_start_i);
        if (fix){
          the_seq[ window_start_i ] = fix_with;
        }
      }

    }

//    if( ( cons[window_start_i + 0 ] == "G" ) && 
//          # Reference must mutate from G
//        ( the_seq[window_start_i + context_indx1 ] %in% c( "A", "G" ) ) && 
//          # Context position 1 must match R = [AG] in query
//        ( the_seq[window_start_i + context_indx2 ] %in% c( "A", "G", "T" ) ) ){ 
//          # Context position 2 must match D = [AGT] in query
//        num_potential_mut <- num_potential_mut + 1;
//        hyper_muted <- as.character( the_seq[ window_start_i + 0 ] ) == "A"
//        all_mut_pos <- rbind(all_mut_pos,
//          data.frame(pos = window_start_i,
//                     base_in_query = as.character( the_seq[ window_start_i + 0 ] ),
//                     full_seq = paste(as.character( the_seq[ window_start_i + c(0, context_indx1, context_indx2) ] ), 
//                                      sep = '', collapse = ''),
//                     type = 'mut',
//                     muted = hyper_muted,
//                     stringsAsFactors = F)
//          )
//        if( hyper_muted ) { 
//            # If G -> A mutation occurred
//            num_mut <- num_mut + 1;
//            if ( fix_with != FALSE ) {
//              the_seq[ window_start_i ] <- fix_with
//            }
//        }
//    }


  
  }

//  // hacky slashy for dev
//  all_mut_pos.push_back(1);
//  all_mut_pos.push_back(2);
//  all_mut_pos.push_back(3);
//  all_mut_pos.push_back(4);
//  all_mut_pos.push_back(5);

  Rcpp::List result;
  result = Rcpp::List::create(
    Rcpp::Named("num_mut") = num_mut,
    Rcpp::Named("num_potential_mut") = num_potential_mut,
    Rcpp::Named("num_control") = num_control,
    Rcpp::Named("num_potential_control") = num_potential_control,
    Rcpp::Named("p_value") = 0.1,
    Rcpp::Named("all_mut_pos") = all_mut_pos,
    Rcpp::Named("the_seq") = the_seq
  );

  return result;
}

