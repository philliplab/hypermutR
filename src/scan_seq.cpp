
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

Rcpp::List rcpp_scan_seq_int(CharacterVector r_cons, CharacterVector r_the_seq, 
  CharacterVector r_the_pattern, CharacterVector r_fix_with){

  int num_potential_mut = 0;
  int num_mut = 0;
  int num_potential_control = 0;
  int num_control = 0;
  
  std::vector< int > amp_pos;
  std::vector< char > amp_base_in_query;
  std::vector< std::string > amp_full_seq;
  std::vector< std::string > amp_type;
  std::vector< std::string > amp_muted;

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
      amp_pos.push_back(window_start_i);
      amp_base_in_query.push_back(the_seq[ window_start_i ]);
      amp_full_seq.push_back( the_seq.substr( window_start_i, 3) );
      amp_type.push_back( "mut" );

      if ( the_seq[ window_start_i + 0 ] == 'A'){
        hypermuted = true;
        num_mut++;
        amp_muted.push_back( "yes" );
        if (fix){
          the_seq[ window_start_i ] = fix_with;
        }
      } else {
        amp_muted.push_back( "no" );
      }
    }

    // Check for control spots
    if (( cons[ window_start_i + 0 ] == 'G' ) && 
      // Reference must mutate from G
      ((( the_seq[ window_start_i + context_indx1 ] == 'C' ) || 
      ( the_seq[ window_start_i + context_indx1 ] == 'T' )) && 
      // Option 1 Context position 1 must match Y = [CT] in query
      (( the_seq[ window_start_i + context_indx2 ] == 'A' ) || 
      ( the_seq[ window_start_i + context_indx2 ] == 'C' ) || 
      ( the_seq[ window_start_i + context_indx2 ] == 'G' ) || 
      ( the_seq[ window_start_i + context_indx2 ] == 'T' ))) || 
      // Option 1 Context position 2 must match N = [ACGT] in query
      ((( the_seq[ window_start_i + context_indx1 ] == 'A' ) || 
      ( the_seq[ window_start_i + context_indx1 ] == 'G' )) && 
      // Option 2 Context position 1 must match R = [AG] in query
      ( the_seq[ window_start_i + context_indx2 ] == 'C' ))){ 
      // Option 2 Context position 2 must match C

      num_potential_control++;
      amp_pos.push_back(window_start_i);
      amp_base_in_query.push_back(the_seq[ window_start_i ]);
      amp_full_seq.push_back( the_seq.substr( window_start_i, 3) );
      amp_type.push_back( "mut" );
  
      if ( the_seq[ window_start_i + 0 ] == 'A' ){
        num_control++;
        amp_muted.push_back( "yes" );
      } else {
        amp_muted.push_back( "no" );
      }
    }
  }

  Rcpp::List all_mut_pos;
  all_mut_pos = Rcpp::List::create(
    Rcpp::Named("pos") = amp_pos,
    Rcpp::Named("base_in_query") = amp_base_in_query,
    Rcpp::Named("full_seq") = amp_full_seq,
    Rcpp::Named("type") = amp_type,
    Rcpp::Named("muted") = amp_muted
  );

  Rcpp::List result;
  result = Rcpp::List::create(
    Rcpp::Named("num_mut") = num_mut,
    Rcpp::Named("num_potential_mut") = num_potential_mut,
    Rcpp::Named("num_control") = num_control,
    Rcpp::Named("num_potential_control") = num_potential_control,
    Rcpp::Named("all_mut_pos") = all_mut_pos,
    Rcpp::Named("the_seq") = the_seq
  );

  return result;
}

