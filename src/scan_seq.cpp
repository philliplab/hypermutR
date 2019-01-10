
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

Rcpp::List rcpp_scan_seq(CharacterVector cons, CharacterVector the_seq, 
  CharacterVector the_pattern, LogicalVector fix_with){

  int num_potential_mut = 0;
  int num_mut = 0;
  int num_potential_control = 0;
  int num_control = 0;
  std::vector< int > all_mut_pos;

  // hacky slashy for dev
  all_mut_pos.push_back(1);
  all_mut_pos.push_back(2);
  all_mut_pos.push_back(3);
  all_mut_pos.push_back(4);
  all_mut_pos.push_back(5);

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

