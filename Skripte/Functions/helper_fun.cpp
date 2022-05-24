#include <Rcpp.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
using namespace Rcpp;



// [[Rcpp::export]]
std::vector<int> combinations(int x) {
  std::vector<int> out;
  for(int a = 1; a <=x; a = a + 1 ) {
    for(int b = 1; b <= x; b = b + 1 ) {
      if( b!= a) {
        out.push_back(a);
        out.push_back(b);
      }
      
    }
  }
  return out;
}

// [[Rcpp::export]]
int findrow_old(NumericMatrix a, NumericVector b){
  int length_n;
  bool c = true;
  length_n = a.nrow();
  for (int row=0; row<length_n; row++)
  {
    c = (a(row,0) == b(0)) & (a(row,1) == b(1));
    if(c){
      return row + 1;
      break;
    }
  }
  return 0;
}

// [[Rcpp::export]]
int findrow2(NumericVector b,NumericMatrix a){
  int length_n;
  bool c = true;
  length_n = a.nrow();
  for (int row=0; row<length_n; row++)
  {
    c = (a(row,0) == b(0)) & (a(row,1) == b(1));
    if(c){
      return row + 1;
      break;
    }
  }
  return 0;
}


/*** R



*/