#include <Rcpp.h>
#include <string>
#include <iostream>
using namespace Rcpp; // To not write Rcpp:: everytime

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
NumericVector cpp_get_indices(NumericMatrix& data, const IntegerMatrix& RP, const IntegerMatrix& I, IntegerVector& bootsample, int d) {
  
  // instantiation
  const int N = RP.nrow();
  const int cRP = RP.ncol();
  const double cst_sqrt = sqrt(2.0);
  int i_1, i_2;
  double val1, val2, valsum, valprod, terma, termb, termc, termean, termean2, iquot, iiquot;

  // output
  NumericVector out(cRP, 0.0);

  // estimate indices
  for(int k = 0; k < cRP; k++){

    terma = 0.0;
    termb = 0.0;
    termc = 0.0;
    termean = 0.0;
    
    for (int i = 0; i < N; i++) {
      
      // reorder given RP
      i_1 = bootsample[i] - 1;
      i_2 = RP(i_1,k) - 1;
      val1 = data(i_1,0);
      val2 = data(i_2,1);
      valsum = val1 + val2;
      valprod = val1 * val2;
      iquot = 1.0/(double)(i+1);
      iiquot = i * iquot;
      
      // recursive formula of Janon estimator (see doc)
      
      //update cov and var terms
      terma *= iiquot;
      termb *= iiquot;
      termc *= iiquot;
      termean2 = std::pow(termean, 2.0);
      terma += iiquot * termean2;
      termb += iiquot * termean2;
      
      //update mean
      termean *= iiquot;
      termean += valsum/(double)(2*(i+1));
      termean2 = std::pow(termean, 2.0);
      
      //update cov and var terms
      terma += iquot * valprod - termean2;
      termb += iquot * std::pow(valsum/cst_sqrt, 2.0) - termean2;
      termc += iquot * valprod;
    }
    
    // quotient evaluation
    out[k] = terma/(termb-termc);
    
    // deal with second-order indices
    if (k >= d){
      out[k] += - out[ I(k-d,0) - 1] - out[ I(k-d,1) - 1 ];
    }
  }

  return out;
}


// [[Rcpp::export]]
NumericVector cpp_get_total_indices(NumericMatrix& data, IntegerVector& bootsample) {
  
  // instantiation
  const int N = data.nrow();
  const int d = data.ncol()-1;
  const double cst_sqrt = sqrt(2.0);
  int i_1;
  double val1, val2, valsum, valprod, terma, termb, termc, termean, termean2, iquot, iiquot;
  
  // output
  NumericVector out(d, 0.0);
  
  // estimate indices
  for(int k = 0; k < d; k++){
    
    terma = 0.0;
    termb = 0.0;
    termc = 0.0;
    termean = 0.0;
    
    for (int i = 0; i < N; i++) {
      
      // reorder given RP
      i_1 = bootsample[i] - 1;
      val1 = data(i_1,0);
      val2 = data(i_1,k+1);
      valsum = val1 + val2;
      valprod = val1 * val2;
      iquot = 1.0/(double)(i+1);
      iiquot = i * iquot;
      
      // recursive formula of Janon estimator (see doc)
      
      //update cov and var terms
      terma *= iiquot;
      termb *= iiquot;
      termc *= iiquot;
      termean2 = std::pow(termean, 2.0);
      termb += iiquot * termean2;
      
      //update mean
      termean *= iiquot;
      termean += valsum/(double)(2*(i+1));
      termean2 = std::pow(termean, 2.0);
      
      //update cov and var terms
      terma += 1/2.0 * iquot * std::pow(val1-val2, 2.0);
      termb += iquot * std::pow(valsum/cst_sqrt, 2.0) - termean2;
      termc += iquot * valprod;
    }
    
    // quotient evaluation
    out[k] = terma/(termb-termc);
  }
  
  return out;
}
