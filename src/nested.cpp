// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
using namespace std;


int ceil_of_div(double i, double j){
  
  return std::ceil(i/j);
}


// [[Rcpp::export]] 
IntegerVector nested_permu_cplus(IntegerVector layers){
  
  const int u = layers.size();
  const int n = std::accumulate(layers.begin(), layers.end(), 1, std::multiplies<int>());  
  int p = 0;
  IntegerVector m(u+1, 0);
  IntegerVector t(u, 0);
  IntegerVector res(n, 0);
  
  std::partial_sum(layers.begin(), layers.end(), m.begin()+1, std::multiplies<int>());
  std::transform(m.begin()+1, m.end(), t.begin(),std::bind1st(std::divides<int>(), n));
  
  for (int i=0;i<u;i++){
    
    std::vector<int> diff;
    IntegerVector C(m(i), 0);
    IntegerVector perm(m(i+1)-m(i), 0);
    IntegerVector v = Rcpp::seq(1, m(i+1));
    std::transform(res.begin(), res.begin()+m(i), C.begin(), std::bind(ceil_of_div,std::placeholders::_1,t(i)));
    std::sort(C.begin(),C.end());
    std::set_difference(v.begin(), v.end(), C.begin(), C.end(), std::inserter(diff, diff.begin()));
    perm = RcppArmadillo::sample(diff, m(i+1)-m(i), false, NumericVector::create());
    
    if (i==u){
      std::copy(perm.begin(), perm.end(), res.begin()+m(i+1));
      return(res);
    } 
    
    for (int j=m(i);j<m(i+1);j++){

      p = perm(j-m(i));
      IntegerVector w = Rcpp::seq((p-1)*t(i)+1, p*t(i));
      res(j) = Rcpp::as<int>(wrap(RcppArmadillo::sample(w, 1, false, NumericVector::create())));
    }    
  } 
  return(res);
}


//[[Rcpp::export]] 
IntegerMatrix nested_lhs_cplus(int d, IntegerVector layers){
  
  const int n = std::accumulate(layers.begin(), layers.end(), 1, std::multiplies<int>());
  IntegerMatrix out(n,d);

  for (int i=0;i<d;i++){
    out(_,i) = nested_permu_cplus(layers);
  }
  return(out);
}
  
  
  
  
