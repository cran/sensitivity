#include <Rcpp.h>
using namespace Rcpp;
using namespace std;


double eucl_dist_op(double i, double j){
  
  return std::pow((i-j),2);
}


// [[Rcpp::export]] 
double maximin_cpp(NumericMatrix X){
  
  const int d = X.ncol();
  const int l = X.nrow();
  double s = 0.0;
  NumericVector Xiptr(d);
  NumericVector Xkptr(d);
  NumericVector LocalVecXi_k(d);
  NumericVector Stock(l);
  NumericVector LocalMin(l);

  for (int i=0;i<l;i++){

      Xiptr = X(i,_);
      Stock(i) = std::sqrt(d);

      for (int k=0;k<l;k++){
        if(i!=k){
        Xkptr = X(k,_);
        std::transform(Xiptr.begin(),Xiptr.end(), Xkptr.begin(), LocalVecXi_k.begin(),eucl_dist_op);
        Stock(k) = std::sqrt(std::accumulate(LocalVecXi_k.begin(), LocalVecXi_k.end(), 0.0));
        }
      }
  LocalMin(i)=*std::min_element(Stock.begin(),Stock.end());
  }
  s = *std::min_element(LocalMin.begin(),LocalMin.end());
  return s;
}
