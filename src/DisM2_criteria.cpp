#include <Rcpp.h>
using namespace Rcpp;
using namespace std;


void DisM2_perElement_Max(NumericVector Xi,NumericVector Xk,NumericVector LocalVecXi_k, int size, int i, int k){
      
      for (int elem=0;elem<size;elem++){
          LocalVecXi_k[elem] = 2.0 - std::max(Xi[i*size+elem],Xk[k*size+elem]);
      }
  }


// [[Rcpp::export]] 
double DisM2_Crossprod(NumericVector X, int d){
  
    int n = X.size();
    int l = n/d;
    double s = 0.0;
    
    NumericVector Xiptr = clone(X);
    NumericVector Xkptr = clone(X);
    NumericVector LocalVecFinal(d);

    for (int i=0;i<l;i++){
//         cout << "i= :" << i << endl;

        Xkptr = clone(X);

        for (int k=i;k<l;k++){
//            cout << "k= :" << k << endl;
            if(k!=i){
                DisM2_perElement_Max(Xiptr,Xkptr,LocalVecFinal, d, i, k);            
                s += 2.0*std::accumulate(LocalVecFinal.begin(), LocalVecFinal.end(), 1.0, std::multiplies<double>());
            } else {
                DisM2_perElement_Max(Xiptr,Xkptr,LocalVecFinal, d, i, k);            
                s += std::accumulate(LocalVecFinal.begin(), LocalVecFinal.end(), 1.0, std::multiplies<double>());
            }
        }
    }
    return s;
}

// [[Rcpp::export]]
double DisM2_Rowprod(NumericVector x, int d){
  int n = x.size();
  int l = n/d;
  double out = 0.0;
  
//  std::cout << n << ' ' << l << std::endl;
//  NumericVector out(l);

  for (int i = 0; i < l; i++) {
    out += std::accumulate(x.begin()+i*d, x.begin()+(i+1)*d, 1.0, std::multiplies<double>());
  }
  return(out);
}
