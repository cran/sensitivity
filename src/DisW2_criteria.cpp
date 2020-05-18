#include <Rcpp.h>
using namespace Rcpp;
using namespace std;


void DisW2_perElement_AbsDiff(NumericVector Xi,NumericVector Xk,NumericVector LocalVecXi_k, int size, int i, int k){
      
      for (int elem=0;elem<size;elem++){
          LocalVecXi_k[elem] = 1.5 - abs(Xi[i*size+elem]-Xk[k*size+elem]) * (1.0 - abs(Xi[i*size+elem]-Xk[k*size+elem]));
      }
  }


// [[Rcpp::export]] 
double DisW2_Crossprod(NumericVector X, int d){
  
    int n = X.size();
    int l = n/d;
    double s = 0.0;
    
    NumericVector Xiptr = clone(X);
    NumericVector Xkptr = clone(X);
    NumericVector LocalVecXi_k(d);

    for (int i=0;i<l;i++){
//         cout << "i= :" << i << endl;

        Xkptr = clone(X);

        for (int k=i;k<l;k++){
//            cout << "k= :" << k << endl;
            if (i!=k){
                DisW2_perElement_AbsDiff(Xiptr,Xkptr,LocalVecXi_k, d, i, k);
                s += 2.0*std::accumulate(LocalVecXi_k.begin(), LocalVecXi_k.end(), 1.0, std::multiplies<double>());
            } else {
                DisW2_perElement_AbsDiff(Xiptr,Xkptr,LocalVecXi_k, d, i, k);
                s += std::accumulate(LocalVecXi_k.begin(), LocalVecXi_k.end(), 1.0, std::multiplies<double>());
            }
        }
    }
    return s;
}

