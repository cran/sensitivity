#include <Rcpp.h>
using namespace Rcpp;
using namespace std;


double  DisC2_super_Operation(double i){
        return 1.0 + abs(i - 0.5)*0.5;
    }
    
double  DisC2_super_Operation_Magic(double i){
      return 0.5 * abs(i - 0.5);
    }

void  DisC2_perElement_AbsDiff(NumericVector Xi,NumericVector Xk,NumericVector LocalVecXi_k, int size, int i, int k){

      for (int elem=0;elem<size;elem++){
          LocalVecXi_k[elem] = abs(Xi[i*size+elem] -Xk[k*size+elem]);
      }
  }

void  DisC2_perElement_Final(NumericVector Xi,NumericVector Xk,NumericVector LocalVecXi_k, NumericVector LocalVecFinal, int size){

      for (int elem=0;elem<size;elem++){
          LocalVecFinal[elem] = Xi[elem] + Xk[elem] - 0.5*LocalVecXi_k[elem];
      }
  }

// [[Rcpp::export]] 
double DisC2_Crossprod(NumericVector X, int d){
  
    int n = X.size();
    int l = n/d;
    double s = 0.0;
    
    NumericVector Xiptr = clone(X);
    NumericVector Xkptr = clone(X);
    NumericVector LocalVecXi(d);
    NumericVector LocalVecXk(d);
    NumericVector LocalVecXi_k(d);
    NumericVector LocalVecFinal(d);

    for (int i=0;i<l;i++){
//         cout << "i= :" << i << endl;

        std::transform(Xiptr.begin()+i*d,Xiptr.begin()+(i+1)*d, LocalVecXi.begin(), DisC2_super_Operation);        
        Xkptr = clone(X);

        for (int k=i;k<l;k++){
//            cout << "k= :" << k << endl;
            if(k!=i){
                std::transform(Xkptr.begin()+k*d,Xkptr.begin()+(k+1)*d, LocalVecXk.begin(), DisC2_super_Operation_Magic);
                DisC2_perElement_AbsDiff(Xiptr,Xkptr,LocalVecXi_k, d, i, k);            
                DisC2_perElement_Final(LocalVecXi , LocalVecXk , LocalVecXi_k, LocalVecFinal, d);
                s += 2.0*std::accumulate(LocalVecFinal.begin(), LocalVecFinal.end(), 1.0, std::multiplies<double>());
            } else {
                std::transform(Xkptr.begin()+k*d,Xkptr.begin()+(k+1)*d, LocalVecXk.begin(), DisC2_super_Operation_Magic);
                DisC2_perElement_AbsDiff(Xiptr,Xkptr,LocalVecXi_k, d, i, k);            
                DisC2_perElement_Final(LocalVecXi , LocalVecXk , LocalVecXi_k, LocalVecFinal, d);
                s += std::accumulate(LocalVecFinal.begin(), LocalVecFinal.end(), 1.0, std::multiplies<double>());
            }
        }
    }
    return s;
}

// [[Rcpp::export]]
double DisC2_Rowprod(NumericVector x, int d){
  int n = x.size();
  int l = n/d;
  double out = 0.0;

  for (int i = 0; i < l; i++) {
    out += std::accumulate(x.begin()+i*d, x.begin()+(i+1)*d, 1.0, std::multiplies<double>());
  }
  return(out);
}
