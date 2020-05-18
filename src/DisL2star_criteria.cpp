#include <Rcpp.h>
using namespace Rcpp;
using namespace std;


double DisL2star_super_Operation(double i){
        return 1.0 - i;
    }
    
double DisL2star_super_Operation_Magic(double i){
        return 1.0 - i*i;
    }

void DisL2star_perElement_Max(NumericVector Xi,NumericVector Xk,NumericVector LocalVecXi_k, int size, int i, int k){
      
      for (int elem=0;elem<size;elem++){
          LocalVecXi_k[elem] = 1.0 - std::max(Xi[i*size+elem],Xk[k*size+elem]);
      }
  }


// [[Rcpp::export]] 
double DisL2star_Crossprod(NumericVector X, int d){
  
    int n = X.size();
    int l = n/d;
    double s = 0.0;
    double t1;
    double t2;
    double t;
    
    NumericVector Xiptr = clone(X);
    NumericVector Xkptr = clone(X);
    NumericVector LocalVecXi_k(d);
    NumericVector v1(d);
    NumericVector v2(d);

    for (int i=0;i<l;i++){
//         cout << "i= :" << i << endl;

        Xkptr = clone(X);

        for (int k=i;k<l;k++){
//            cout << "k= :" << k << endl;
            if (i!=k){
                DisL2star_perElement_Max(Xiptr,Xkptr,LocalVecXi_k, d, i, k);
                t = 2.0*std::accumulate(LocalVecXi_k.begin(), LocalVecXi_k.end(), 1.0, std::multiplies<double>())/(l*l);
            } else {
                std::transform(Xiptr.begin()+i*d,Xiptr.begin()+(i+1)*d, v1.begin(), DisL2star_super_Operation);
                t1 = std::accumulate(v1.begin(), v1.end(), 1.0, std::multiplies<double>());
                std::transform(Xiptr.begin()+i*d,Xiptr.begin()+(i+1)*d, v2.begin(), DisL2star_super_Operation_Magic);
                t2 = std::accumulate(v2.begin(), v2.end(), 1.0, std::multiplies<double>());
                t = t1/(l*l) - std::pow(2,1 - d)/l * t2;
            }
            
            s += t;

        }
    }
    return s;
}
