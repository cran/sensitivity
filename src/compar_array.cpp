#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
bool Compar_array(IntegerMatrix X, IntegerMatrix Y) {
  
    int cX = X.ncol();
    int cY = Y.ncol();
      
    for(int i=0;i<cX;i++){
      
      int* bX = X(_,i).begin();
      int* eX = X(_,i).end();
        
      for(int j=0;j<cY;j++){

        if(std::equal(bX,eX,Y(_,j).begin())){
          
          return true;
        }
      }     
    }
    
    return false;
}
