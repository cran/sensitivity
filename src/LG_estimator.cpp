#include <iostream>
#include <stdio.h>
#include <math.h>
#include <R.h>

using namespace std;


extern "C"
{
  void LG_estimator(const double* Y, const double* mean, const int*d, const int* N, const int* ind, const int* ind2, double* a, double *b, double* c, double* S)
  { 
    int i_1;
    int i_2;
    //cout<<"Sorting the array"<<endl;
    for(int k = 0; k < *d; k++){
      //cout<<"k="<< k <<endl;
      for (int i = 0; i < *N; i++) {
        //cout<<"i="<< i <<endl;
        i_1 = i+(ind[k]-1)*(*N);
        i_2 = i+(ind2[k]-1)*(*N);
        a[k] += (Y[i_1]-mean[k])*(Y[i_2]-mean[k]);
        b[k] += pow((Y[i_1]+Y[i_2])/sqrt(2.0)-(sqrt(2.0)+1)*mean[k],2.0);
        c[k] += Y[i_1]*Y[i_2];
      }
        S[k] = a[k]/(b[k]-c[k]);
    }
  }
}

// extern "C"
// {
//   void LG_estimator(double* Y, double* mean, int*d, int* N, int* ind, int* ind2, double* a, double* S)
//   {
//     //cout<<"Sorting the array"<<endl;
//     for(int k = 0; k < *d; k++){
//       //cout<<"k="<< k <<endl;
//       for (int i = 0; i < *N; i++) {
//         //cout<<"i="<< i <<endl;
//         a[k] += pow((Y[i+(ind[k]-1)*(*N)]+Y[i+(ind2[k]-1)*(*N)])/sqrt(2.0)-(sqrt(2.0)+1)*mean[k],2.0);
//       }
//       S[k] = a[k];
//     }
//   }
// }
