#include <iostream>
#include <cstdlib>
#include <climits>
#include <R.h>

using namespace std;

extern "C"
{
  int LG_compare(const void * a, const void * b)
  {
    if (*(double*)a > *(double*)b) return 1;
    else if (*(double*)a < *(double*)b) return -1;
    else return 0;  
  }
}

extern "C"
{
  void LG_rowsort(int* N, double* x, int*d, double*inter, double* out)
  {
    //cout<<"Sorting the array"<<endl;
    for (int i = 0; i < *N; i++) {
      //cout<<"i="<< i <<endl;
      for(int k = 0; k < *d; k++){
        //cout<<"pos="<< i*(*d)+k <<endl;
        inter[k] = x[i*(*d)+k];
      }
      //cout << "here" << endl;
      std::qsort(inter,*d,sizeof(double),LG_compare);
      for(int k = 0; k < *d; k++){
        out[i*(*d)+k] = inter[k];
      }
    }
  }
}

// extern "C"
// {
//   void rowsort(int* length_x, double* x, double* out)
//   {
//     std::qsort(x,*length_x,sizeof(double),compare);
//     for(int k = 0; k < *length_x; k++){
//         out[k] = x[k];
//     }
//   }
// }


