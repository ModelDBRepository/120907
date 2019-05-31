/* This file is for my own math-helper functions which seems to be missing/
   or not functioning as expected. */
#include "my_math.h"

double min(double a, double b){
  if(a < b){
    return a;
  }else{
    return b;
  }
}

int my_round(double x){
  if(ceil(x)-x > x-floor(x)){
    return floor(x);
  }else{
    return ceil(x);
  }
}

unsigned long factorial(int N){
  int  n;
  unsigned long N_fact;
  n = N;
  N_fact = 1;
  while(n>1){
    if(N_fact > ULONG_MAX/n){
      return ULONG_MAX;
    }
    N_fact = N_fact*n;
    n = n-1;
  }

  return N_fact;
}
