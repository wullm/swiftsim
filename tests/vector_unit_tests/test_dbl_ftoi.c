#include "vector_tests.h"

int main(){

  vector v;
  for(int i = 0; i < VEC_SIZE; i++){
    v.d[i] = (double)i;
  }

  vector result;
  result.m = vec_dbl_ftoi(v.vd);
  for(int i = 0; i < VEC_SIZE; i++){
    if(result.i[i] != i){
      return 1;
    }
  }

    return 0;
}


