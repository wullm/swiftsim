#include "vector_tests.h"

int main(){

  vector v;
  for(int i = 0; i < VEC_SIZE; i++){
    v.f[i] = (float)i;
  }

  vector result;
  result.m = vec_ftoi(v.v);
  for(int i = 0; i < VEC_SIZE; i++){
    if(result.i[i] != i){
      return 1;
    }
  }

    return 0;
}


