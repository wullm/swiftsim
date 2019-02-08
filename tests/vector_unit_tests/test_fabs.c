#include "vector_tests.h"

int main(){

  vector v, result;
  for(int i = 0; i < VEC_SIZE; i++){
    v.f[i] = -1.0f * (float)i;
  }
  result.v = vec_fabs(v.v);
  for(int i = 0; i < VEC_SIZE; i++){
    if(result.f[i] != (float)i){
        return 1;
    }
  }


    return 0;
}


