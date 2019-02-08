#include "vector_tests.h"

int main(){

  vector v1, v2;
  for(int i= 0; i < VEC_SIZE; i++){
    v1.f[i] = 0.0f;
    v2.f[i] = (float)i;
  }

  vector result;
  result.v = vec_fmax(v1.v, v2.v);
  for(int i = 0; i < VEC_SIZE; i++){
    if(result.f[i] != (float)i){
        return 1;
    }
  }

    return 0;
}


