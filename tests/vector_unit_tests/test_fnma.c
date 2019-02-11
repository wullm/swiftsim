#include "vector_tests.h"

int main(){


    vector v1, v2,v3, result;
    for(int i = 0; i < VEC_SIZE; i++){
      v1.f[i] = (float)i;
      v2.f[i] = (float)i;
      v3.f[i] = (float)i;
    }
    result.v = vec_fnma(v1.v, v2.v, v3.v);

    for(int i = 0; i < VEC_SIZE; i++){
      float temp = (float)i;
      if(result.f[i] != (temp*temp)-temp){
          return 1;
      }
    }

    return 0;
}


