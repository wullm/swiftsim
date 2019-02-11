#include "vector_tests.h"

int main(){

  vector v1, v2;
  for(int i= 0; i < VEC_SIZE/2; i++){
    v1.d[i] = 0.0f;
    v2.d[i] = (double)i;
  }

  vector result;
  result.vd = vec_dbl_fmin(v1.vd, v2.vd);
  for(int i = 0; i < VEC_SIZE/2; i++){
    if(result.d[i] != 0.0f){
        return 1;
    }
  }

    return 0;
}


