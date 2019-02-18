#include "vector_tests.h"

int main(){

  mask_t result;
  vector v;
  vector zero;
  for(int i = 0; i < VEC_SIZE; i++){
    if(i & 1){
      v.f[i] = -1.0;
    }else{
      v.f[i] = 1.0;
    }
    zero.f[i] = -1.0;
  }

  result.v = vec_cmp_lte(v.v, zero.v);
  for(int i = 0; i < VEC_SIZE; i++){
    if( i & 1 && result.i[i] != -1 ){
        return 1;
    }else if( i & 1 == 0 && result.i[i] != 0){
        return 1;
    }
  } 

  return 0;
}
