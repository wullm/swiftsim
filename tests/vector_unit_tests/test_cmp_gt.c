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
    zero.f[i] = 0.0;
  }
  vec_create_mask(result, vec_cmp_gt(v.v, zero.v));

  for(int i = 0; i < VEC_SIZE; i++){
    if(i & 1 && vec_is_mask_bit_true(result, i)){
      return 1;
    }else if(!(i & 1) && !vec_is_mask_bit_true(result, i)){
      return 1;
    }
  }

    return 0;
}
