#include "vector_tests.h"

int main(){

  mask_t mask;
  int pad = 2;
  vec_init_mask_true(mask);

  vec_pad_mask(mask, pad);

  for(int i = 0; i < VEC_SIZE; i++){
    if(i < VEC_SIZE-pad && !vec_is_mask_bit_true(mask, i)){
      return 1;
    }else if( i >= VEC_SIZE - pad && vec_is_mask_bit_true(mask, i)){
      return 1; 
    }
  }

  return 0;
}
