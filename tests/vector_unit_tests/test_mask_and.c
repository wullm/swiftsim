#include "vector_tests.h"

int main(){

  mask_t result;
  mask_t mask, mask2;

  vec_init_mask_true(mask2);

#ifdef HAVE_AVX512_F
  mask = 0xAAAA;
#else
  for(int i = 0; i < VEC_SIZE; i++){
    if(i & 1){
      mask.i[i] = 0xFFFFFFFF;
    }else{
      mask.i[i] = 0;
    }
  }
#endif

  vec_create_mask(result, vec_mask_and(mask2, mask));

  for(int i = 0; i < VEC_SIZE; i++){
    if(i & 1 && !vec_is_mask_bit_true(result, i)){
      return 1;
    }else if( !(i&1) && vec_is_mask_bit_true(result, i)){
      return 1;
    }
  }

  return 0;
}
