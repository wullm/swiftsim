#include "vector_tests.h"

int main(){

  mask_t m1, m2;
  vec_init_mask_true(m1);

#ifdef HAVE_AVX512_F
  m2 = 1;
#else
  for(int i = 0; i < VEC_SIZE; i++){
    if(i == 0){
      m2.i[i] = 0xFFFFFFFF;
    }else{
      m2.i[i] = 0;
    }
  }
#endif

  vec_combine_masks(m1, m2);

  for(int i = 0; i < VEC_SIZE; i++){
    
    if(i == 0 && !vec_is_mask_bit_true(m1, i)){
      return 1;
    }else if( i > 0 && vec_is_mask_bit_true(m1, i)){
      return 1;
    }

  }

  return 0;
}
