#include "vector_tests.h"

int main(){


#ifdef HAVE_AVX512_F
  printf("Not yet defined mask tests for AVX512\n");
  return 0;
#else
  mask_t mask;
  vec_init_mask_true(mask);

  for(int i = 0; i < VEC_SIZE; i++){
    if(mask.i[i] != 0xFFFFFFFF){
         return 1;
    }
  }


    return 0;
#endif
}


