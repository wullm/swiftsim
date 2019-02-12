#include "vector_tests.h"

int main(){

#ifdef HAVE_AVX512_F
  /* VEC_FORM_PACKED_MASK is not needed in AVX512 as an explicit instruction exists for VEC_LEFT_PACK so a packed mask is not needed.*/
  return 0;
#else

  int mask = 0xAAAA;
  mask_t packed_mask;

  VEC_FORM_PACKED_MASK(mask, packed_mask);
  
  for(int i = 0; i < VEC_SIZE; i++){
    if(i < VEC_SIZE / 2 && packed_mask.i[i] != 2*i + 1){
      printf("Lower mask fail\n");
      return 1;
    }else if(i >= VEC_SIZE / 2 && packed_mask.i[i] != 0){
      printf("Upper mask fail\n");
      return 2;
    }
  }

  return 0;
#endif
} 
