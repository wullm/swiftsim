#include "vector_tests.h"

int main(){

#ifdef HAVE_AVX512_F
printf("No mask tests for AVX512\n");
return 0;
#else

  mask_t m;
  vec_zero_mask(m);
  for(int i = 0; i < VEC_SIZE; i++){
    if(m.i[i] != 0){
        return 1;
    }
  }


    return 0;
#endif
}


