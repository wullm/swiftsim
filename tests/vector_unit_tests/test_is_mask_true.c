#include "vector_tests.h"

int main(){

#ifdef HAVE_AVX512_F
    printf("No mask tests for AVX512 yet.\n");
    return 0;
#else

    vector v1,v2; 
    for(int i = 0; i < VEC_SIZE; i++){
      v1.f[i] = (float)i;
      if(i & 1){
          v2.f[i] = (float)i;
      }else{
          v2.f[i] = -1.0;
      }
    }
    int result = vec_is_mask_true(vec_cmp_lte(v1.v, v2.v));
    for(int i = 0; i < VEC_SIZE; i++){
      if((i & 1) && !(result & (1 << i)) ){
          return 1;
      }else if( !(i&1) && result & (1<<i)){
          return 1;
      }
    }


    return 0;
#endif
}


