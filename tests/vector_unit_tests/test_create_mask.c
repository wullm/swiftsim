#include "vector_tests.h"

int main(){

#ifdef HAVE_AVX512_F
  printf("test_create_mask not yet implemented for AVX512_F\n");
  return 1;
#else
  mask_t masky;
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
  vec_create_mask(masky, vec_cmp_gt(v.v, zero.v));
 for(int i = 0; i < VEC_SIZE; i++){
    if( i & 1 && masky.i[i] != 0 ){
        return 1;
    }else if( !(i & 1) && masky.i[i] != -1){
        return 1;
    }
 } 

    return 0;
#endif
}


