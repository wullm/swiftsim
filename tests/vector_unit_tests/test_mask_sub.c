#include "vector_tests.h"

int main(){

  vector v1, v2;
  mask_t mask;

  for(int i = 0; i < VEC_SIZE; i++){
    v1.f[i] = (float)i;
    v2.f[i] = (float)i;
#ifdef HAVE_AVX512_F
  }
  mask = 0x7FFF;
#else
    if(i < VEC_SIZE-1){
      mask.i[i] = 0xFFFFFFFF;
    }else{
      mask.i[i] = 0;
    }
  }
#endif

  vector z;
  z.v = vec_mask_sub(v1.v,v2.v, mask);

  for(int i = 0; i < VEC_SIZE; i++){
    if(i < VEC_SIZE-1 && z.f[i] != 0.0f){
      return 1;
    }
    if(i == VEC_SIZE-1 && z.f[i] != (float)i){
      return 2;
    }
  }


  return 0;
}
