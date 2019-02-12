#include "vector_tests.h"

int main(){

  vector v1, v2;
  mask_t mask;

  for(int i = 0; i < VEC_SIZE; i++){
    v1.i[i] = 0xFFFF0000;
    v2.i[i] = 0x0000FFFF; 
#ifdef HAVE_AVX512_F
  } 
  mask = 0xAAAA;
#else
   if(i & 1){
      mask.i[i] = 0xFFFFFFFF;
    }else{
      mask.i[i] = 0;
    }
  }
#endif

  vector result;
  result.v = vec_blend(mask,v1.v, v2.v);
  
  for(int i = 0; i < VEC_SIZE; i++){
    if( i&1 && result.i[i] != 0x0000FFFF){
      return 1;
    }else if ( !(i&1) && result.i[i] != 0xFFFF0000){
      return 1;
    }
  }

  return 0;
}
