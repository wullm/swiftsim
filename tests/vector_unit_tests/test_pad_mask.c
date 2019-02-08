#include "vector_tests.h"

int main(){

#ifdef HAVE_AVX512_F
  printf("No mask tests for AVX512\n");
  return 0;
#else

mask_t m;
int pad = 2;
for(int i = 0; i < VEC_SIZE; i++){
  m.i[i] = 0xFFFFFFFF;
}
vec_pad_mask(m, pad);

for(int i = 0; i < VEC_SIZE; i++){
  if(i < VEC_SIZE-pad && m.i[i] != 0xFFFFFFFF){
     return 1;
  }else if( i >= VEC_SIZE - pad && m.i[i] != 0){
     return 1; 
  }
}

    return 0;
#endif
}


