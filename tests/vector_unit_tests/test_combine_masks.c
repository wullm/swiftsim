#include "vector_tests.h"

int main(){

#ifdef HAVE_AVX512_F
  printf("No mask tests for AVX512\n");
  return 0;
#else

mask_t m1, m2;
for(int i = 0; i < VEC_SIZE; i++){
   m1.i[i] = 0xFFFFFFFF;
  if(i == 0){
    m2.i[i] = 0xFFFFFFFF;
  }else{
    m2.i[i] = 0;
  }
}

vec_combine_masks(m1, m2);

for(int i = 0; i < VEC_SIZE; i++){

    if(i == 0 && m1.i[i] != 0xFFFFFFFF){

        return 1;
    }else if( i > 0 && m1.i[i] != 0){
        return 1;
    }

}

    return 0;
#endif
}


