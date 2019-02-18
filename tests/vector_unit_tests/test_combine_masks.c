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

    if(i == 0 && m1.i[i] != 0xFFFFFFFF){
        return 1;
    }else if( i > 0 && m1.i[i] != 0){
        return 1;
    }

  }

  return 0;
}
