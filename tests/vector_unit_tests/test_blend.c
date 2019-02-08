#include "vector_tests.h"

int main(){

#ifdef HAVE_AVX512_F
printf("No mask tests for AVX512\n");
return 0;
#else

    vector v1, v2;
    mask_t m;

    for(int i = 0; i < VEC_SIZE; i++){
       v1.i[i] = 0xFFFF0000;
       v2.i[i] = 0x0000FFFF;
       if(i & 1){
       m.i[i] = 0xFFFFFFFF;
       }else{
        m.i[i] = 0;
       }
    }
    vector result;
    result.v = vec_blend(m,v1.v, v2.v);
    for(int i = 0; i < VEC_SIZE; i++){
      if( i&1 && result.i[i] != 0x0000FFFF){
        return 1;
      }else if ( !(i&1) && result.i[i] != 0xFFFF0000){
        return 1;
      }
    }




    return 0;
#endif
}


