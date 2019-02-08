#include "vector_tests.h"

int main(){

#ifdef HAVE_AVX512_F
    printf("Not yet implemented mask tests for AVX512\n");
    return 0;
#else

    vector v1, result;
    mask_t mask;
    for(int i = 0; i < VEC_SIZE; i++){
        v1.f[i] = (float)i;
        if(i & 1){
            mask.i[i] = 0xFFFFFFFF;
        }else{
            mask.i[i] = 0;
        }
    }

    result.v = vec_and_mask(v1.v, mask);

    for(int i = 0; i < VEC_SIZE; i++){
        if(i & 1 && result.f[i] != (float)i){
            return 1;
        }else if( !(i&1) && result.i[i] != 0){
            return 1;
        }
    }


    return 0;
#endif
}


