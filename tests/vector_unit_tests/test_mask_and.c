#include "vector_tests.h"

int main(){

#ifdef HAVE_AVX512_F
    printf("Not yet implemented mask tests for AVX512\n");
    return 0;
#else

    vector result;
    mask_t mask, mask2;
    for(int i = 0; i < VEC_SIZE; i++){
        mask2.i[i] = 0xFFFFFFFF;
        if(i & 1){
            mask.i[i] = 0xFFFFFFFF;
        }else{
            mask.i[i] = 0;
        }
    }

    result.v = vec_mask_and(mask2, mask);

    for(int i = 0; i < VEC_SIZE; i++){
        if(i & 1 && result.i[i] != 0xFFFFFFFF){
            return 1;
        }else if( !(i&1) && result.i[i] != 0){
            return 1;
        }
    }


    return 0;
#endif
}


