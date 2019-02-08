#include "vector_tests.h"

int main(){

    vector v;
    for(int i = 0; i < VEC_SIZE; i++){
        v.f[i] = 1.0/(float)(i+1);
    }
    vector result;
    result.v = vec_rcp(v.v);
    for(int i = 0; i < VEC_SIZE; i++){
        if((result.f[i] <= (float)(i)) || (result.f[i] >= (float)(i+2))){
            printf("%i %f\n", i, result.f[i]);
            return 1;
        }
    }



    return 0;
}


