#include "vector_tests.h"

int main(){

    vector v1, v2, result;
    for(int i = 0; i < VEC_SIZE; i++){
        v1.f[i] = (float)i;
        if(i & 1){
            v2.f[i] = 1.0*(float)i;
        }else{
            v2.f[i] = -1.0*(float)i;
        }
    }

    result.v = vec_and(v1.v, v2.v);

    for(int i = 0; i < VEC_SIZE; i++){
        if(result.f[i] != (float)i){
	    printf("%f %f\n", result.f[i], (float)i);
            return 1;
        }
    }


    return 0;
}
