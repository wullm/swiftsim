#include "vector_tests.h"

int main(){

    vector v;
    for(int i = 0; i < VEC_SIZE/2; i++){
        v.d[i] = 1.0/(double)(i+1);
    }
    vector result;
    result.vd = vec_dbl_rcp(v.vd);
    for(int i = 0; i < VEC_SIZE/2; i++){
        if((result.d[i] <= (double)(i)) || (result.d[i] >= (double)(i+2))){
            printf("%i %f\n", i, result.d[i]);
            return 1;
        }
    }



    return 0;
}


