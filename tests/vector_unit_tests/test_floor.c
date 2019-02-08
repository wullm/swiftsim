#include "vector_tests.h"

int main(){

    vector v, result;
    for(int i = 0; i < VEC_SIZE; i++){
        if(i & 1){
          v.f[i] = ((float)i) + 0.27f;
        }else{
          v.f[i] = ((float)i) + 0.72f;
        }
    }

    result.v = vec_floor(v.v);
    for(int i = 0; i < VEC_SIZE; i++){
        if(result.v[i] != (float)i){
            return 1;
        }
    }



    return 0;
}


