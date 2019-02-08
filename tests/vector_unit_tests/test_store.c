#include "vector_tests.h"

int main(){

vector v1;
for(int i = 0; i < VEC_SIZE; i++){
  v1.i[i] = i;
}
vector z1;
vec_store(v1.v, &(z1.f));

for(int i = 0; i < VEC_SIZE; i++){
    if(z1.i[i] != i){
        return 1;
    }
}


    return 0;
}


