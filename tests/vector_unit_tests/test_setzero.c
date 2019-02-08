#include "vector_tests.h"

int main(){

vector v;
v.v = vec_setzero();
for(int i = 0; i < VEC_SIZE; i++){
  if(v.f[0] != 0.0f){
    return 1;
  }
}


    return 0;
}


