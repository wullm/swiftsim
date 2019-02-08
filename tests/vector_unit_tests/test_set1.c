#include "vector_tests.h"

int main(){
vector v;

v.v = vec_set1(2.0f);
for(int i = 0; i < VEC_SIZE; i++){
  if(v.f[i] != 2.0f){
    return 1;
  }
}



    return 0;
}


