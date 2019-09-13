#include "vector_tests.h"

int main(){

  vector_array v;
  for(int i = 0; i < VEC_SIZE; i++){
    v.i[i] = i;
  }

  //VEC_FLOAT x;
  //x = vec_load(v.f);
  VEC_DBL y;
  y = vec_load(v.d);
  
  vector_array z;
  //vec_store(x,z.f);
  vec_store(y,z.d);

  for(int i = 0; i < VEC_SIZE; i++){
     if(z.i[i] != i){
       printf("failed\n");
       return 1;
     }
  }

  printf("passed\n");
  return 0;
}


