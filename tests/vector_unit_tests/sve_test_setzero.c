#include "vector_tests.h"

int main(){

  VEC_FLOAT z;
  z = vec_setzero();

  vector_array v;
  vec_store(z,v.f);
  //v.v = vec_setzero();
  for(int i = 0; i < VEC_SIZE; i++){
    if(v.f[0] != 0.0f){
      printf("failed\n");
      return 1;
    }
  }

  printf("passed\n");
  return 0;
}


