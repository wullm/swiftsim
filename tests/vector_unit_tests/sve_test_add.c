#include "vector_tests.h"

int main(){

  vector_array v1, v2;
  for(int i = 0; i < VEC_SIZE; i++){
    v1.f[i] = (float)i;
    v2.f[i] = (float)i;
  }
  
  VEC_FLOAT z1, z2, z3;
  z1 = vec_load(v1.f);
  z2 = vec_load(v2.f);
  z3 = vec_add(z1,z2);
  
  vector_array z;
  vec_store(z3, z.f);
  
  for(int i = 0; i < VEC_SIZE; i++){
    if(z.f[i] != (float)i+(float)i){
      printf("failed\n");
      return 1;
    }
  }

  printf("passed\n");
  return 0;
}


