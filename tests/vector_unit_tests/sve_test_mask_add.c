#include "vector_tests.h"

int main(){

  vector_array v1, v2;
  mask_t mask;

  for(int i = 0; i < VEC_SIZE; i++){
    v1.f[i] = (float)i;
    v2.f[i] = (float)i;
  }

  for (int i = 0; i < VEC_SIZE; i+=svcntw()) mask = svwhilelt_b32(i, VEC_SIZE-1);
  
  VEC_FLOAT z1, z2, z3;
  z1 = vec_load(v1.f);
  z2 = vec_load(v2.f);
  z3 = vec_mask_add(z1,z2,mask);
  
  vector_array z;
  vec_store(z3, z.f);
  for (int i = 0; i < VEC_SIZE; i++) printf("%d %f\n", i, z.f[i]);

  for(int i = 0; i < VEC_SIZE; i++){
    if(i < VEC_SIZE-1 && z.f[i] != (float)i+(float)i){
      printf("failed 1\n");
      return 1;
    }
    if(i == VEC_SIZE-1 && z.f[i] != (float)i){
      printf("failed 2 i %d z.f[i] %f\n", i, z.f[i]);
      return 2;
    }
  }

  printf("passed\n");
  return 0;
}
