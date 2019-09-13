#include "vector_tests.h"

int main(){


  vector_array v1, v2, v3, result;
  for(int i = 0; i < VEC_SIZE; i++){
    v1.f[i] = (float)i;
    v2.f[i] = (float)i;
    v3.f[i] = (float)i;
  }
  
  VEC_FLOAT z1, z2, z3, z;
  z1 = vec_load(v1.f);
  z2 = vec_load(v2.f);
  z3 = vec_load(v3.f);
  z = vec_fnma(z1,z2,z3);

  vec_store(z, result.f);

  for(int i = 0; i < VEC_SIZE; i++){
    float temp = (float)i;
    if(result.f[i] != -(temp*temp)+temp){
	printf("failed\n");
        return 1;
    }
  }

  printf("passed\n");
  return 0;
}


