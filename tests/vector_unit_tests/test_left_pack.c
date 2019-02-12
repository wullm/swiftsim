#include "vector_tests.h"

int main(){

  int mask = 0xAAAA;
  mask_t packed_mask;
  vector v1;
  float result[VEC_SIZE];

  for(int i = 0; i < VEC_SIZE; i++){
    v1.f[i] = (float)i;
    result[i] = 0.f;
  }

  VEC_FORM_PACKED_MASK(mask, packed_mask);
  VEC_LEFT_PACK(v1.v, packed_mask, result);

  for(int i = 0; i < VEC_SIZE; i++){
    if(i < VEC_SIZE / 2 && result[i] != 2*i + 1.0){
      return 1;
    }else if(i >= VEC_SIZE / 2 && result[i] != 0){
      return 2;
    }
  }

  return 0;
} 
