#include "vector_tests.h"

int main(){

  mask_t mask;
  vec_init_mask_true(mask);

  for(int i = 0; i < VEC_SIZE; i++){
    if(!vec_is_mask_bit_true(mask, i)){
         return 1;
    }
  }

    return 0;
}
