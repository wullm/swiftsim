#include "vector_tests.h"

int main(){

  vector v1,v2; 
  for(int i = 0; i < VEC_SIZE; i++){
    v1.f[i] = (float)i;
    if(i & 1){
      v2.f[i] = (float)i;
    }else{
      v2.f[i] = -1.0;
    }
  }
  int result = vec_cmp_result(vec_cmp_lte(v1.v, v2.v));
  for(int i = 0; i < VEC_SIZE; i++){
    if((i & 1) && !(result & (1 << i)) ){
      return 1;
    }else if( !(i&1) && result & (1<<i)){
      return 1;
    }
  }

  return 0;
}
