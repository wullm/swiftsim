#include "vector_tests.h"

int main(){

#ifdef __ARM_NEON
  printf("no vec_set test yet for ARM\n");
  return 0;
#else

vector v;
#if VEC_SIZE == 16
v.v = vec_set(0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0);
#endif

#if VEC_SIZE == 8
v.v = vec_set(0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0);
#endif

#if VEC_SIZE == 4
v.v = vec_set(0.0,1.0,2.0,3.0);
#endif

for(int i = 0; i < VEC_SIZE; i++){
  if(v.f[i] != (float)i){
    return 1;
  }
}
    return 0;

#endif
}
