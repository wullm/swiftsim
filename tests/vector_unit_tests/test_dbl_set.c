#include "vector_tests.h"

int main(){

vector v;
#if VEC_SIZE==16
v.vd = vec_dbl_set(0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0);
#endif

#if VEC_SIZE == 8
v.vd = vec_dbl_set(0.0, 1.0, 2.0, 3.0);
#endif

#if VEC_SIZE == 4
v.vd = vec_dbl_set(0.0,1.0);
#endif

for(int i = 0; i < VEC_SIZE/2; i++){
  if(v.d[i] != (double)i){
    return 1;
  }
}
    return 0;
}


