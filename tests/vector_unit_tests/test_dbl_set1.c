#include "vector_tests.h"

int main(){
vector v;

v.vd = vec_dbl_set1(2.0);
for(int i = 0; i < VEC_SIZE/2; i++){
  if(v.d[i] != 2.0){
    return 1;
  }
}



    return 0;
}


