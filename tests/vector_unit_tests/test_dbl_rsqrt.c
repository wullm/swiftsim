#include "vector_tests.h"

int main(){

vector v;
for(int i = 0; i < VEC_SIZE; i++){
  v.d[i] = 1.0/(double)((i+1)*(i+1));
}
vector result;
result.vd = vec_dbl_rsqrt(v.vd);

for(int i = 0; i < VEC_SIZE; i++){
  if(result.d[i] < sqrt((double)((i+1)*(i+1)))-1.0 || result.d[i] > sqrt((double)((i+1)*(i+1)) + 1.0)){
    return 1;
  }
}


    return 0;
}


