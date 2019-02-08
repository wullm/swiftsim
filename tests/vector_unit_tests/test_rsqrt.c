#include "vector_tests.h"

int main(){

vector v;
for(int i = 0; i < VEC_SIZE; i++){
  v.f[i] = 1.0/(float)((i+1)*(i+1));
}
vector result;
result.v = vec_rsqrt(v.v);

for(int i = 0; i < VEC_SIZE; i++){
  if(result.f[i] < sqrt((float)((i+1)*(i+1)))-1.0 || result.f[i] > sqrt((float)((i+1)*(i+1)) + 1.0)){
    return 1;
  }
}


    return 0;
}


