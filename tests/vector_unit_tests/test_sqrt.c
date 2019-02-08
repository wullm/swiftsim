#include "vector_tests.h"

int main(){

vector v;
for(int i = 0; i < VEC_SIZE; i++){
  v.f[i] = (float)((i+1)*(i+1));
}
vector result;
result.v = vec_sqrt(v.v);

for(int i = 0; i < VEC_SIZE; i++){
  if(result.f[i] != sqrt((float)((i+1)*(i+1)))){
    return 1;
  }
}


    return 0;
}


