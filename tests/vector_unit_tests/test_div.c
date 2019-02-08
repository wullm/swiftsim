#include "vector_tests.h"

int main(){

vector v1, v2;
for(int i = 0; i < VEC_SIZE; i++){
  v1.f[i] = (float)(2*(i+1));
  v2.f[i] = (float)(i+1);
}
vector z;
z.v = vec_div(v1.v,v2.v);

for(int i = 0; i < VEC_SIZE; i++){
  if(z.f[i] != 2.0f){
    return 1;
  }
}


    return 0;
}


