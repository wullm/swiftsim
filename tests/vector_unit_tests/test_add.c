#include "vector_tests.h"

int main(){

vector v1, v2;
for(int i = 0; i < VEC_SIZE; i++){
  v1.f[i] = (float)i;
  v2.f[i] = (float)i;
}
vector z;
z.v = vec_add(v1.v,v2.v);

for(int i = 0; i < VEC_SIZE; i++){
  if(z.f[i] != (float)i+(float)i){
    return 1;
  }
}


    return 0;
}


