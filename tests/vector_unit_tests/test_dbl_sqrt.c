#include "vector_tests.h"

int main(){

vector v;
for(int i = 0; i < VEC_SIZE/2; i++){
  v.d[i] = (double)((i+1)*(i+1));
}
vector result;
result.vd = vec_dbl_sqrt(v.vd);

for(int i = 0; i < VEC_SIZE/2; i++){
  if(result.d[i] != sqrt((float)((i+1)*(i+1)))){
    return 1;
  }
}


    return 0;
}


