#include "vector_tests.h"

int main(){

vector v;
for(int i = 0; i < VEC_SIZE; i++){
  v.i[i] = i;
}
vector z;
z.vd = vec_dbl_load(v.d);
for(int i = 0; i < VEC_SIZE; i++){
   if(z.i[i] != i){
       return 1;
   }
}


    return 0;
}


