#include "vector_tests.h"

int main(){

   int max = VEC_SIZE-1;
   float expected = (float)max;
   vector v;
   float result;
   for(int i = 0; i < VEC_SIZE; i++){
    v.f[i] = (float)i;
   }
   VEC_HMAX(v,result);
   if(result != expected){
        return 1;
   }


    return 0;
}


