#include "vector_tests.h"

int main(){

   vector v1;
   float expected = 0.0f;
   for(int i = 0; i < VEC_SIZE; i++){
     v1.f[i] = (float)i;
     expected += v1.f[i];
   }

   float result;
   VEC_HADD(v1, result);
  
   if(result != expected){
    return 1;
   }

    return 0;
}


