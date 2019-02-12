#include "vector_tests.h"

int main(){

  vector v = FILL_VEC(0.5f);

  for(int i = 0; i < VEC_SIZE; i++){
    if(v.f[i] != 0.5f){
      return 1;
    }
  }


    return 0;
}


