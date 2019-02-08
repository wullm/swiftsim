#include "vector_tests.h"

int main(){

vector v;
v.m = vec_setint1(2);

for(int i = 0; i < VEC_SIZE; i++){
  if(v.i[i] != 2){
    return 1;
  }
}


    return 0;
}


