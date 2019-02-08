#include "vector_tests.h"

int main(){

#ifdef HAVE_AVX512_F
printf("Mask tests not available yet for AVX_512F\n");
return 0;
#else


vector v1, v2;
mask_t m;


for(int i = 0; i < VEC_SIZE; i++){
  v1.f[i] = (float)i;
  v2.f[i] = (float)i;
  if(i < VEC_SIZE-1){
      m.i[i] = 0xFFFFFFFF;
  }else{
      m.i[i] = 0;
  }
}
vector z;
z.v = vec_mask_sub(v1.v,v2.v, m);

for(int i = 0; i < VEC_SIZE; i++){
  if(i < VEC_SIZE-1 && z.f[i] != 0.0f){
    return 1;
  }
  if(i == VEC_SIZE-1 && z.f[i] != (float)i){
    return 2;
  }
}


    return 0;
#endif
}


