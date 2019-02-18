#!/bin/bash
let fails=0
#for i in ./*.c; do
for i in test_cmp_lte.c; do
  gcc -march=native -mcpu=native $i -lm
  if ./a.out; then
     :
#    echo $i success
  else
    echo $i failed
    let fails=$fails+1
  fi
else
  echo $i failed to compile
  let fails=$fails+1
fi
done

echo "$fails tests failed."
