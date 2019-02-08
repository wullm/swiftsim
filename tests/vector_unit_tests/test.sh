#!/bin/bash
let fails=0
for i in ./*.c; do
  gcc -mavx2 $i -lm
  if ./a.out; then
    echo $i success
  else
    echo $i failed
    let fails=$fails+1
  fi
done

echo "$fails tests failed."
