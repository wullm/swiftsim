#!/bin/bash
let fails=0
for i in ./*.c; do
 if  gcc -mavx2 $i -lm; then
  if ./a.out; then
     :
#    echo $i success
  else
    echo $i failed
    let fails=$fails+1
  fi
else
  let fails=$fails+1
fi
done

echo "$fails tests failed."
