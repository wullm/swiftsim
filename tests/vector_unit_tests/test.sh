#!/bin/bash
let fails=0
for i in ./*.c; do
 if  gcc -march=native -mtune=native $i -lm; then
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
