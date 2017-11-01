#! /bin/bash

for scheme in gadget2 gizmo gadget2 hopkins
do
  echo "Testing " $scheme
  cd ../../..
  source hydro_tests/Noh/1D/"$scheme"_config.sh
  make -j 16
  cd hydro_tests/Noh/1D/
  for n in 100 200 400 800 1600 3200
  do
    folder="$scheme"_"$n"
    mkdir $folder
    python makeIC.py $n
    ../../../examples/swift -s -t 1 noh.yml 2>&1 | tee noh.log
    mv noh_0001.hdf5 $folder/
    mv timesteps_1.txt $folder/
    mv noh.log $folder/
    python analyze.py $folder
  done
done
