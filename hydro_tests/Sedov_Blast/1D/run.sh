#! /bin/bash

for scheme in gizmo gadget2 hopkins
do
  echo "Testing " $scheme
  cd ../../..
  source hydro_tests/Sedov_Blast/1D/"$scheme"_config.sh
  make -j 16
  cd hydro_tests/Sedov_Blast/1D/
  for n in 100 200 400 800 1600 3200
  do
    folder="$scheme"_"$n"
    mkdir $folder
    python makeIC.py $n
    ../../../examples/swift -s -t 1 sedov.yml
    mv sedov_0005.hdf5 $folder/
    mv timesteps_1.txt $folder/
    python analyze.py $folder
  done
done
