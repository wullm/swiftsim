#! /bin/bash

if [ ! -e glassCube_64.hdf5 ]
then
  echo "Fetching initial glass files for the 3D Noh tests..."
  ./getGlass.sh
fi

for scheme in gizmo gadget2 hopkins
do
  echo "Testing " $scheme
  cd ../../..
  source hydro_tests/Noh/3D/"$scheme"_config.sh
  make -j 16
  cd hydro_tests/Noh/3D/
  for n in $(python setups.py)
  do
    folder="$scheme"_"$n"
    mkdir $folder
    python makeIC.py $n
    ../../../examples/swift -s -t 16 noh.yml 2>&1 | tee noh.log
    mv noh_0001.hdf5 $folder/
    mv timesteps_16.txt $folder/
    mv noh.log $folder/
    python analyze.py $folder
  done
done
