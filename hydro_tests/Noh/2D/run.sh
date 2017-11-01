#! /bin/bash

if [ ! -e glassPlane_128.hdf5 ]
then
  echo "Fetching initial glass files for the 2D Noh tests..."
  ./getGlass.sh
fi

for scheme in gadget2 gizmo gadget2 hopkins
do
  echo "Testing " $scheme
  cd ../../..
  source hydro_tests/Noh/2D/"$scheme"_config.sh
  make -j 16
  cd hydro_tests/Noh/2D/
  for n in $(python setups.py)
  do
    folder="$scheme"_"$n"
    mkdir $folder
    python makeIC.py $n
    ../../../examples/swift -s -t 4 noh.yml 2>&1 | tee noh.log
    mv noh_0001.hdf5 $folder/
    mv timesteps_4.txt $folder/
    mv noh.log $folder/
    python analyze.py $folder
  done
done
