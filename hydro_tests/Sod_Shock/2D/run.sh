#! /bin/bash

if [ ! -e glassPlane_128.hdf5 ]
then
  echo "Fetching initial glass files for the 2D Sod convergence test..."
  ./getGlass.sh
fi

for scheme in gizmo gadget2 hopkins
do
  echo "Testing " $scheme
  cd ../../..
  source hydro_tests/Sod_Shock/2D/"$scheme"_config.sh
  make -j 16
  cd hydro_tests/Sod_Shock/2D/
  for n in $(python setups.py)
  do
    folder="$scheme"_"$n"
    mkdir $folder
    python makeIC.py $n
    ../../../examples/swift -s -t 4 sodShock.yml 2>&1 | tee sodShock.log
    mv sodShock_0001.hdf5 $folder/
    mv timesteps_4.txt $folder/
    mv sodShock.log $folder/
    python analyze.py $folder
  done
done
