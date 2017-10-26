#! /bin/bash

if [ ! -e glassPlane_128.hdf5 ]
then
  echo "Fetching initial glass files for the 2D Sedov blast tests..."
  ./getGlass.sh
fi

for scheme in gizmo gadget2 hopkins
do
  echo "Testing " $scheme
  cd ../../..
  source hydro_tests/Sedov_Blast/2D/"$scheme"_config.sh
  make -j 16
  cd hydro_tests/Sedov_Blast/2D/
  for n in $(python setups.py)
  do
    folder="$scheme"_"$n"
    mkdir $folder
    python makeIC.py $n
    ../../../examples/swift -s -t 4 sedov.yml
    mv sedov_0005.hdf5 $folder/
    mv timesteps_4.txt $folder/
    python analyze.py $folder
  done
done
