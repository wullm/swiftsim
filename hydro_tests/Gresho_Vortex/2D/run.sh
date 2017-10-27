#! /bin/bash

if [ ! -e glassPlane_128.hdf5 ]
then
  echo "Fetching initial glass files for the 2D Gresho vortex tests..."
  ./getGlass.sh
fi

for scheme in gizmo gadget2 hopkins
do
  echo "Testing " $scheme
  cd ../../..
  source hydro_tests/Gresho_Vortex/2D/"$scheme"_config.sh
  make -j 16
  cd hydro_tests/Gresho_Vortex/2D/
  for n in $(python setups.py)
  do
    folder="$scheme"_"$n"
    mkdir $folder
    python makeIC.py $n
    ../../../examples/swift -s -t 4 gresho.yml
    mv gresho_0011.hdf5 $folder/
    mv timesteps_4.txt $folder/
    python analyze.py $folder
  done
done
