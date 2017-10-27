#! /bin/bash

if [ ! -e glassCube_64.hdf5 ]
then
  echo "Fetching initial glass files for the 3D disc patch tests..."
  ./getGlass.sh
fi

for scheme in gizmo gadget2 hopkins
do
  echo "Testing " $scheme
  # initial conditions
  cd ../../..
  source hydro_tests/Disc_Patch/3D/"$scheme"_config_ic.sh
  make -j 16
  cd hydro_tests/Disc_Patch/3D/
  for n in $(python setups.py)
  do
    folder="$scheme"_"$n"
    mkdir $folder
    python makeIC.py $n
    ../../../examples/swift -g -s -t 16 disc-patch-ic.yml 2>&1 \
      | tee disc-patch-ic.log
    mv Disc-Patch_0001.hdf5 $folder/Disc-Patch_ic.hdf5
    mv timesteps_16.txt $folder/timesteps_16_ic.txt
    mv disc-patch-ic.log $folder/
    mv
  done
  # stability runs
  cd ../../..
  source hydro_tests/Disc_Patch/3D/"$scheme"_config_run.sh
  make -j 16
  cd hydro_tests/Disc_Patch/3D/
  for n in $(python setups.py)
  do
    folder="$scheme"_"$n"
    cp $folder/Disc-Patch_ic.hdf5 Disc-Patch.hdf5
    ../../../examples/swift -g -s -t 16 disc-patch-run.yml 2>&1 \
      | tee disc-patch-run.log
    mv Disc-Patch_0001.hdf5 $folder/Disc-Patch_run.hdf5
    mv timesteps_16.txt $folder/timesteps_16_run.txt
    mv disc-patch-run.log $folder/
    python analyze.py $folder
  done
done
