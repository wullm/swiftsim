#! /bin/bash --login

git clone https://gitlab.cosma.dur.ac.uk/swift/swiftsim.git
mkdir outputs

cd swiftsim
./autogen.sh
# Change this if you would like to use a different branch.
# git checkout master
cd ..

for scheme in gizmo gadget2 hopkins
do
  echo "Testing " $scheme
  cd swiftsim
  make clean
  source ../"$scheme"_config.sh
  make -j 8
  cd ../outputs
  for n in 100 150 200 300 400 500 600 700 800 900 1000 1200 1400 1600 2000 2400 2800 3200
  do
    folder="$scheme"_"$n"
    mkdir $folder
    python3 ../../makeIC.py $n
    ../swiftsim/examples/swift -s -t 16 ../../sodShock.yml 2>&1 | tee sodShock.log
    mv sodShock_0001.hdf5 $folder/
    mv timesteps_*.txt $folder/
    mv sodShock.log $folder/
  done
  cd ..
done

python3 make_convergence_plot.py


