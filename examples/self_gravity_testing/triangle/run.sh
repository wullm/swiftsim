#!/bin/bash

#m_orders=(1 2 3 4 5)
minus_log_eta=(0 1 2 3 4 5)
#cd ../../../
#configure --disable-vec --enable-no-gravity-below-id=11
#make -j 8
#cd -
for i in "${minus_log_eta[@]}"; do 
    cd minus_log_eta_$i
    cp ../makeIC.py ./
    cp ../make_energy_output.py ./
    python makeIC.py
    ../../../swift -G -e -t 1 triangle.yml
    python make_energy_output.py 
    rm *.hdf5
    cd ..
done
python eta_test_plot.py
echo "Done"



