#!/bin/bash -l

module purge
module load swift
module load swift/c5/intel/intelmpi/5.0.3
module load utils

m_order=(1 2 3 4 5) 
for i in "${m_order[@]}"; do 
    cd ../../../
    configure --disable-vec --enable-no-gravity-below-id=3 --with-multipole-order=$i
    make -j 8
    cd -
    python makeIC.py 2 3.0
    ../../swift -G -e -t 1 particle_line_theta_0p7.yml
    python write_position_output.py 2 3.0
    python write_E_L_output.py 2 3.0
done

../../swift -G -t 1 particle_line_theta_0.yml
python write_position_output.py 2 3.0
python write_E_L_output.py 2 3.0

python make_error_plot.py 2 3.0



