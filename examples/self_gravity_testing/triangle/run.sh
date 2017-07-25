#!/bin/bash

#rm plots/xy_error_plot_*
#rm plots/orbit_plot_*
#rm plots/total_error_plot_*
cd ../../../
configure --disable-vec --enable-no-gravity-below-id=4 --with-multipole-order=1
make -j 8
cd -
m_order=(1 2 3 4 5) 
for i in "${m_order[@]}"; do 
    cd ../../../
    configure --disable-vec --enable-no-gravity-below-id=4 --with-multipole-order=$i
    make -j 8
    cd -
    python makeIC.py
    ../../swift -G -e -t 1 triangle.yml
    python write_position_output.py
done




