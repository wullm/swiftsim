#!/bin/bash


#cd ../../../
#configure --disable-vec --enable-no-gravity-below-id=11
#make -j 8
#cd -
#rm plots/xy_error_plot_*
#rm plots/orbit_plot_*
#rm plots/total_error_plot_*

m_order=(1 2 3 4 5) 
for i in "${m_order[@]}"; do 
    cd ../../../
    configure --disable-vec --enable-no-gravity-below-id=3 --with-multipole-order=$i
    make -j 8
    cd -
    python makeIC.py 3 3.0
    ../../swift -G -e -t 1 particle_line.yml
    python write_position_output.py
done




