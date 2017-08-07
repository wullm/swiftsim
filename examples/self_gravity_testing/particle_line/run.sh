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
    configure --disable-vec --enable-no-gravity-below-id=11 --with-multipole-order=$i
    make -j 8
    cd -
    python makeIC.py 10 3.0
    ../../swift -G -e -t 1 particle_line_theta_1.yml
    python write_position_output.py 3.0
    python write_E_L_output.py 3.0
done

../../swift -G -t 1 particle_line_theta_0.yml
python write_position_output.py 3.0
python write_E_L_output.py 3.0

python make_error_plot.py 3.0



