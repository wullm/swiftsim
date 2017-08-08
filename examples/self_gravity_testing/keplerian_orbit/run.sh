#!/bin/bash -l 

if [ ! -e ./data/ ]
then
    mkdir data
fi

if [ ! -e ./plots/ ]
then
    mkdir plots
fi

python makeIC_self_grav.py 1.0 0.5

mle_vals=(1 2 3)
for mle in "${mle_vals[@]}"; do
    ../../swift -G -e -t 1 point_mass_self_grav_mle_$mle.yml 2>&1 | tee output.log
    mv timesteps_1.txt timesteps_mle_$mle.txt
    python write_position_output.py 1.0 0.5
    python write_E_L_output.py 1.0 0.5 
done

python make_error_plot.py 1.0 0.5