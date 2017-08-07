#!/bin/bash -l 

#cd ../../../
#configure --disable-vec --enable-no-gravity-below-id=2
#make -j 8
#cd -

#
# Batch script for bash users
#
#BSUB -L /bin/bash
#BSUB -n 1
#BSUB -J swift_point_mass_const_dt
#BSUB -oo ./out.%J
#BSUB -eo ./err.%J
#BSUB -q cosma5
#BSUB -P dp004
#BSUB -W 00:20


module purge
module load intel_comp/c4/2015 intel_mpi/5.0.3 parallel_hdf5/1.8.14-mt metis/5.1.0 doxygen
module load utils

# Run the program

python makeIC_self_grav.py 1.0 0.5

mle_vals=(1 2 3)
for mle in "${mle_vals[@]}"; do
    ../../swift -G -e -t 1 point_mass_self_grav_mle_$mle.yml 2>&1 | tee output.log
    mv timesteps_1.txt timesteps_mle_$mle.txt
    python write_position_output.py 1.0 0.5
    python write_E_L_output.py 1.0 0.5 
done

python make_error_plot.py 1.0 0.5