#!/bin/bash

#cd ../../../
#configure --disable-vec --enable-no-gravity-below-id=2
#make -j 8
#cd -
python makeIC_self_grav.py 3.0 0.0
../../swift -G -e -t 1 point_mass_self_grav.yml
python write_position_output.py 3.0 0.0