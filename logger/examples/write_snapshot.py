#!/usr/bin/env python3
"""
Read a logger file by using an index file and then write a snapshot.
Example: ./write_snapshot.py ../../examples/SedovBlast_3D/index 0.1 output.hdf5
"""
import sys
import numpy as np
sys.path.append("../.libs/")

import liblogger as logger

# Get filenames
if len(sys.argv) != 4:
    print("WARNING missing arguments. Will use the default ones")
    basename = "../../examples/HydroTests/SedovBlast_3D/index"
    time = 0.05
    output = "output.hdf5"
else:
    basename = sys.argv[-3]
    time = sys.argv[-2]
    output = sys.argv[-3]

# read dump
data = logger.writeSnapshot(basename, time, output)
