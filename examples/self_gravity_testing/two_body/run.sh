#!/bin/bash

# Generate the initial conditions if they are not present.
echo "Generating initial conditions for the two body system..."
python makeIC.py 1.0 0.0

# Run SWIFT with self gravity
../swift -G -e -t 1 self_gravity_test.yml 2>&1 | tee output.log

# Save a plot of the test particle's orbit
echo "Making a plot of the orbit"
python plot_orbit.py 1.0 0.0
