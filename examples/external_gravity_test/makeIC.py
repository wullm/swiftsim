###############################################################################
 # This file is part of SWIFT.
 # Copyright (c) 2016 Stefan Arridge (stefan.arridge@durham.ac.uk)
 # 
 # This program is free software: you can redistribute it and/or modify
 # it under the terms of the GNU Lesser General Public License as published
 # by the Free Software Foundation, either version 3 of the License, or
 # (at your option) any later version.
 # 
 # This program is distributed in the hope that it will be useful,
 # but WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 # GNU General Public License for more details.
 # 
 # You should have received a copy of the GNU Lesser General Public License
 # along with this program.  If not, see <http://www.gnu.org/licenses/>.
 # 
 ##############################################################################

import h5py
import sys
import numpy as np
import math

#constant
CONST_G_CGS = 6.672e-8
# units
const_unit_length_in_cgs = 1.0
const_unit_mass_in_cgs = 1.0
const_unit_velocity_in_cgs = 2.583215e-4

print "UnitMass_in_cgs:     ", const_unit_mass_in_cgs 
print "UnitLength_in_cgs:   ", const_unit_length_in_cgs
print "UnitVelocity_in_cgs: ", const_unit_velocity_in_cgs

#derived quantities
const_unit_time_in_cgs = (const_unit_length_in_cgs / const_unit_velocity_in_cgs)
print "UnitTime_in_cgs:     ", const_unit_time_in_cgs
const_G                = ((CONST_G_CGS*const_unit_mass_in_cgs*const_unit_time_in_cgs*const_unit_time_in_cgs/(const_unit_length_in_cgs*const_unit_length_in_cgs*const_unit_length_in_cgs)))
print 'G=', const_G

# read list of eccentricities
e = np.array(sys.argv)
e = e[1:]
e = e.astype(float)


# Parameters
periodic= 1            # 1 For periodic box
boxSize = 10.          

# we will have one particle on each orbit
n = len(e)

# initally put the particles at apoapsis on the x-axis, semi-major axis is 1.0
coords = np.zeros((n,3))
coords[:,0] = 1. + e

# give them suitable velocity in positive y-direction, so that they stay on these orbits.
# we use units so that G = 1.0, and set M and a (mass and semi-major axis) to be 1.0
vels = np.zeros((n,3))
vels[:,1] = np.sqrt((1.-e**2)/(1.+e)**2)

# Create the file
filename = "external_gravity_test.hdf5"
file = h5py.File(filename, 'w')

#Units
grp = file.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = const_unit_length_in_cgs
grp.attrs["Unit mass in cgs (U_M)"] = const_unit_mass_in_cgs 
grp.attrs["Unit time in cgs (U_t)"] = const_unit_length_in_cgs / const_unit_velocity_in_cgs
grp.attrs["Unit current in cgs (U_I)"] = 1.
grp.attrs["Unit temperature in cgs (U_T)"] = 1.

# Runtime parameters
grp = file.create_group("/RuntimePars")
grp.attrs["PeriodicBoundariesOn"] = periodic

coords += np.full((n,3),boxSize/2.)

# Header
grp = file.create_group("/Header")
grp.attrs["BoxSize"] = boxSize
grp.attrs["NumPart_Total"] =  [0 ,0 , 0, 0, n, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [0, 0, 0, 0, n, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = [0, 0, 0, 0, 0, 0]
grp.attrs["Dimension"] = 3

# Particle group
grp = file.create_group("/PartType4")

ds = grp.create_dataset('Coordinates', (n, 3), 'd')
ds[()] = coords
coords = np.zeros(1)

ds = grp.create_dataset('Velocities', (n, 3), 'f')
ds[()] = vels
vels = np.zeros(1)

# All particles of equal mass
m = np.full((n,),1.0)
ds = grp.create_dataset('Masses', (n, ), 'f')
ds[()] = m
m = np.zeros(1)



# Particle IDs
ids = 1 + np.linspace(0, n, n, endpoint=False)
ds = grp.create_dataset('ParticleIDs', (n, ), 'L')
ds[()] = ids

file.close()
