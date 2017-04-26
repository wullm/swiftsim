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
const_unit_length_in_cgs =  1.495978707e13 # astronomical units
const_unit_mass_in_cgs = 1.98855e33 # solar mass
const_unit_velocity_in_cgs = 4.74372e5 # astronomical units per year

print "UnitMass_in_cgs:     ", const_unit_mass_in_cgs 
print "UnitLength_in_cgs:   ", const_unit_length_in_cgs
print "UnitVelocity_in_cgs: ", const_unit_velocity_in_cgs

#derived quantities
const_unit_time_in_cgs = (const_unit_length_in_cgs / const_unit_velocity_in_cgs)
print "UnitTime_in_cgs:     ", const_unit_time_in_cgs
const_G                = ((CONST_G_CGS*const_unit_mass_in_cgs*const_unit_time_in_cgs*const_unit_time_in_cgs/(const_unit_length_in_cgs*const_unit_length_in_cgs*const_unit_length_in_cgs)))
print 'G=', const_G


# Parameters
periodic= 0            # 1 For periodic box
boxSize = 12.          


# distance between particles is l, set by command line argument
l = float(sys.argv[1])
drift_speed = float(sys.argv[2])
coords = np.zeros((2,3))
coords[0,:] = [1.0,boxSize/2.,boxSize/2.] 
coords[1,:] = [1.0 + l,boxSize/2.,boxSize/2.] 

# particles are on a circular orbit
vels = np.zeros((2,3))
speed = np.sqrt(const_G/(2.*l))
vels[0,:] = [drift_speed,-speed,0]
vels[1,:] = [drift_speed,speed,0]


# both particles have unit mass
mass = np.zeros(2)
mass[0] = 1.0
mass[1] = 1.0
# Create the file
filename = "self_gravity_test.hdf5"
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


# Header
grp = file.create_group("/Header")
grp.attrs["BoxSize"] = boxSize
grp.attrs["NumPart_Total"] =  [0 , 2, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [0, 2, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = [0, 0, 0, 0, 0, 0]
grp.attrs["Dimension"] = 3

# Particle group
grp = file.create_group("/PartType1")

ds = grp.create_dataset('Coordinates', (2, 3), 'd')
ds[()] = coords
coords = np.zeros(1)

ds = grp.create_dataset('Velocities', (2, 3), 'f')
ds[()] = vels
vels = np.zeros(1)

ds = grp.create_dataset('Masses', (2, ), 'f')
ds[()] = mass
mass = np.zeros(1)


# Particle IDs
ids = 1 + np.linspace(0, 2, 2, endpoint=False)
ds = grp.create_dataset('ParticleIDs', (2, ), 'L')
ds[()] = ids

file.close()
