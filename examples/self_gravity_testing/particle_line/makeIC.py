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
const_G                = CONST_G_CGS*const_unit_mass_in_cgs*const_unit_time_in_cgs**2/const_unit_length_in_cgs**3
print 'G=', const_G


# Parameters
periodic= 0           # 1 For periodic box
boxSize = 12.          

# number of particles in the line
n = int(sys.argv[1])

#radius of orbit
r = float(sys.argv[2])

# place particles at evenly spaced intervals on the x-axis between 0.25 and 0.75
box_centre = np.full((n+1,3),boxSize/2.)
coords = np.zeros((n+1,3))
if (n > 1):
    coords[:n,0] = np.linspace(0.25,0.75,n)

# if just one particle, have it at 0.5
else:
    coords[0,0] = 0.5

# put the orbiting particle on the x-axis at r + 0.5
coords[n,0] = r + 0.5

# move to centre of box
coords += box_centre

print "Coordinates of particles are: \n", coords

# the mass of the line of particles totals 1
masses = np.zeros(n+1)
masses[:n] = 1./n

# mass of the orbiting particle is 1.0e-6
masses[n] = 1.0e-6

print "Masses of particles are: \n" , masses
# line of particles are stationary, orbiting particle has velocity equal to that of a particle in a circular orbit of radius 1 around 
# unit mass at the origin, in the positive y-direction

vels = np.zeros((n+1,3))
speed = np.sqrt(const_G/r)
vels[n,1] = speed

print "Velocities of particles are: \n" , vels
# Create the file
filename = "particle_line.hdf5"
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
grp.attrs["NumPart_Total"] =  [0 , n+1, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [0, n+1, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = [0, 0, 0, 0, 0, 0]
grp.attrs["Dimension"] = 3

# Particle group
grp = file.create_group("/PartType1")

ds = grp.create_dataset('Coordinates', (n+1, 3), 'd')
ds[()] = coords
coords = np.zeros(1)

ds = grp.create_dataset('Velocities', (n+1, 3), 'f')
ds[()] = vels
vels = np.zeros(1)

ds = grp.create_dataset('Masses', (n+1, ), 'f')
ds[()] = masses
masses = np.zeros(1)


# Particle IDs
ids = 1 + np.linspace(0, n+1, n+1, endpoint=False)
ds = grp.create_dataset('ParticleIDs', (n+1, ), 'L')
ds[()] = ids

file.close()
