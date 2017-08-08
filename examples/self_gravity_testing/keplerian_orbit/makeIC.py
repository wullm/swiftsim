###############################################################################
 # This file is part of SWIFT.
 # Copyright (c) 2017 Stefan Arridge (stefan.arridge@durham.ac.uk)
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
CONST_G_CGS =  6.67408e-8
# units
const_unit_length_in_cgs =   1.49597870700e13 # astronomical unit
const_unit_mass_in_cgs = 1.9885e33 # solar mass
const_unit_time_in_cgs =  3.15569252e7 # year
const_unit_velocity_in_cgs = const_unit_length_in_cgs / const_unit_time_in_cgs

const_G = CONST_G_CGS*const_unit_mass_in_cgs*const_unit_time_in_cgs**2/(const_unit_length_in_cgs**3)


# Hard-coded parameters
periodic= 0            # non-periodic box
boxSize = 12.          
particle_mass = 1.0

# Read orbital parameters from command line
sma = float(sys.argv[1])
eccen = float(sys.argv[2])

# Particles are on an orbit determined by the orbital parameters, initially at periapsis
coords = np.zeros((2,3))
coords[1,0] = sma*(1. - eccen)
coords += np.full((3,),boxSize/2.)

vels = np.zeros((2,3))
vels[1,1] = np.sqrt((const_G*particle_mass/sma)*((1.+eccen)/(1.-eccen)))

# set both masses equal to particle_mass, set a few lines above

mass = np.zeros(2)
mass[0] = particle_mass
mass[1] = particle_mass


# Create the file
filename = "point_mass_self_grav.hdf5"
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
