################################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
#               2017 Bert Vandenbroucke (bv7@st-andrews.ac.uk)
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
################################################################################

import h5py
from numpy import *
import sys
import setups

# Generates a swift IC file for the Gresho-Chan vortex in a periodic box

# Parameters
gamma = 5./3.     # Gas adiabatic index
rho0 = 1           # Gas density
P0 = 0.           # Constant additional pressure (should have no impact on the dynamics)
fileOutputName = "greshoVortex.hdf5" 

if len(sys.argv) < 2:
  print "No number of particles specified!"
  exit()

numPart_request = int(sys.argv[1])

if not numPart_request in setups.setups:
  print "No setup found for requested number of particles!"
  exit()

num_copy = setups.setups[numPart_request]["num_copy"]
fileGlass = "glassPlane_{0}.hdf5".format(
              setups.setups[numPart_request]["glass_name"])

#---------------------------------------------------

# Get position and smoothing lengths from the glass
fileInput = h5py.File(fileGlass, 'r')
pos = fileInput["/PartType0/Coordinates"][:,:]
h = fileInput["/PartType0/SmoothingLength"][:]
boxSize = fileInput["/Header"].attrs["BoxSize"][0]
fileInput.close()

# Make copies
pos /= num_copy
h /= num_copy

# x direction
pos_all = pos
h_all = h
for i in range(1, num_copy):
  pos_all = append(pos_all, pos + array([i * 1. / num_copy, 0., 0.]), axis=0)
  h_all = append(h_all, h)

# y direction
pos = pos_all
h = h_all
for i in range(1, num_copy):
  pos_all = append(pos_all, pos + array([0., i * 1. / num_copy, 0.]), axis=0)
  h_all = append(h_all, h)

pos = pos_all
h = h_all

numPart = size(h)

# Now generate the rest
m = ones(numPart) * rho0 * boxSize**2 / numPart
ids = linspace(1, numPart, numPart)
u = zeros(numPart)
v = zeros((numPart, 3))

for i in range(numPart):
    
    x = pos[i,0]
    y = pos[i,1]

    r2 = (x - boxSize / 2)**2 + (y - boxSize / 2)**2
    r = sqrt(r2)

    v_phi = 0.
    if r < 0.2:
        v_phi = 5.*r
    elif r < 0.4:
        v_phi = 2. - 5.*r
    else:
        v_phi = 0.
    v[i,0] = -v_phi * (y - boxSize / 2) / r
    v[i,1] =  v_phi * (x - boxSize / 2) / r
    v[i,2] = 0.

    P = P0
    if r < 0.2:
        P = P + 5. + 12.5*r2
    elif r < 0.4:
        P = P + 9. + 12.5*r2 - 20.*r + 4.*log(r/0.2)
    else:
        P = P + 3. + 4.*log(2.)
    u[i] = P / ((gamma - 1.)*rho0)
    


#File
fileOutput = h5py.File(fileOutputName, 'w')

# Header
grp = fileOutput.create_group("/Header")
grp.attrs["BoxSize"] = [boxSize, boxSize, 0.2]
grp.attrs["NumPart_Total"] =  [numPart, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFileOutputsPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = [0, 0, 0, 0, 0, 0]
grp.attrs["Dimension"] = 2

#Runtime parameters
grp = fileOutput.create_group("/RuntimePars")
grp.attrs["PeriodicBoundariesOn"] = 1

#Units
grp = fileOutput.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = 1.
grp.attrs["Unit mass in cgs (U_M)"] = 1.
grp.attrs["Unit time in cgs (U_t)"] = 1.
grp.attrs["Unit current in cgs (U_I)"] = 1.
grp.attrs["Unit temperature in cgs (U_T)"] = 1.

#Particle group
grp = fileOutput.create_group("/PartType0")
ds = grp.create_dataset('Coordinates', (numPart, 3), 'd')
ds[()] = pos
ds = grp.create_dataset('Velocities', (numPart, 3), 'f')
ds[()] = v
ds = grp.create_dataset('Masses', (numPart, 1), 'f')
ds[()] = m.reshape((numPart,1))
ds = grp.create_dataset('SmoothingLength', (numPart,1), 'f')
ds[()] = h.reshape((numPart,1))
ds = grp.create_dataset('InternalEnergy', (numPart,1), 'f')
ds[()] = u.reshape((numPart,1))
ds = grp.create_dataset('ParticleIDs', (numPart,1), 'L')
ds[()] = ids.reshape((numPart,1))

fileOutput.close()
