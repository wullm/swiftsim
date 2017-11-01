################################################################################
# This file is part of SWIFT.
# Copyright (c) 2017 Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
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

# Computes the analytical solution of the 2D Noh test

# Parameters
rho0 = 1.          # Background Density
P0 = 1.e-6         # Background Pressure
v0 = 1.            # Initial radial velocity
gas_gamma = 5./3.  # Gas polytropic index

# ---------------------------------------------------------------
# Don't touch anything after this.
# ---------------------------------------------------------------

import numpy as np
import h5py
import sys

folder = sys.argv[1]
snap = 1

# Read the simulation data
sim = h5py.File("{folder}/noh_{snap:04d}.hdf5".format(
                 folder = folder, snap = snap), "r")
boxSize = sim["/Header"].attrs["BoxSize"][0]
time = sim["/Header"].attrs["Time"][0]
scheme = sim["/HydroScheme"].attrs["Scheme"]
kernel = sim["/HydroScheme"].attrs["Kernel function"]
neighbours = sim["/HydroScheme"].attrs["Kernel target N_ngb"]
eta = sim["/HydroScheme"].attrs["Kernel eta"]
git = sim["Code"].attrs["Git Revision"]

pos = sim["/PartType0/Coordinates"][:,:]
x = np.array(pos[:,0]) - 0.5 * boxSize
y = np.array(pos[:,1]) - 0.5 * boxSize
vel = sim["/PartType0/Velocities"][:,:]
r = np.sqrt(x**2 + y**2)
v_r = (x * vel[:,0] + y * vel[:,1]) / r
u = sim["/PartType0/InternalEnergy"][:]
S = sim["/PartType0/Entropy"][:]
P = sim["/PartType0/Pressure"][:]
rho = sim["/PartType0/Density"][:]

N = len(rho)

# Now, work out the solution....

rs = 0.5 * (gas_gamma - 1.) * v0 * time

rho_s = np.where(r < rs,
                 np.ones(N) * rho0 * ((gas_gamma + 1.) / (gas_gamma - 1.))**2,
                 rho0 * (1. + v0 * time / r))
v_s = np.where(r < rs, np.zeros(N), -np.ones(N) * v0)
P_s = np.where(r < rs,
               np.ones(N) * 0.5 * rho0 * v0**2 * \
                 (gas_gamma + 1.)**2 / (gas_gamma - 1.),
               np.ones(N) * P0)

rho_xi2_tot_array = (rho - rho_s)
rho_xi2_tot = sum( rho_xi2_tot_array**2 ) / N

v_xi2_tot_array = (v_r - v_s)
v_xi2_tot = sum( v_xi2_tot_array**2 ) / N

P_xi2_tot_array = (P - P_s)
P_xi2_tot = sum( P_xi2_tot_array**2 ) / N

print "rho:"
print "xi2 total:", rho_xi2_tot
       
print "v:"
print "xi2 total:", v_xi2_tot
       
print "P:"
print "xi2 total:", P_xi2_tot

times = np.loadtxt("{folder}/timesteps_4.txt".format(folder = folder))

time_total = sum(times[:,6])

print "Total time (ms):", time_total

file = open("{folder}/summary.txt".format(folder = folder), 'w')

file.write("{{\"rho_xi2_tot\": {val},\n".format(val = rho_xi2_tot))

file.write("\"v_xi2_tot\": {val},\n".format(val = v_xi2_tot))

file.write("\"P_xi2_tot\": {val},\n".format(val = P_xi2_tot))

file.write("\"time_total\": {val}}}\n".format(val = time_total))

import matplotlib
matplotlib.use("Agg")
import pylab as pl

pl.rcParams["figure.figsize"] = (10, 8)
pl.rcParams["text.usetex"] = True

rrange = np.arange(0., 0.5 * boxSize, 0.0005 * boxSize)
rhorange = np.where(rrange < rs,
                    np.ones(1000) * rho0 * \
                      ((gas_gamma + 1.) / (gas_gamma - 1.))**2,
                    np.ones(1000) * rho0 * (1. + v0 * time / rrange))
vrange = np.where(rrange < rs, np.zeros(1000),
                  -np.ones(1000) * v0)
Prange = np.where(rrange < rs,
                  np.ones(1000) * 0.5 * rho0 * v0**2 * \
                    (gas_gamma + 1.)**2 / (gas_gamma - 1.),
                  np.ones(1000) * P0)

fig, ax = pl.subplots(2, 3, sharex = True)
ax[0][0].plot(r, rho, "k.")
ax[0][0].plot(rrange, rhorange, "r-")
ax[0][1].plot(r, v_r, "k.")
ax[0][1].plot(rrange, vrange, "r-")
ax[0][2].plot(r, P, "k.")
ax[0][2].plot(rrange, Prange, "r-")
ax[1][0].plot(r, rho_xi2_tot_array, "k.")
ax[1][1].plot(r, v_xi2_tot_array, "k.")
ax[1][2].plot(r, P_xi2_tot_array, "k.")
ax[0][0].set_title("density")
ax[0][1].set_title("velocity")
ax[0][2].set_title("pressure")
ax[0][0].set_ylabel("value")
ax[1][0].set_ylabel("absolute difference")
ax[1][0].set_xlabel("radius")
ax[1][1].set_xlabel("radius")
ax[1][2].set_xlabel("radius")

pl.tight_layout()
pl.subplots_adjust(top = 0.9)
pl.suptitle("{0}, {1} particles".format(scheme, N))
pl.savefig("{folder}/result.png".format(folder = folder))
