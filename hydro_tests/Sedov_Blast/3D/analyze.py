################################################################################
# This file is part of SWIFT.
# Copyright (c) 2015 Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
#                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

# Computes the analytical solution of the 2D Sedov blast wave.
# The script works for a given initial box and dumped energy and computes the
# solution at a later time t.

# Parameters
rho_0 = 1.          # Background Density
P_0 = 1.e-6         # Background Pressure
E_0 = 1.            # Energy of the explosion
gas_gamma = 5./3.   # Gas polytropic index

# ---------------------------------------------------------------
# Don't touch anything after this.
# ---------------------------------------------------------------

import numpy as np
import h5py
import sys

folder = sys.argv[1]
snap = 1

# Read the simulation data
sim = h5py.File("{folder}/sedov_{snap:04d}.hdf5".format(
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
z = np.array(pos[:,2]) - 0.5 * boxSize
vel = sim["/PartType0/Velocities"][:,:]
r = np.sqrt(x**2 + y**2 + z**2)
v_r = (x * vel[:,0] + y * vel[:,1] + z * vel[:,2]) / r
u = sim["/PartType0/InternalEnergy"][:]
S = sim["/PartType0/Entropy"][:]
P = sim["/PartType0/Pressure"][:]
rho = sim["/PartType0/Density"][:]

npart = len(rho)

# Now, work out the solution....

from scipy.special import gamma as Gamma

def calc_a(g,nu=3):
    """ 
    exponents of the polynomials of the sedov solution
    g - the polytropic gamma
    nu - the dimension
    """
    a = [0]*8
   
    a[0] = 2.0 / (nu + 2)
    a[2] = (1-g) / (2*(g-1) + nu)
    a[3] = nu / (2*(g-1) + nu)
    a[5] = 2 / (g-2)
    a[6] = g / (2*(g-1) + nu)
   
    a[1] = ( ((nu+2)*g)/(2.0+nu*(g-1.0)) ) * \
           ( (2.0*nu*(2.0-g))/(g*(nu+2.0)**2) - a[2])
    a[4] = a[1]*(nu+2) / (2-g)
    a[7] = (2 + nu*(g-1))*a[1]/(nu*(2-g))
    return a

def calc_beta(v, g, nu=3):
    """ 
    beta values for the sedov solution (coefficients of the polynomials of the
    similarity variables) 
    v - the similarity variable
    g - the polytropic gamma
    nu- the dimension
    """

    beta = (nu+2) * (g+1) * np.array((0.25, (g/(g-1))*0.5,
            -(2 + nu*(g-1))/2.0 / ((nu+2)*(g+1) -2*(2 + nu*(g-1))),
     -0.5/(g-1)), dtype=np.float64)

    beta = np.outer(beta, v)

    beta += (g+1) * np.array((0.0,  -1.0/(g-1),
                           (nu+2) / ((nu+2)*(g+1) -2.0*(2 + nu*(g-1))),
                           1.0/(g-1)), dtype=np.float64).reshape((4,1))

    return beta


def sedov(t, E0, rho0, g, n=1000, nu=3):
    """ 
    solve the sedov problem
    t - the time
    E0 - the initial energy
    rho0 - the initial density
    n - number of points (10000)
    nu - the dimension
    g - the polytropic gas gamma
    """
    # the similarity variable
    v_min = 2.0 / ((nu + 2) * g)
    v_max = 4.0 / ((nu + 2) * (g + 1))

    v = v_min + np.arange(n) * (v_max - v_min) / (n - 1.0)

    a = calc_a(g, nu)
    beta = calc_beta(v, g=g, nu=nu)
    lbeta = np.log(beta)
    
    r = np.exp(-a[0] * lbeta[0] - a[2] * lbeta[1] - a[1] * lbeta[2])
    rho = ((g + 1.0) / (g - 1.0)) * \
          np.exp(a[3] * lbeta[1] + a[5] * lbeta[3] + a[4] * lbeta[2])
    p = np.exp(nu * a[0] * lbeta[0] + (a[5] + 1) * lbeta[3] + \
        (a[4] - 2 * a[1]) * lbeta[2])
    u = beta[0] * r * 4.0 / ((g + 1) * (nu + 2))
    p *= 8.0 / ((g + 1) * (nu + 2) * (nu + 2))

    # we have to take extra care at v=v_min, since this can be a special point.
    # It is not a singularity, however, the gradients of our variables (wrt v) are.
    # r -> 0, u -> 0, rho -> 0, p-> constant

    u[0] = 0.0; rho[0] = 0.0; r[0] = 0.0; p[0] = p[1]

    # volume of an n-sphere
    vol = (np.pi ** (nu / 2.0) / Gamma(nu / 2.0 + 1)) * r**nu

    # note we choose to evaluate the integral in this way because the
    # volumes of the first few elements (i.e near v=vmin) are shrinking 
    # very slowly, so we dramatically improve the error convergence by 
    # finding the volumes exactly. This is most important for the
    # pressure integral, as this is on the order of the volume.

    # (dimensionless) energy of the model solution
    de = rho * u * u * 0.5 + p / (g - 1)
    # integrate (trapezium rule)
    q = np.inner(de[1:] + de[:-1], np.diff(vol)) * 0.5

    # the factor to convert to this particular problem
    fac = (q * (t ** nu) * rho0 / E0) ** (-1.0 / (nu + 2))

    # shock speed
    shock_speed = fac * (2.0 / (nu + 2))
    rho_s = ((g + 1) / (g - 1)) * rho0
    r_s = shock_speed * t * (nu + 2) / 2.0
    p_s = (2.0 * rho0 * shock_speed * shock_speed) / (g + 1)
    u_s = (2.0 * shock_speed) / (g + 1)

    r *= fac * t
    u *= fac
    p *= fac * fac * rho0
    rho *= rho0
    return r, p, rho, u, r_s, p_s, rho_s, u_s, shock_speed

# The main properties of the solution
r_s, P_s, rho_s, v_s, r_shock, _, _, _, _ = \
  sedov(time, E_0, rho_0, gas_gamma, 1000, 3)

# Append points for after the shock
r_s = np.insert(r_s, np.size(r_s), [r_shock, r_shock*1.5])
rho_s = np.insert(rho_s, np.size(rho_s), [rho_0, rho_0])
P_s = np.insert(P_s, np.size(P_s), [P_0, P_0])
v_s = np.insert(v_s, np.size(v_s), [0, 0])

# filter out duplicate values in the solution arrays
r_s, u_i = np.unique(r_s, return_index = True)
rho_s = np.array([rho_s[i] for i in u_i])
v_s = np.array([v_s[i] for i in u_i])
P_s = np.array([P_s[i] for i in u_i])
# make sure the edge of the shock front is set to the unshocked values
rho_s[-2] = rho_0
v_s[-2] = 0.
P_s[-2] = P_0

# interpolate through the analytic solution
import scipy.interpolate as interpol
rhospline = interpol.interp1d(r_s, rho_s, bounds_error=False, fill_value=rho_0)
vspline = interpol.interp1d(r_s, v_s, bounds_error=False, fill_value=0.)
Pspline = interpol.interp1d(r_s, P_s, bounds_error=False, fill_value=P_0)

N = len(r)
rho_s = rhospline(r)
v_s = vspline(r)
P_s = Pspline(r)

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

times = np.loadtxt("{folder}/timesteps_16.txt".format(folder = folder))

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

fig, ax = pl.subplots(2, 3, sharex = True)
ax[0][0].plot(r, rho, "k.")
ax[0][0].plot(r_s, rhospline(r_s), "r-")
ax[0][1].plot(r, v_r, "k.")
ax[0][1].plot(r_s, vspline(r_s), "r-")
ax[0][2].plot(r, P, "k.")
ax[0][2].plot(r_s, Pspline(r_s), "r-")
ax[1][0].plot(r, rho_xi2_tot_array, "k.")
ax[1][1].plot(r, v_xi2_tot_array, "k.")
ax[1][2].plot(r, P_xi2_tot_array, "k.")
ax[0][0].set_title("density")
ax[0][1].set_title("velocity")
ax[0][2].set_title("pressure")
pl.suptitle("{0}, {1} particles".format(scheme, npart))

pl.savefig("{folder}/result.png".format(folder = folder))
