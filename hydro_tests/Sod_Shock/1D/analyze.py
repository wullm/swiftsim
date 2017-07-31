################################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
#               2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

# Computes the analytic solution and the xi2 difference between the simulation
# result and the analytic solution for the given folder, both for the total
# solution and for the three characteristic waves present in the solution.

# Also computes the total runtime of the simulation

# Generates the analytical  solution for the Sod shock test case
# The script works for a given left (x<0) and right (x>0) state and computes the
# solution at a later time t. This follows the solution given in (Toro, 2009).

# Parameters
gas_gamma = 5./3.      # Polytropic index
rho_L = 1.             # Density left state
rho_R = 0.125          # Density right state
v_L = 0.               # Velocity left state
v_R = 0.               # Velocity right state
P_L = 1.               # Pressure left state
P_R = 0.1              # Pressure right state

# sub-intervals containing the different characteristic waves
rarefaction_interval = [-0.5, 0.]
contact_interval = [0., 0.25]
shock_interval = [0.25, 0.5]

import numpy as np
import h5py
import sys
import matplotlib
#matplotlib.use("Agg")
import pylab as pl

folder = sys.argv[1]
snap = 1

# Read the simulation data
sim = h5py.File("{folder}/sodShock_{snap:04d}.hdf5".format(
                 folder = folder, snap = snap), "r")
boxSize = sim["/Header"].attrs["BoxSize"][0]
time = sim["/Header"].attrs["Time"][0]
scheme = sim["/HydroScheme"].attrs["Scheme"]
kernel = sim["/HydroScheme"].attrs["Kernel function"]
neighbours = sim["/HydroScheme"].attrs["Kernel target N_ngb"]
eta = sim["/HydroScheme"].attrs["Kernel eta"]
git = sim["Code"].attrs["Git Revision"]

x = np.array(sim["/PartType0/Coordinates"])[:,0]
v = sim["/PartType0/Velocities"][:,0]
u = sim["/PartType0/InternalEnergy"][:]
S = sim["/PartType0/Entropy"][:]
P = sim["/PartType0/Pressure"][:]
rho = sim["/PartType0/Density"][:]

N = len(x)  # Number of points
x_min = -1.
x_max = 1.

x += x_min

# selection of characteristic wave sub intervals
idx_rarefaction = np.intersect1d(np.where(x > rarefaction_interval[0])[0],
                                 np.where(x < rarefaction_interval[1])[0])
N_rarefaction = len(idx_rarefaction)
idx_contact = np.intersect1d(np.where(x > contact_interval[0])[0],
                             np.where(x < contact_interval[1])[0])
N_contact = len(idx_contact)
idx_shock = np.intersect1d(np.where(x > shock_interval[0])[0],
                           np.where(x < shock_interval[1])[0])
N_shock = len(idx_shock)

# ---------------------------------------------------------------
# Don't touch anything after this.
# ---------------------------------------------------------------

c_L = np.sqrt(gas_gamma * P_L / rho_L)   # Speed of the rarefaction wave
c_R = np.sqrt(gas_gamma * P_R / rho_R)   # Speed of the shock front

# Helpful variable
Gama = (gas_gamma - 1.) / (gas_gamma + 1.)
beta = (gas_gamma - 1.) / (2. * gas_gamma)

# Characteristic function and its derivative, following Toro (2009)
def compute_f(P_3, P, c):
    u = P_3 / P
    if u > 1:
        term1 = gas_gamma*((gas_gamma+1.)*u + gas_gamma-1.)
        term2 = np.sqrt(2./term1)
        fp = (u - 1.)*c*term2
        dfdp = c*term2/P + \
               (u - 1.)*c/term2*(-1./term1**2)*gas_gamma*(gas_gamma+1.)/P
    else:
        fp = (u**beta - 1.)*(2.*c/(gas_gamma-1.))
        dfdp = 2.*c/(gas_gamma-1.)*beta*u**(beta-1.)/P
    return (fp, dfdp)

# Solution of the Riemann problem following Toro (2009) 
def RiemannProblem(rho_L, P_L, v_L, rho_R, P_R, v_R):
    P_new = ((c_L + c_R + (v_L - v_R)*0.5*(gas_gamma-1.))/ \
             (c_L / P_L**beta + c_R / P_R**beta))**(1./beta)
    P_3 = 0.5*(P_R + P_L)
    f_L = 1.
    while abs(P_3 - P_new) > 1e-6:
        P_3 = P_new
        (f_L, dfdp_L) = compute_f(P_3, P_L, c_L)
        (f_R, dfdp_R) = compute_f(P_3, P_R, c_R)
        f = f_L + f_R + (v_R - v_L)
        df = dfdp_L + dfdp_R
        dp =  -f/df
        prnew = P_3 + dp
    v_3 = v_L - f_L
    return (P_new, v_3)

# Solve Riemann problem for post-shock region
(P_3, v_3) = RiemannProblem(rho_L, P_L, v_L, rho_R, P_R, v_R)

# Check direction of shocks and wave
shock_R = (P_3 > P_R)
shock_L = (P_3 > P_L)

# Velocity of shock front and and rarefaction wave
if shock_R:
    v_right = v_R + c_R**2*(P_3/P_R - 1.)/(gas_gamma*(v_3-v_R))
else:
    v_right = c_R + 0.5*(gas_gamma+1.)*v_3 - 0.5*(gas_gamma-1.)*v_R

if shock_L:
    v_left = v_L + c_L**2*(P_3/p_L - 1.)/(gas_gamma*(v_3-v_L))
else:
    v_left = c_L - 0.5*(gas_gamma+1.)*v_3 + 0.5*(gas_gamma-1.)*v_L

# Compute position of the transitions
x_23 = -abs(v_left) * time
if shock_L :
    x_12 = -abs(v_left) * time
else:
    x_12 = -(c_L - v_L) * time

x_34 = v_3 * time

x_45 = abs(v_right) * time
if shock_R:
    x_56 = abs(v_right) * time
else:
    x_56 = (c_R + v_R) * time

# Prepare arrays
rho_s = np.zeros(N)
P_s = np.zeros(N)
v_s = np.zeros(N)

# Compute solution in the different regions
for i in range(N):
    if x[i] >= -0.5 and x[i] <= 0.5:
        x_s = x[i]
    else:
        if x[i] < -0.5:
            x_s = -1. - x[i]
        else:
            x_s = 1. - x[i]
    if x_s <= x_12:
        rho_s[i] = rho_L
        P_s[i] = P_L
        v_s[i] = v_L
    if x_s >= x_12 and x_s < x_23:
        if shock_L:
            rho_s[i] = rho_L*(Gama + P_3/P_L)/(1. + Gama * P_3/P_L)
            P_s[i] = P_3
            v_s[i] = v_3
        else:
            rho_s[i] = rho_L*(Gama * (0. - x_s)/(c_L * time) + \
                              Gama * v_L/c_L + (1.-Gama))**(2./(gas_gamma-1.))
            P_s[i] = P_L*(rho_s[i] / rho_L)**gas_gamma
            v_s[i] = (1.-Gama)*(c_L -(0. - x_s) / time) + Gama*v_L
    if x_s >= x_23 and x_s < x_34:
        if shock_L:
            rho_s[i] = rho_L*(Gama + P_3/P_L)/(1+Gama * P_3/p_L)
        else:
            rho_s[i] = rho_L*(P_3 / P_L)**(1./gas_gamma)
        P_s[i] = P_3
        v_s[i] = v_3
    if x_s >= x_34 and x_s < x_45:
        if shock_R:
            rho_s[i] = rho_R*(Gama + P_3/P_R)/(1. + Gama * P_3/P_R)
        else:
            rho_s[i] = rho_R*(P_3 / P_R)**(1./gas_gamma)
        P_s[i] = P_3
        v_s[i] = v_3
    if x_s >= x_45 and x_s < x_56:
        if shock_R:
            rho_s[i] = rho_R
            P_s[i] = P_R
            v_s[i] = v_R
        else:
            rho_s[i] = rho_R*(Gama*(x_s)/(c_R*time) - \
                              Gama*v_R/c_R + (1.-Gama))**(2./(gas_gamma-1.))
            P_s[i] = p_R*(rho_s[i]/rho_R)**gas_gamma
            v_s[i] = (1.-Gama)*(-c_R - (-x_s)/time) + Gama*v_R
    if x_s >= x_56:
        rho_s[i] = rho_R
        P_s[i] = P_R
        v_s[i] = v_R
        
    if x[i] < -0.5 or x[i] > 0.5:
        v_s[i] = -v_s[i]

rho_xi2_tot = sum( ( (rho-rho_s)/(rho+rho_s) )**2 ) / N
rho_xi2_rar = \
  sum( ( (rho[idx_rarefaction]-rho_s[idx_rarefaction])/ \
         (rho[idx_rarefaction]+rho_s[idx_rarefaction]) )**2 ) / N_rarefaction
rho_xi2_con = \
  sum( ( (rho[idx_contact]-rho_s[idx_contact])/ \
         (rho[idx_contact]+rho_s[idx_contact]) )**2 ) / N_contact
rho_xi2_sho = \
  sum( ( (rho[idx_shock]-rho_s[idx_shock])/ \
         (rho[idx_shock]+rho_s[idx_shock]) )**2 ) / N_shock

v_xi2_tot = sum( ( (v-v_s)/(v+v_s) )**2 ) / N
v_xi2_rar = \
  sum( ( (v[idx_rarefaction]-v_s[idx_rarefaction])/ \
         (v[idx_rarefaction]+v_s[idx_rarefaction]) )**2 ) / N_rarefaction
v_xi2_con = \
  sum( ( (v[idx_contact]-v_s[idx_contact])/ \
         (v[idx_contact]+v_s[idx_contact]) )**2 ) / N_contact
v_xi2_sho = \
  sum( ( (v[idx_shock]-v_s[idx_shock])/ \
         (v[idx_shock]+v_s[idx_shock]) )**2 ) / N_shock

P_xi2_tot = sum( ( (P-P_s)/(P+P_s) )**2 ) / N
P_xi2_rar = \
  sum( ( (P[idx_rarefaction]-P_s[idx_rarefaction])/ \
         (P[idx_rarefaction]+P_s[idx_rarefaction]) )**2 ) / N_rarefaction
P_xi2_con = \
  sum( ( (P[idx_contact]-P_s[idx_contact])/ \
         (P[idx_contact]+P_s[idx_contact]) )**2 ) / N_contact
P_xi2_sho = \
  sum( ( (P[idx_shock]-P_s[idx_shock])/ \
         (P[idx_shock]+P_s[idx_shock]) )**2 ) / N_shock

print "rho:"
print "xi2 total:", rho_xi2_tot
print "xi2 rarefaction:", rho_xi2_rar
print "xi2 contact:", rho_xi2_con
print "xi2 shock:", rho_xi2_sho
       
print "v:"
print "xi2 total:", v_xi2_tot
print "xi2 rarefaction:", v_xi2_rar
print "xi2 contact:", v_xi2_con
print "xi2 shock:", v_xi2_sho
       
print "P:"
print "xi2 total:", P_xi2_tot
print "xi2 rarefaction:", P_xi2_rar
print "xi2 contact:", P_xi2_con
print "xi2 shock:", P_xi2_sho

times = np.loadtxt("{folder}/timesteps_1.txt".format(folder = folder))

time_total = sum(times[:,6])

print "Total time (ms):", time_total

file = open("{folder}/summary.txt".format(folder = folder), 'w')

file.write("{{\"rho_xi2_tot\": {val},\n".format(val = rho_xi2_tot))
file.write("\"rho_xi2_rar\": {val},\n".format(val = rho_xi2_rar))
file.write("\"rho_xi2_con\": {val},\n".format(val = rho_xi2_con))
file.write("\"rho_xi2_sho\": {val},\n".format(val = rho_xi2_sho))

file.write("\"v_xi2_tot\": {val},\n".format(val = v_xi2_tot))
file.write("\"v_xi2_rar\": {val},\n".format(val = v_xi2_rar))
file.write("\"v_xi2_con\": {val},\n".format(val = v_xi2_con))
file.write("\"v_xi2_sho\": {val},\n".format(val = v_xi2_sho))

file.write("\"P_xi2_tot\": {val},\n".format(val = P_xi2_tot))
file.write("\"P_xi2_rar\": {val},\n".format(val = P_xi2_rar))
file.write("\"P_xi2_con\": {val},\n".format(val = P_xi2_con))
file.write("\"P_xi2_sho\": {val},\n".format(val = P_xi2_sho))

file.write("\"time_total\": {val}}}\n".format(val = time_total))
