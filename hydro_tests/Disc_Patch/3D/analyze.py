################################################################################
# This file is part of SWIFT.
# Copyright (c) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

##
# This script plots the Disc-Patch_*.hdf5 snapshots.
# It takes two (optional) parameters: the counter value of the first and last
# snapshot to plot (default: 0 21).
##

import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import pylab as pl
import glob
import sys

if len(sys.argv) < 2:
  print "No input folder given!"
  exit()

folder = sys.argv[1]

# Parameters
surface_density = 10.
scale_height = 100.
x_disc = 400.
x_trunc = 300.
x_max = 350.
utherm = 20.2678457288
gamma = 5. / 3.

# Get the analytic solution for the density
def get_analytic_density(x):
  return 0.5 * surface_density / scale_height / \
           np.cosh( (x - x_disc) / scale_height )**2

# Get the analytic solution for the (isothermal) pressure
def get_analytic_pressure(x):
  return (gamma - 1.) * utherm * get_analytic_density(x)

# Get the data fields to plot from the snapshot file with the given name:
#  hydro scheme, x-coord, density, pressure, velocity norm
def get_data(name):
  file = h5py.File(name, "r")
  coords = np.array(file["/PartType0/Coordinates"])
  rho = np.array(file["/PartType0/Density"])
  u = np.array(file["/PartType0/InternalEnergy"])
  v = np.array(file["/PartType0/Velocities"])
  scheme = file["/HydroScheme"].attrs["Scheme"]

  P = (gamma - 1.) * rho * u

  vtot = np.sqrt( v[:,0]**2 + v[:,1]**2 + v[:,2]**2 )

  return scheme, coords[:,0], rho, P, vtot

def analyze_file(label):
  scheme, x, rho, P, v = get_data("{0}/Disc-Patch_{1}.hdf5".format(
                                    folder, label))
  numPart = len(x)

  # get the analytic solution
  rho_s = get_analytic_density(x)
  v_s = np.zeros(numPart)
  P_s = get_analytic_pressure(x)

  pl.rcParams["figure.figsize"] = (10, 8)
  pl.rcParams["text.usetex"] = True

  fig, ax = pl.subplots(2, 3, sharex = True)

  xrange = np.arange(x_disc - x_max, x_disc + x_max, 0.002 * x_max)
  ax[0][0].plot(x, rho, "r.")
  ax[0][0].plot(xrange, get_analytic_density(xrange), "k-")
  ax[0][0].plot([x_disc - x_max, x_disc - x_max], [0, 10], "k--", alpha=0.5)
  ax[0][0].plot([x_disc + x_max, x_disc + x_max], [0, 10], "k--", alpha=0.5)
  ax[0][0].plot([x_disc - x_trunc, x_disc - x_trunc], [0, 10], "k--", alpha=0.5)
  ax[0][0].plot([x_disc + x_trunc, x_disc + x_trunc], [0, 10], "k--", alpha=0.5)
  ax[0][0].set_ylim(0., 1.2 * get_analytic_density(x_disc))
  ax[1][0].plot(x, (rho - rho_s) / (rho + rho_s), "k.")

  ax[0][1].plot(x, v, "r.")
  ax[0][1].plot(xrange, np.zeros(len(xrange)), "k-")
  ax[0][1].plot([x_disc - x_max, x_disc - x_max], [0, 10], "k--", alpha=0.5)
  ax[0][1].plot([x_disc + x_max, x_disc + x_max], [0, 10], "k--", alpha=0.5)
  ax[0][1].plot([x_disc - x_trunc, x_disc - x_trunc], [0, 10], "k--", alpha=0.5)
  ax[0][1].plot([x_disc + x_trunc, x_disc + x_trunc], [0, 10], "k--", alpha=0.5)
  ax[0][1].set_ylim(-0.5, 10.)
  ax[1][1].plot(x, (v - v_s) / (v + v_s), "k.")

  ax[0][2].plot(x, P, "r.")
  ax[0][2].plot(xrange, get_analytic_pressure(xrange), "k-")
  ax[0][2].plot([x_disc - x_max, x_disc - x_max], [0, 10], "k--", alpha=0.5)
  ax[0][2].plot([x_disc + x_max, x_disc + x_max], [0, 10], "k--", alpha=0.5)
  ax[0][2].plot([x_disc - x_trunc, x_disc - x_trunc], [0, 10], "k--", alpha=0.5)
  ax[0][2].plot([x_disc + x_trunc, x_disc + x_trunc], [0, 10], "k--", alpha=0.5)
  ax[0][2].set_xlim(0., 2. * x_disc)
  ax[0][2].set_ylim(0., 1.2 * get_analytic_pressure(x_disc))
  ax[1][2].plot(x, (P - P_s) / (P + P_s), "k.")

  ax[0][0].set_title("density")
  ax[0][1].set_title("velocity")
  ax[0][2].set_title("pressure")

  ax[0][0].set_ylabel("value")
  ax[1][0].set_ylabel("absolute difference")
  ax[1][0].set_xlabel("position")
  ax[1][1].set_xlabel("position")
  ax[1][2].set_xlabel("position")

  pl.tight_layout()
  pl.subplots_adjust(top = 0.9)
  pl.suptitle("{0}: {1}, {2} particles".format(label, scheme, numPart))
  pl.savefig("{0}/result_{1}.png".format(folder, label))
  pl.close()

  rho_xi2_tot_array = (rho - rho_s)
  rho_xi2_tot = sum( rho_xi2_tot_array**2 ) / numPart

  v_xi2_tot_array = (v - v_s)
  v_xi2_tot = sum( v_xi2_tot_array**2 ) / numPart

  P_xi2_tot_array = (P - P_s)
  P_xi2_tot = sum( P_xi2_tot_array**2 ) / numPart

  print "## {0} ##".format(label)

  print "rho:"
  print "xi2 total:", rho_xi2_tot
       
  print "v:"
  print "xi2 total:", v_xi2_tot

  print "P:"
  print "xi2 total:", P_xi2_tot

  times = np.loadtxt("{0}/timesteps_16_{1}.txt".format(folder, label))

  time_total = sum(times[:,6])
  print "Total time (ms):", time_total

  file = open("{0}/summary_{1}.txt".format(folder, label), 'w')

  file.write("{{\"rho_xi2_tot\": {val},\n".format(val = rho_xi2_tot))

  file.write("\"v_xi2_tot\": {val},\n".format(val = v_xi2_tot))

  file.write("\"P_xi2_tot\": {val},\n".format(val = P_xi2_tot))

  file.write("\"time_total\": {val}}}\n".format(val = time_total))

  file.close()

analyze_file("ic")
analyze_file("run")
