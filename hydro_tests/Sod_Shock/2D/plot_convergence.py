################################################################################
# This file is part of SWIFT.
# Copyright (c) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

# Plots the convergence rate for the 1D Sod shock test using different numbers
# of cells
# Note that the indicated number of cells is the number of cells in the high
# density region; the total number of cells is larger as there also is a low
# density region.
# We plot the total convergence rate and the convergence rates for the three
# characteristic waves present in the solution.

import numpy as np
import matplotlib
matplotlib.use("Agg")
import pylab as pl

sims = {"gizmo": ["GIZMO", 'r'],
        "gadget2": ["Gadget2", 'g'],
        "hopkins": ["Pressure-entropy SPH", 'b']
       }

ncell = np.array([1, 2, 3, 4, 5])
ncell_real = (128**2 * 2 + 48**2 * 2) * ncell**2

def get_data(sim):
  file = open("{sim}/summary.txt".format(sim = sim), 'r')
  return eval(file.read())

data = {}
for sim in sims:
  data[sim] = {}
  for n in ncell:
    data[sim][n] = get_data("{sim}_{n}".format(sim = sim, n = n))

fig, ax = pl.subplots(2, 2, sharex = True)

for sim in sims:
  ax[0][0].semilogy(ncell_real, [data[sim][n]["rho_xi2_tot"] for n in ncell],
                    "o-", color = sims[sim][1], label = sims[sim][0])
  ax[0][1].semilogy(ncell_real, [data[sim][n]["rho_xi2_rar"] for n in ncell],
                    "o-", color = sims[sim][1], label = sims[sim][0])
  ax[1][0].semilogy(ncell_real, [data[sim][n]["rho_xi2_con"] for n in ncell],
                    "o-", color = sims[sim][1], label = sims[sim][0])
  ax[1][1].semilogy(ncell_real, [data[sim][n]["rho_xi2_sho"] for n in ncell],
                    "o-", color = sims[sim][1], label = sims[sim][0])

ax[0][0].set_title("Total")
ax[0][1].set_title("Rarefaction")
ax[1][0].set_title("Contact")
ax[1][1].set_title("Shock")
dncell = 0.1 * (ncell_real[-1] - ncell_real[0])
ax[0][0].set_xlim(ncell_real[0] - dncell, ncell_real[-1] + dncell)
ax[1][0].set_xlabel("Number of particles")
ax[1][1].set_xlabel("Number of particles")
ax[0][0].set_ylabel(r"$\langle{}\chi{}^2\rangle{}$")
ax[1][0].set_ylabel(r"$\langle{}\chi{}^2\rangle{}$")
ax[0][0].legend(loc = "best")
pl.suptitle("2D Sod shock convergence")
pl.savefig("SodShock_2D_convergence.png")
