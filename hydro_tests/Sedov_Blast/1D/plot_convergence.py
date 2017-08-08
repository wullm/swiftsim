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

ncell = np.array([100, 200, 400])
ncell_real = ncell + np.floor(0.125*ncell)

def get_data(sim):
  file = open("{sim}/summary.txt".format(sim = sim), 'r')
  return eval(file.read())

data = {}
for sim in sims:
  data[sim] = {}
  for n in ncell:
    data[sim][n] = get_data("{sim}_{n}".format(sim = sim, n = n))

fig, ax = pl.subplots(1, 1, sharex = True)

for sim in sims:
  ax.semilogy(ncell_real, [data[sim][n]["rho_xi2_tot"] for n in ncell],
                    "o-", color = sims[sim][1], label = sims[sim][0])

ax.set_title("Total")
dncell = 0.1 * (ncell_real[-1] - ncell_real[0])
ax.set_xlim(ncell_real[0] - dncell, ncell_real[-1] + dncell)
ax.set_xlabel("Number of particles")
ax.set_ylabel(r"$\langle{}\chi{}^2\rangle{}$")
ax.legend(loc = "best")
ax.set_xticks([0, 1000, 2000, 3000])
pl.suptitle("1D Sedov blast convergence")
pl.savefig("SedovBlast_1D_convergence.png")
