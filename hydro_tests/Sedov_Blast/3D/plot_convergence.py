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

# Plots the convergence rate for the 3D Sedov blast test using different numbers
# of cells

import numpy as np
import matplotlib
matplotlib.use("Agg")
import pylab as pl
import setups

pl.rcParams["figure.figsize"] = (4, 4)
pl.rcParams["text.usetex"] = True

sims = {"gizmo": ["GIZMO", "#d7191c"],
        "gadget2": ["Gadget2", "#fdae61"],
        "hopkins": ["Pressure-entropy SPH", "#2c7bb6"]
       }

ncell = sorted(setups.setups.keys())

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
  ax.semilogy(ncell, [data[sim][n]["rho_xi2_tot"] for n in ncell],
              "o-", color = sims[sim][1], label = sims[sim][0])

ax.set_title("3D Sedov blast convergence")
ax.set_xlabel("Number of particles")
ax.set_ylabel(r"$\langle{}\chi{}^2\rangle{}$")
ax.legend(loc = "best")
pl.tight_layout()
pl.savefig("SedovBlast_3D_convergence.png")
