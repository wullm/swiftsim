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

# Plots the timings for the 3D disc patch test using different numbers of cells

import numpy as np
import matplotlib
matplotlib.use("Agg")
import pylab as pl
import setups

pl.rcParams["figure.figsize"] = (10, 5)
pl.rcParams["text.usetex"] = True

sims = {"gizmo": ["GIZMO", "#d7191c"],
        "gadget2": ["Gadget2", "#fdae61"],
        "hopkins": ["Pressure-entropy SPH", "#2c7bb6"]
       }

ncell = sorted(setups.setups.keys())

def get_total_time(sim, label):
  return sum(np.loadtxt("{sim}/timesteps_16_{label}.txt".format
                          sim = sim, label = label))[:,6])

times = {}
for sim in sims:
  times[sim] = {"ic": np.zeros(len(ncell), "run": np.zeros(len(ncell))}
  for i in range(len(ncell)):
    times[sim]["ic"][i] = get_total_time("{0}_{1}".format(sim, ncell[i]), "ic")
    times[sim]["run"][i] = \
      get_total_time("{0}_{1}".format(sim, ncell[i]), "run")

fig, ax = pl.subplots(1, 2, sharex = True)

for sim in sims:
  ax[0].semilogy(ncell, times[sim]["ic"], "o-", color = sims[sim][1],
                 label = sims[sim][0])
  ax[1].semilogy(ncell, times[sim]["run"], "o-", color = sims[sim][1],
                 label = sims[sim][0])

ax[0].set_xlabel("Number of particles")
ax[1].set_xlabel("Number of particles")
ax[0].set_ylabel("Runtime (ms)")
ax[0].legend(loc = "best")
ax[0].set_title("initial condition")
ax[1].set_title("run")
pl.tight_layout()
pl.subplots_adjust(top = 0.9)
pl.suptitle("3D Sedov blast timings")
pl.savefig("SedovBlast_3D_timings.png")
pl.close()
