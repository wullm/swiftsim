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

# Plots the timings for the 1D Noh test using different numbers of particles

import numpy as np
import matplotlib
matplotlib.use("Agg")
import pylab as pl

pl.rcParams["figure.figsize"] = (4, 4)
pl.rcParams["text.usetex"] = True

sims = {"gizmo": ["GIZMO", "#d7191c"],
        "gadget2": ["Gadget2", "#fdae61"],
        "hopkins": ["Pressure-entropy SPH", "#2c7bb6"]
       }

ncell = np.array([100, 200, 400, 800, 1600, 3200])

def get_total_time(sim):
  return sum(np.loadtxt("{sim}/timesteps_1.txt".format(sim = sim))[:,6])

times = {}
for sim in sims:
  times[sim] = np.zeros(len(ncell))
  for i in range(len(ncell)):
    times[sim][i] = get_total_time("{0}_{1}".format(sim, ncell[i]))

fig, ax = pl.subplots(1, 1)

for sim in sims:
  ax.semilogy(ncell, times[sim], "o-", color = sims[sim][1],
              label = sims[sim][0])

ax.set_xlabel("Number of particles")
ax.set_ylabel("Runtime (ms)")
ax.legend(loc = "best")
ax.set_title("1D Noh timings")
pl.tight_layout()
pl.savefig("Noh_1D_timings.png")
pl.close()
