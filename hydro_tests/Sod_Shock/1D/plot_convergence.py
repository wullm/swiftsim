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
import pylab as pl

ncell = np.array([100, 200, 400, 800, 1600, 3200])

def get_data(sim):
  file = open("{sim}/summary.txt".format(sim = sim), 'r')
  return eval(file.read())

gizmo = {}
gadget2 = {}
minimal = {}
hopkins = {}
for n in ncell:
  gizmo[n] = get_data("gizmo_{n}".format(n = n))
  gadget2[n] = get_data("gadget2_{n}".format(n = n))
  minimal[n] = get_data("minimal_{n}".format(n = n))
  hopkins[n] = get_data("hopkins_{n}".format(n = n))

fig, ax = pl.subplots(2, 2, sharex = True)

ax[0][0].semilogy(ncell, [gizmo[n]["rho_xi2_tot"] for n in ncell], "ro-")
ax[0][1].semilogy(ncell, [gizmo[n]["rho_xi2_rar"] for n in ncell], "ro-")
ax[1][0].semilogy(ncell, [gizmo[n]["rho_xi2_con"] for n in ncell], "ro-")
ax[1][1].semilogy(ncell, [gizmo[n]["rho_xi2_sho"] for n in ncell], "ro-")
ax[0][0].semilogy(ncell, [gadget2[n]["rho_xi2_tot"] for n in ncell], "go-")
ax[0][1].semilogy(ncell, [gadget2[n]["rho_xi2_rar"] for n in ncell], "go-")
ax[1][0].semilogy(ncell, [gadget2[n]["rho_xi2_con"] for n in ncell], "go-")
ax[1][1].semilogy(ncell, [gadget2[n]["rho_xi2_sho"] for n in ncell], "go-")
ax[0][0].semilogy(ncell, [minimal[n]["rho_xi2_tot"] for n in ncell], "yo-")
ax[0][1].semilogy(ncell, [minimal[n]["rho_xi2_rar"] for n in ncell], "yo-")
ax[1][0].semilogy(ncell, [minimal[n]["rho_xi2_con"] for n in ncell], "yo-")
ax[1][1].semilogy(ncell, [minimal[n]["rho_xi2_sho"] for n in ncell], "yo-")
ax[0][0].semilogy(ncell, [hopkins[n]["rho_xi2_tot"] for n in ncell], "bo-")
ax[0][1].semilogy(ncell, [hopkins[n]["rho_xi2_rar"] for n in ncell], "bo-")
ax[1][0].semilogy(ncell, [hopkins[n]["rho_xi2_con"] for n in ncell], "bo-")
ax[1][1].semilogy(ncell, [hopkins[n]["rho_xi2_sho"] for n in ncell], "bo-")

ax[0][0].set_title("Total")
ax[0][1].set_title("Rarefaction")
ax[1][0].set_title("Contact")
ax[1][1].set_title("Shock")
ax[0][0].set_xlim(0., 3300.)
pl.show()
