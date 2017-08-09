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

# contains the available setups for this test as a dictionary with desired
# numbers of cells and glass file combinations to generate them
# this file can either be imported in another script to use the dictionary
# directly, or it can be called from the command line, in which case it will
# just print a list of available particle numbers (that can be used in bash
# scripts) and exit.

# supported setups:
# the key is the total number of particles in the setup, the underlying
# dictionary tells us which glass files to use and how many copies of them to
# make
setups = {
  9216: {"small_glass": 8, "large_glass": 16, "num_copy": 1},
  73728: {"small_glass": 16, "large_glass": 32, "num_copy": 1},
  589824: {"small_glass": 32, "large_glass": 64, "num_copy": 1},
  4718592: {"small_glass": 64, "large_glass": 128, "num_copy": 1},
#  37748736: {"small_glass": 64, "large_glass": 128, "num_copy": 2}
}

# main method. If the script is called from the command line, we just print the
# list of available setup sizes
if __name__ == "__main__":
  for key in sorted(setups):
    print key,
