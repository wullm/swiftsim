/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Willem Elbers (whe@willemelbers.com)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#ifndef SWIFT_DEFAULT_FERMI_DIRAC_H
#define SWIFT_DEFAULT_FERMI_DIRAC_H

/* Standard headers */
#include <math.h>
#include <stdint.h>

/* Faster exponential */
#include "exp.h"

/**
 * @brief Calculate the neutrino density at the particle's location in phase
 * space, according to the 0th order background model: f_0(x,p,t).
 *
 * @param z Argument of the FD function: z = p/kbT
 */
INLINE static float fermi_dirac_density(float z) {
  return 1.0 / (optimized_expf(z) + 1.0);
}

double neutrino_seed_to_fermi_dirac(uint64_t seed);
void neutrino_seed_to_direction(uint64_t seed, double n[3]);

#endif /* SWIFT_DEFAULT_FERMI_DIRAC_H */
