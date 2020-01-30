/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Willem Elbers (whe@willemelbers.com)
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
#ifndef SWIFT_NEUTRINO_PHASE_SPACE_H
#define SWIFT_NEUTRINO_PHASE_SPACE_H

/* The general cosmology header */
#include "../cosmology.h"

/* Engine header */
#include "../engine.h"

double fermi_dirac_density(const struct engine* engine, double* x, float* v);
double sample_density(const struct engine* engine, double* x, float* v);
double fermi_dirac_momentum(const struct engine* engine, float* v);
double fermi_dirac_energy(const struct engine* engine, float* v);

#endif /* SWIFT_NEUTRINO_PHASE_SPACE_H */
