/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2017 Stefan Arridge (stefan.arridge@durham.ac.uk)
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
#ifndef SWIFT_STAR_FORMATION
#define SWIFT_STAR_FORMATION

/**
 * @file src/star_formation.h
 * @brief Branches between the different cooling functions.
 */

/* Config parameters. */
#include "../config.h"

/* Import the right cooling definition */
#if defined(STAR_FORMATION_NONE)
#include "./star_formation/none/star_formation.h"
#elif defined(STAR_FORMATION_DENS_THRESH)
#include "./star_formation/dens_thresh/star_formation.h"
#else
#error "Invalid choice of cooling function."
#endif

/* Common functions */
void star_formation_init(const struct swift_params* parameter_file,
                  const struct UnitSystem* us,
                  const struct phys_const* phys_const,
                  struct star_formation_data* star_formation);

void star_formation_print(const struct star_formation_data* star_formation);

#endif /* SWIFT_STAR_FORMATION */
