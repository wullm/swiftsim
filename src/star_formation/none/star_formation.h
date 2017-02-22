
/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2017 Stefan Arridge (stefan.arridge@durham.ac.uk)

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
#ifndef SWIFT_STAR_FORMATION_NONE_H
#define SWIFT_STAR_FORMATION_NONE_H

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <float.h>

/* Local includes. */
#include "error.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "units.h"
#include "cell.h"
#include "hydro.h"
#include "stragglers.h"

/**
 * @brief Star Formation Properties
 */
struct star_formation_data {};

/**
 * @brief Does nothing
 *        
 * @param star_formation The proerties of the star formation model.
 * @param phys_const The physical constants in internal units.
 * @param p Pointer to the particle data.
 * @param c Pointer to the cell to which the particle belongs
 */
__attribute__((always_inline)) INLINE static void do_star_formation(
    const struct star_formation_data* restrict star_formation,
    const struct phys_const* restrict phys_const, struct stragglers* restrict stragglers,
    struct part* restrict p, struct cell* restrict c, float dt) {}



/**
 * @brief Initialises the star formation properties in the internal system
 * of units.
 *
 * Nothing to do here.
 *
 * @param parameter_file The parsed parameter file
 * @param phys_const Physical constants in internal units
 * @param us The current internal system of units
 * @param star_formation The star formation properties to initialize
 */
static INLINE void star_formation_init_backend(
    const struct swift_params* parameter_file,
    const struct UnitSystem* us,
    const struct phys_const* phys_const,
    struct star_formation_data* star_formation) {}

/**
 * @brief Prints the properties of the star formation model to stdout.
 *
 * @param star_formation The star formation properties.
 */
static INLINE void star_formation_print_backend(
    const struct star_formation* star_formation) {

  message("Star formation model is 'No star formation'.");
}

#endif /* SWIFT_STAR_FORMATION_NONE_H */
