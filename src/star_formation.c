/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *           (c) 2017 Stefan Arridge (stefan.arridge@durham.ac.uk)
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

/* Config parameters. */
#include "../config.h"

/* This object's header. */
#include "star_formation.h"

/**
 * @brief Initialises the star formation properties.
 *
 * Calls star_formation_init_backend for the chosen star formation model.
 *
 * @param parameter_file The parsed parameter file.
 * @param us The current internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param star_formation The star formation properties to initialize
 */
void star_formation_init(const struct swift_params* parameter_file,
                  const struct UnitSystem* us,
                  const struct phys_const* phys_const,
                  struct star_formation_data* star_formation) {

  star_formation_init_backend(parameter_file, us, phys_const, star_formation);
}

/**
 * @brief Prints the properties of the star_formation model to stdout.
 *
 * Calls star_formation_print_backend for the chosen star_formation function.
 *
 * @param star_formation The properties of the star_formation model.
 */
void star_formation_print(const struct star_formation_data* star_formation) {

  star_formation_print_backend(star_formation);
}
