/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Loic Hausammann (loic.hausammann@epfl.ch)
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

#include "logger_hydro.h"

int hydro_logger_local_to_global[hydro_logger_field_count];
/* Define the size of all the fields. */
#define member_size(type, member) sizeof(((type *)0)->member)

const int hydro_logger_field_size[hydro_logger_field_count] = {
    member_size(struct part, x),        // coordinates
    member_size(struct part, v),        // velocities
    member_size(struct part, a_hydro),  // accelerations
    member_size(struct part, mass),     // massses
    member_size(struct part, h),        // Smoothing Length
    member_size(struct part, entropy),  // Entropy
    member_size(struct part, id),       // IDs
    member_size(struct part, rho),      // density
};
