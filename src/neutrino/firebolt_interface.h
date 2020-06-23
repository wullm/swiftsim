/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Willem Elbers (whe@willemelbers.com).
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

/**
 * firebolt_interface.h  -  call Firebolt for linear theory distribution function
 */

#ifndef SWIFT_FIREBOLT_INTERFACE_H
#define SWIFT_FIREBOLT_INTERFACE_H

#include "../engine.h"

/* The renderer object renders transfer functions onto the grid */
#include "renderer.h"

#include <firebolt_min.h>

int firebolt_init(struct swift_params *params, struct renderer *rend, const struct engine *e);
int firebolt_update(const struct renderer *rend, const struct engine *e);
int firebolt_free(void);

#endif /* SWIFT_FIREBOLT_INTERFACE_H */
