/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2015 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_GRAVITY_H
#define SWIFT_GRAVITY_H

/* Config parameters. */
#include "../config.h"

/* Local headers. */
#include "const.h"
#include "inline.h"
#include "part.h"

/* So far only one model here */
/* Straight-forward import */
#include "./gravity/Default/gravity.h"
#include "./gravity/Default/gravity_iact.h"

struct engine;
struct space;

void gravity_exact_force_ewald_init(double boxSize);
void gravity_exact_force_ewald_free(void);
void gravity_exact_force_ewald_evaluate(double rx, double ry, double rz,
                                        double corr_f[3], double *corr_p);
void gravity_exact_force_compute(struct space *s, const struct engine *e);
void gravity_exact_force_check(struct space *s, const struct engine *e,
                               float rel_tol);

#endif
