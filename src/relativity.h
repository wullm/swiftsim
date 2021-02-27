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
#ifndef SWIFT_RELATIVITY_H
#define SWIFT_RELATIVITY_H

/* Config parameters. */
#include "../config.h"

#include "engine.h"

/* Calculation of sqrt(x^2 + y^2 + z^2), without undue overflow or underflow. */
__attribute__((always_inline, const)) INLINE static double hypot3(double x,
                                                                  double y,
                                                                  double z) {

  return hypot(x, hypot(y, z));
}


/**
 * @brief Calculate the relativistic correction to the 'drift' timestep
 *
 * @param e The engine
 * @param V The scale factor
 */
__attribute__((always_inline)) INLINE static double relat_corr_drift(
    const struct engine *e, const float *V) {

  /* Perform a relativistic correction, see eq. (5.13) in 1604.06065 */
  double a = e->cosmology->a;
  double c = e->physical_constants->const_speed_light_c;
  float v = hypot3(V[0], V[1], V[2]) / (a * c);

  /* When this is enabled, the internal velocity variable is
   * the spatial component of the 4-velocity U^i, multiplied by a^2,
   * which reduces to the usual comoving velocity times a^2 in the
   * limit of c --> infinity.
   */

  return 1. / hypot(v, 1.);
}

/**
 * @brief Calculate the relativistic correction to the 'kick' timestep
 *
 * @param e The engine
 * @param V The scale factor
 */
__attribute__((always_inline)) INLINE static double relat_corr_kick(
    const struct engine *e, const float *V, const float V_i) {

  /* Perform a relativistic correction, see eq. (5.14) in 1604.06065 */
  double a = e->cosmology->a;
  double c = e->physical_constants->const_speed_light_c;
  float v = hypot3(V[0], V[1], V[2]) / (a * c);
  (void) v;

  float v_i = V_i / (a * c);

  /* When this is enabled, the internal velocity variable is
   * the spatial component of the 4-velocity U^i, multiplied by a^2,
   * which reduces to the usual comoving velocity times a^2 in the
   * limit of c --> infinity.
   */

  return (2 * v_i * v_i + 1.) / hypot(v_i, 1.);
}

/**
 * @brief Calculate the relativistic correction to the 3-acceleration, i.e.
 * the second time derivative of the particle coordinates. This is only used
 * in get_gpart_timestep, to correct the Gadget timestep criterion ~1/sqrt(a)
 *
 * @param e The engine
 * @param V The scale factor
 */
__attribute__((always_inline)) INLINE static double relat_corr_3accel(
    const struct engine *e, const float *V, const double V_i) {

  /* Perform a relativistic correction */
  double a = e->cosmology->a;
  double c = e->physical_constants->const_speed_light_c;
  float vv = (V[0] * V[0] + V[1] * V[1] + V[2] * V[2]) / (a * a * c * c);
  (void) vv;

  float v_i = V_i / (a * c);
  float vv_i = v_i * v_i;

  /* When this is enabled, the internal velocity variable is
   * the spatial component of the 4-velocity U^i, multiplied by a^2,
   * which reduces to the usual comoving velocity times a^2 in the
   * limit of c --> infinity.
   */

  /* Derivative of v / sqrt(v^2 + 1), using v' = (2v^2 + 1)a / sqrt(v^2 + 1) */
  return (2 * vv_i + 1.) / ((vv_i + 1.0) * (vv_i + 1.0));
}

#endif /* SWIFT_RELATIVITY_H */
