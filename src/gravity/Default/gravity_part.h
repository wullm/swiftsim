/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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
#ifndef SWIFT_DEFAULT_GRAVITY_PART_H
#define SWIFT_DEFAULT_GRAVITY_PART_H

/* Gravity particle. */
struct gpart {

  /*! Particle ID. If negative, it is the negative offset of the #part with
     which this gpart is linked. */
  long long id_or_neg_offset;

  /*! Particle position. */
  double x[3];

  /*! Offset between current position and position at last tree rebuild. */
  float x_diff[3];

  /*! Particle velocity. */
  float v_full[3];

  /*! Particle acceleration. */
  float a_grav[3];

  /*! Particle mass. */
  float mass;

  /*! Gravitational potential */
  float potential;

  /*! Time-step length */
  timebin_t time_bin;

  /*! Type of the #gpart (DM, gas, star, ...) */
  enum part_type type;

  /*! Particle offset into FOF group id array. */
  size_t offset;

#ifdef SWIFT_DEBUG_CHECKS

  /* Numer of gparts this gpart interacted with */
  long long num_interacted;

  /* Time of the last drift */
  integertime_t ti_drift;

  /* Time of the last kick */
  integertime_t ti_kick;

#endif

#ifdef SWIFT_GRAVITY_FORCE_CHECKS

  /* Brute-force particle acceleration. */
  double a_grav_exact[3];

  /* Brute-force particle potential. */
  double potential_exact;
#endif

} SWIFT_STRUCT_ALIGN;

#endif /* SWIFT_DEFAULT_GRAVITY_PART_H */
