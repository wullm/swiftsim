/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef LOGGER_LOGGER_PARTICLE_H
#define LOGGER_LOGGER_PARTICLE_H

#include "logger_header.h"
#include "logger_time.h"
#include "logger_tools.h"

#include <stdio.h>
#include <stdlib.h>

#if defined(HYDRO_DIMENSION_1D)
#define DIM 1
#elif defined(HYDRO_DIMENSION_2D)
#define DIM 2
#elif defined(HYDRO_DIMENSION_3D)
#define DIM 3
#endif

struct logger_reader;

/**
 * @brief Store the data from a record.
 *
 * This structure contains all the required fields
 * present in a file.
 *
 * As we need only a few particles, no need to keep it small.
 *
 * The particle is initialized with #logger_particle_init
 * and can be updated with a record through #logger_particle_read.
 *
 * In #logger_particle_read, we use #logger_particle_read_field on
 * each field and #logger_particle_interpolate if a linear
 * interpolation is required.
 */
struct logger_particle {
  /* position. */
  double pos[DIM];

  /* velocity. */
  float vel[DIM];

  /* acceleration. */
  float acc[DIM];

  /* entropy. */
  float entropy;

  /* smoothing length. */
  float h;

  /* density. */
  float density;

  /* mass. */
  float mass;

  /* unique id. */
  size_t id;

  /* time of the record. */
  double time;
};

/**
 * @brief Defines the type of interpolation
 */
enum logger_reader_type {
  logger_reader_const, /* Constant interpolation. */
  logger_reader_lin,   /* Linear interpolation. */
};

void logger_particle_print(const struct logger_particle *p);

size_t logger_particle_read(struct logger_particle *part,
                            const struct logger_reader *reader, size_t offset,
                            const double time,
                            const enum logger_reader_type reader_type);

void logger_particle_init(struct logger_particle *part);

void *logger_particle_read_field(struct logger_particle *part, void *map,
                                 const char *field, const size_t size);

void logger_particle_interpolate(struct logger_particle *part_curr,
                                 const struct logger_particle *part_next,
                                 const double time);

#endif  // LOGGER_LOGGER_PARTICLE_H
