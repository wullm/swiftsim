/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2019 Josh Borrow (joshua.borrow@durham.ac.uk)
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

#ifndef SWIFT_IMAGE_H
#define SWIFT_IMAGE_H

/* Config parameters. */
#include "../config.h"

/* Local headers. */
#include "const.h"
#include "engine.h"

/**
 * @brief Image data to pass to threadpool.
 */
struct image_data {
  /*! Pointer to final image */
  float* image;

  /*! Size of resultant image */
  int image_size[2];

  /*! Boxsize */
  double box_size[2];

  /*! Drop to single cell factor */
  float drop_to_single_cell_factor;

  /*! Lock for when merging individual images into larger whole */
  swift_lock_type lock;
};

static float imaging_kernel(float r, float H);

void image_dump_image(struct engine* e);
static void image_add(float* image_from, float* image_to, size_t size);

#endif /* SWIFT_IMAGE_H */