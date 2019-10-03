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
#include "lock.h"

/**
 * @brief Image parameters read from the parameter file.
 */
struct image_props {

  /*! Size of the resultant image (in pixels) */
  int image_size[2];

  /*! HDF5 image filename for images */
  char image_filename[PARSER_MAX_LINE_SIZE];

  /*! Base name inside snapshots for images */
  char image_base_name[PARSER_MAX_LINE_SIZE];
};

/*! Parameters for rendering the image, generally made up of information from
 * existing structs */
struct render_properties {

  /*! Boxsize */
  double box_size[2];
};

/**
 * @brief Image data to pass to threadpool.
 */
struct image_data {

  /*! Pointer to final image */
  float* image;

  /*! Properties of the output image */
  struct image_props image_properties;

  /*! Properties used to render the image */
  struct render_properties render_properties;

  /*! Lock for when merging individual images into larger whole */
  swift_lock_type lock;
};

static float imaging_kernel(float r, float H);

void image_dump_image(struct engine* e);

static void image_add(float* image_from, float* image_to, size_t size);

void image_init(struct swift_params* params, const struct unit_system* us,
                const struct phys_const* phys_const,
                struct image_props* image_properties);

void image_props_struct_dump(const struct image_props* p, FILE* stream);
void image_props_struct_restore(const struct image_props* p, FILE* stream);

#endif /* SWIFT_IMAGE_H */