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

#include <stdio.h>

#include "image.h"
#include "engine.h"
#include "space.h"
#include "part.h"
#include "hdf5.h"

/**
 * @brief dumps a projected image.
 * 
 * @props e, pointer to the engine.
 **/
void image_dump_image(struct engine* e) {
  const struct space* s = e->s;
  const size_t num_parts = s->nr_parts;
  const double box_size[2] = {s->dim[0], s->dim[1]};

  /* Should change these to user-defined parameters... */
  const int number_of_pixels[2] = {128, 128};
  const size_t total_number_of_pixels = (size_t) number_of_pixels[0] * number_of_pixels[1];

  const float pixel_area = (float)(box_size[0] / ((double)number_of_pixels[0]) *
                                   box_size[1] / ((double)number_of_pixels[1]));

  float* image = (float*)malloc( total_number_of_pixels *
                              sizeof(float));

  /* Need to first zero our memory */
  for (size_t pixel = 0; pixel < total_number_of_pixels; pixel++) {
    image[pixel] = 0.f;
  }

  for (size_t particle = 0; particle < num_parts; particle++) {
    /* Grab the particle and extract properties! */
    const struct part* p = &s->parts[particle];
    const int x = number_of_pixels[0] * (p->x[0] / box_size[0]);
    const int y = number_of_pixels[1] * (p->x[1] / box_size[1]);
    const float density_contribution = p->mass / pixel_area;

    /* Now need to figure out which pixel this friendly neighbourhood
     * #part icle belongs in */

    const int pixel = x + number_of_pixels[0] * y;

    /* Add on the density contribution to the individual pixel */
    image[pixel] += density_contribution;
  }

  // BASIC I/O, Don't keep this lol
  char fileName[FILENAME_BUFFER_SIZE];
  snprintf(fileName, FILENAME_BUFFER_SIZE, "%s_%04i.dump", "image",
           e->snapshot_output_count);
  FILE *f = fopen(fileName, "wb");
  fwrite(image, sizeof(float), number_of_pixels[0] * number_of_pixels[1], f);
  fclose(f);

  // HDF5 stuff
  // /* For now, we'll just dump this as a single array to HDF5. */
  // h_file = H5Fcreate("IMAGE_OUT_TEST.hdf5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  // if (h_file < 0) {
  //   error("Error while opening file image file");
  // }

  // h_grp = H5Gcreate(h_file, "/Images", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  // if (h_grp < 0) error("Error while creating group\n");

  // /* Create data space */
  // const hid_t h_space = H5Screate(H5S_SIMPLE);
  // if (h_space < 0)
  //   error("Error while creating data space for field '%s'.", props.name);

  // H5Gclose(h_grp);
  // H5Fclose(h_file);

  free(image);
}