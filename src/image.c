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

#include "engine.h"
#include "hdf5.h"
#include "image.h"
#include "minmax.h"
#include "part.h"
#include "space.h"

/* Taken from Dehnen & Aly 2012 */
#define imaging_kernel_gamma 1.897367
#define imaging_kernel_constant 2.228171
#define NUM_PIXELS 512

/**
 * @brief Kernel for imaging. Returns the kernel in correct units (i.e.
 *        the local number density is the sum of these for all particles).
 *
 *        This is the Wendland-C2 kernel in 2D.
 *
 * @props r, radius
 * @props H, kernel cut-off radius
 */
__attribute__((always_inline)) INLINE static float imaging_kernel(float r,
                                                                  float H) {
  const float inverse_H = 1.f / H;
  const float ratio = r * inverse_H;

  float kernel = 0.0;

  if (ratio < 1.0) {
    const float one_minus_ratio = 1.f - ratio;
    const float one_minus_ratio_2 = one_minus_ratio * one_minus_ratio;
    const float one_minus_ratio_4 = one_minus_ratio_2 * one_minus_ratio_2;

    kernel = max(one_minus_ratio_4 * (1.f + 4.f * ratio), 0.f);

    kernel *= imaging_kernel_constant * inverse_H * inverse_H;
  }

  return kernel;
}

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
  const int number_of_pixels[2] = {NUM_PIXELS, NUM_PIXELS};
  const size_t total_number_of_pixels =
      (size_t)number_of_pixels[0] * number_of_pixels[1];

  /* Allocate main image array */
  float* image = (float*)malloc(total_number_of_pixels * sizeof(float));
  bzero(image, total_number_of_pixels * sizeof(float));

  /* Parameters for the imaging */
  const float pixel_width_x = box_size[0] / number_of_pixels[0];
  const float pixel_width_y = box_size[1] / number_of_pixels[1];
  const float pixel_area = pixel_width_x * pixel_width_y;

  /* Conversion factors to number of pixels. Required as float / int comparisons
   * can cause some trickery */
  const float x_conversion_fac = number_of_pixels[0] / box_size[0];
  const float y_conversion_fac = number_of_pixels[1] / box_size[1];

  const float drop_to_single_cell = pixel_width_x * 0.5;

  /* Perform a scatter over all particles */
  for (size_t particle = 0; particle < num_parts; particle++) {
    /* Grab the particle and extract properties! */
    const struct part* p = &s->parts[particle];

    /* Need to make sure the box is wrapped! Especially clear in e.g. a
     * Kelvin Helmholtz test. */
    const double x_wrapped = box_wrap(p->x[0], 0.0, box_size[0]);
    const double y_wrapped = box_wrap(p->x[1], 0.0, box_size[1]);

    /* What we really need is an integer position for comparison to the pixel
     * grid */
    const int x = x_wrapped * x_conversion_fac;
    const int y = y_wrapped * y_conversion_fac;
    const float kernel_width = kernel_gamma * p->h / box_size[0];

    if (kernel_width <= drop_to_single_cell) {
      /* Simple case; we average our mass over the size of the pixel */
      const int pixel = x + number_of_pixels[0] * y;
      image[pixel] += p->mass / pixel_area;
    } else {
      /* More complex; we need to SPH-smooth the data. This follows the python
       * version in swiftsimio very closely */

      /* May run into problems with non-square grids here... */
      const int cells_spanned = (int)(1.f + kernel_width * number_of_pixels[0]);

      /* Ensure our loop bounds stay within the x, y grid */
      const int starting_x = max(0, x - cells_spanned);
      const int ending_x = max(x + cells_spanned, number_of_pixels[0]);
      const int starting_y = max(0, y - cells_spanned);
      const int ending_y = max(y + cells_spanned, number_of_pixels[1]);

      /* Now loop over all cells (the square) that our #part covers in the
       * final image */
      for (int cell_x = starting_x; cell_x < ending_x; cell_x++) {
        /* Precompute properties that are constant in x */
        const float distance_x = (cell_x + 0.5) * pixel_width_x - x_wrapped;
        const float distance_x_2 = distance_x * distance_x;

        for (int cell_y = starting_y; cell_y < ending_y; cell_y++) {
          const float distance_y = (cell_y + 0.5) * pixel_width_y - y_wrapped;
          const float distance_y_2 = distance_y * distance_y;

          const float r = sqrtf(distance_x_2 + distance_y_2);
          const float kernel_eval = imaging_kernel(r, kernel_width);

          const int pixel = cell_x + number_of_pixels[0] * cell_y;
          image[pixel] += kernel_eval * p->mass;
        }
      }
    }
  }

  // BASIC I/O, Don't keep this lol
  char fileName[FILENAME_BUFFER_SIZE];
  snprintf(fileName, FILENAME_BUFFER_SIZE, "%s_%04i.dump", "image",
           e->snapshot_output_count);
  FILE* f = fopen(fileName, "wb");
  fwrite(image, sizeof(float), number_of_pixels[0] * number_of_pixels[1], f);
  fclose(f);

  // HDF5 stuff
  // /* For now, we'll just dump this as a single array to HDF5. */
  // h_file = H5Fcreate("IMAGE_OUT_TEST.hdf5", H5F_ACC_TRUNC, H5P_DEFAULT,
  // H5P_DEFAULT); if (h_file < 0) {
  //   error("Error while opening file image file");
  // }

  // h_grp = H5Gcreate(h_file, "/Images", H5P_DEFAULT, H5P_DEFAULT,
  // H5P_DEFAULT); if (h_grp < 0) error("Error while creating group\n");

  // /* Create data space */
  // const hid_t h_space = H5Screate(H5S_SIMPLE);
  // if (h_space < 0)
  //   error("Error while creating data space for field '%s'.", props.name);

  // H5Gclose(h_grp);
  // H5Fclose(h_file);

  free(image);
}