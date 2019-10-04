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

#ifdef HAVE_HDF5
#include <hdf5.h>
#endif

#ifdef WITH_MPI
#include <mpi.h>
#endif

#include "engine.h"
#include "image.h"
#include "minmax.h"
#include "part.h"
#include "space.h"

/* Taken from Dehnen & Aly 2012 */
#define imaging_kernel_gamma 1.897367
#define imaging_kernel_constant 2.228171

/**
 * @brief Kernel for imaging. Returns the kernel in correct units (i.e.
 *        the local number density is the sum of these for all particles).
 *
 *        This is the Wendland-C2 kernel in 2D.
 *
 * @props r, radius
 * @props H, kernel cut-off radius
 */
__attribute__((always_inline)) INLINE float imaging_kernel(float r, float H) {
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
 * @brief Adds two images together
 *
 * @param image_from, the image that remains fixed and values are taken from
 * @param image_to, the image that is modified and has things added to it
 * @param size, the size of the images.
 */
__attribute__((always_inline)) INLINE void image_add(float* image_from,
                                                     float* image_to,
                                                     size_t size) {
  for (size_t i = 0; i < size; i++) {
    image_to[i] += image_from[i];
  }

  return;
}

/**
 * @brief Mapper function to turn part of the #parts into an image, stored
 *        in the extra data.
 *
 *        Creates a threadlocal image and smooths the given #parts onto it.
 *        At the end, we lock the main image and write back to main memory.
 *
 * @param map_data The array of #part.
 * @param num_parts The number of #parts.
 * @param extra_data Pointer to the image and some extra data.
 */
void create_projected_image_threadpool_mapper(void* map_data, int num_parts,
                                              void* extra_data) {

  /* TODO: Do this on a cell-by-cell basis and calculate geometry so
   * we need much smaller threadlocal allocations; at the moment this is
   * untenable for very large images. */

  /* TODO: correctly handle periodic boundary conditions */

  struct part* restrict parts = (struct part*)map_data;
  struct image_data* image_data = (struct image_data*)extra_data;

  /* First unpack everything to local variables */
  float* image = image_data->image;

  const size_t total_number_of_pixels =
      (size_t)image_data->image_properties.image_size[0] *
      image_data->image_properties.image_size[1];

  /* Create our own threadlocal image to avoid collision on the main
   * copy of the image in memory. */
  float* thread_image = (float*)malloc(total_number_of_pixels * sizeof(float));
  bzero(thread_image, total_number_of_pixels * sizeof(float));

  /* Parameters for the imaging */
  const float pixel_width_x = image_data->render_properties.box_size[0] /
                              image_data->image_properties.image_size[0];
  const float pixel_width_y = image_data->render_properties.box_size[1] /
                              image_data->image_properties.image_size[1];
  const float pixel_area = pixel_width_x * pixel_width_y;

  /* Conversion factors to number of pixels. Required as float / int comparisons
   * can cause some trickery */
  const float x_conversion_fac = 1.f / pixel_width_x;
  const float y_conversion_fac = 1.f / pixel_width_y;

  /* May cause some problems when dealing with non-square boxes */
  /* When to ignore the smoothing kernel and just go into MC mode */
  const float drop_to_single_cell = pixel_width_x * 0.5f;

  /* Perform a scatter over all particles */
  for (int particle = 0; particle < num_parts; particle++) {
    /* Grab the particle and extract properties! */
    const struct part* p = &parts[particle];

    /* Need to make sure the box is wrapped! Especially clear in e.g. a
     * Kelvin Helmholtz test. */
    /* TODO: Allow for offset x, y, z factors so that we can visualise any
     * arbritary area of the box, like py-sphviewer */
    const double x_wrapped =
        box_wrap(p->x[0], 0.0, image_data->render_properties.box_size[0]);
    const double y_wrapped =
        box_wrap(p->x[1], 0.0, image_data->render_properties.box_size[1]);

    /* What we really need is an integer position for comparison to the pixel
     * grid */
    const int x = x_wrapped * x_conversion_fac;
    const int y = y_wrapped * y_conversion_fac;
    /* This probably won't work for non-square boxes, we need to be smarter here
     */
    const float kernel_width =
        kernel_gamma * p->h / image_data->render_properties.box_size[0];

    if (kernel_width <= drop_to_single_cell) {
      /* Simple case; we average our mass over the size of the pixel */
      const int pixel = x + image_data->image_properties.image_size[0] * y;
      thread_image[pixel] += p->mass / pixel_area;
    } else {
      /* More complex; we need to SPH-smooth the data. This follows the python
       * version in swiftsimio very closely */

      /* May run into problems with non-square grids here... */
      const int cells_spanned =
          (int)(1.f +
                kernel_width * image_data->image_properties.image_size[0]);

      /* Ensure our loop bounds stay within the x, y grid */
      const int starting_x = max(0, x - cells_spanned);
      /* Remove 1 because we want the maximal _index_, and image_size is the
       * size of the image */
      const int ending_x = min(x + cells_spanned,
                               image_data->image_properties.image_size[0] - 1);
      const int starting_y = max(0, y - cells_spanned);
      const int ending_y = min(y + cells_spanned,
                               image_data->image_properties.image_size[1] - 1);

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
          const float density_addition = kernel_eval * p->mass;

          /* We can't write to the main image here; we need to create a
           * thread-local copy and then write the whole image at one whilst
           * the main image is locked at the end. */
          const int pixel =
              cell_x + image_data->image_properties.image_size[0] * cell_y;
          thread_image[pixel] += density_addition;
        }
      }
    }
  }

  /* Now we can merge the images from all threads, assuming we have the lock */
  if (lock_lock(&image_data->lock) == 0)
    image_add(thread_image, image, total_number_of_pixels);
  if (lock_unlock(&image_data->lock) != 0)
    error("Failed to unlock image_data.");

  /* Be a good thread and tidy up, will you? */
  free(thread_image);
}

/**
 * @brief Creates a projected image of all #parts that are in this
 *        copy of the #engine.
 *
 * @props e, pointer to the engine
 * @props image_data, pointer to image data.
 */
void create_projected_image(struct engine* e, struct image_data* image_data) {
  threadpool_map(&e->threadpool, create_projected_image_threadpool_mapper,
                 e->s->parts, e->s->nr_parts, sizeof(struct part), 0,
                 image_data);
}

/**
 * @brief dumps a projected image.
 *
 * @props e, pointer to the engine.
 **/
void image_dump_image(struct engine* e) {
  /* Be nice if we've been asked to be verboose */
  if (e->verbose) {
    if (e->policy & engine_policy_cosmology)
      message("Dumping an image at a=%e",
              exp(e->ti_current * e->time_base) * e->cosmology->a_begin);
    else
      message("Dumping an image at t=%e",
              e->ti_current * e->time_base + e->time_begin);
  }

#ifndef HAVE_HDF5

  /* There's no point making the image if we cannot dump it to file */
  message(
      "Warning: unable to create and dump image due to lack of HDF5 library. "
      "Please compile with HDF5 support to create image output.");

  /* Increment output count so we can continue on with the run */
  e->image_output_count++;
  return;

#else

  /* Extract properties from structs */
  const double box_size[2] = {e->s->dim[0], e->s->dim[1]};
  const int image_size[2] = {e->image_properties->image_size[0],
                             e->image_properties->image_size[1]};

  /* Allocate main image array */
  const size_t total_number_of_pixels = image_size[0] * image_size[1];
  float* image = (float*)malloc(total_number_of_pixels * sizeof(float));
  bzero(image, total_number_of_pixels * sizeof(float));

  /* Package this in our struct for passing to a threadpool */
  struct image_data image_data;
  image_data.image = image;
  /* There has to be a better way of doing these assignments */
  image_data.image_properties = *e->image_properties;
  image_data.render_properties.box_size[0] = box_size[0];
  image_data.render_properties.box_size[1] = box_size[1];

  /* Initialise the lock; we can't have everyone writing to the main image
   * all at once! */
  lock_init(&image_data.lock);

  /* Actually make the image! */
  create_projected_image(e, &image_data);

#ifdef WITH_MPI
  /* So we now have an image on each node, that we need to collapse. This is
   * fine, though, as we can just allow rank 0 to do all of the writes as they
   * should be (relatively) small */

  int myrank, res = 0;

  if ((res = MPI_Comm_rank(MPI_COMM_WORLD, &myrank)) != MPI_SUCCESS) {
    error("Call to MPI_Comm_rank failed with error %i.", res)
  }

  /* Sum everything into our local buffer image */
  MPI_Reduce(&image, &image, total_number_of_pixels, MPI_FLOAT, MPI_SUM, 0,
             MPI_COMM_WORLD);

  if (myrank == 0) {
  /* Only the root rank writes! */
#endif

    /* The following assumes that the image file already exists on disk. We
     * should have asserted this during the images_init call by creating an
     * empty HDF5 file. */

    hid_t h_file =
        H5Fopen(e->image_properties->image_filename, H5F_ACC_RDWR, H5P_DEFAULT);

    if (h_file < 0) {
      /* If we've failed we're having a very bad day */
      error("Failed to open image HDF5 file at %s.\n",
            e->image_properties->image_filename);
    }

    /* Do the 'if-exists-open else create' dance for our group */
    hid_t h_grp;
    if (H5Lexists(h_file, "/Images", H5P_DEFAULT)) {
      h_grp = H5Gopen(h_file, "/Images", H5P_DEFAULT);
    } else {
      h_grp =
          H5Gcreate(h_file, "/Images", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

    if (h_grp < 0) {
      error("Error while creating group /Images in %s.\n",
            e->image_properties->image_filename);
    }

    /* Now we want to create the datasets that we're going to use to store
     * our image, as well as write out some metadata properties. */

    char dataset_name[FILENAME_BUFFER_SIZE];
    snprintf(dataset_name, FILENAME_BUFFER_SIZE, "%s%04i",
             e->image_properties->image_base_name, e->image_output_count);

    /* Create a data space for our dataset to live in (a 2D square array) */
    hsize_t image_dims[2] = {(hsize_t)image_size[0], (hsize_t)image_size[1]};
    hid_t h_space = H5Screate_simple(2, image_dims, NULL);

    /* Create our dataset in which we'll store our pixel grid. If it already
     * exists, we delete the LINK to that dataset and start fresh (as it may
     * have a size that is incompatible with our own).
     *
     * Unfortunately, the HDF5 library as-of-yet does not provide a way to
     * actually remove data from a dataset, or reclaim that space. So this will
     * make the file larger than it needs to be. */
    if (H5Lexists(h_grp, dataset_name, H5P_DEFAULT)) {
      H5Ldelete(h_grp, dataset_name, H5P_DEFAULT);
    }

    /* Open up our fresh shiny dataset. */
    hid_t h_data = H5Dcreate(h_grp, dataset_name, H5T_NATIVE_FLOAT, h_space,
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /* Actually write out the data */
    herr_t write = H5Dwrite(h_data, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
                            H5P_DEFAULT, image);

    if (write < 0) {
      error("Failed to write image data to %s %s",
            e->image_properties->image_filename, dataset_name);
    }

    /* Now we can get started on our metadata! */
    const double output_time = e->ti_current * e->time_base;
    io_write_attribute_d(h_data, "Time", output_time);

    if (e->policy & engine_policy_cosmology) {
      io_write_attribute_d(h_data, "Scale-factor", e->cosmology->a);
      io_write_attribute_d(h_data, "Redshift", e->cosmology->z);
    }

    /* Close dataset properties */
    H5Dclose(h_data);
    H5Sclose(h_space);
    H5Gclose(h_grp);
    H5Fclose(h_file);

#ifdef WITH_MPI
    /* Close the bracket from only the root rank writing */
  }
#endif

  /* Increment the output count, otherwise we'll just overwrite our data next
   * time! */
  e->image_output_count++;

  free(image);

#endif /* HAVE_HDF5 */
}

void image_init(struct swift_params* params, const struct unit_system* us,
                const struct phys_const* phys_const,
                struct image_props* image_properties) {

  image_properties->image_size[0] =
      parser_get_param_int(params, "Images:image_size_x");
  image_properties->image_size[1] =
      parser_get_param_int(params, "Images:image_size_y");

  parser_get_opt_param_string(params, "Images:filename",
                              image_properties->image_filename, "images.hdf5");

  parser_get_opt_param_string(params, "Images:basename",
                              image_properties->image_base_name, "PixelGrid");

  /* We must create an empty HDF5 file for our images as we need to open it
   * and append later in the run. Unfortunately, this is the cleanest way of
   * doing this without resorting to nasty OS-dependent tricks.
   *
   * This also means that we overwrite all of our images every time we start
   * or restart a run, which perhaps isn't the best strategy. However, as we
   * do not call this function when restarting, it means that the images
   * HDF5 should be well preserved. */
  hid_t images_hdf5_file = H5Fcreate(image_properties->image_filename,
                                     H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  H5Fclose(images_hdf5_file);
}

/**
 * @brief Write an image_props struct to the given FILE as a stream of bytes.
 *
 * @param p the struct
 * @param stream the file stream
 */
void image_props_struct_dump(const struct image_props* p, FILE* stream) {
  restart_write_blocks((void*)p, sizeof(struct image_props), 1, stream,
                       "image props", "image function");
}

/**
 * @brief Restore a image_props struct from the given FILE as a stream of
 * bytes.
 *
 * @param p the struct
 * @param stream the file stream
 */
void image_props_struct_restore(const struct image_props* p, FILE* stream) {
  restart_read_blocks((void*)p, sizeof(struct image_props), 1, stream, NULL,
                      "image props");
}
