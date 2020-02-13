/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Willem Elbers (whe@willemelbers.com)
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

/* Config parameters. */
#include "../../config.h"

/* This object's header. */
#include "renderer.h"

/* We use GSL for accelerated 2D interpolation */
#ifdef HAVE_LIBGSL
#include <gsl/gsl_spline2d.h>

/* GSL interpolation objects */
const gsl_interp2d_type *interp_type;
gsl_interp_accel *k_acc;
gsl_interp_accel *tau_acc;
gsl_spline2d *spline;
#endif

// Array index (this is the row major format)
inline int box_idx(int N, int x, int y, int z) { return z + N * (y + N * x); }

// Row major index for a half-complex array with N*N*(N/2+1) complex entries
inline int half_box_idx(int N, int x, int y, int z) {
  return z + (N / 2 + 1) * (y + N * x);
}

/* Quick and dirty write binary boxes */
inline void write_floats(char *fname, float *floats, int nfloats) {
  FILE *f = fopen(fname, "wb");
  fwrite(floats, sizeof(float), nfloats, f);
  fclose(f);
}

/* Quick and dirty write binary boxes */
inline void write_doubles_as_floats(char *fname, double *doubles, int nfloats) {
  /* Convert to floats */
  float *floats = (float *)malloc(sizeof(float) * nfloats);
  for (int i = 0; i < nfloats; i++) {
    floats[i] = (float)doubles[i];
  }

  FILE *f = fopen(fname, "wb");
  fwrite(floats, sizeof(float), nfloats, f);
  fclose(f);
  free(floats);
}

void rend_init(struct renderer *rend, struct swift_params *params,
               const struct engine *e) {
  // The user-specified number of k-bins used in the power spectrum calculation
  rend->num_of_k_bins =
      parser_get_opt_param_int(params, "Boltzmann:k_bins", BOLTZ_DEFAULT_BINS);

  // Read the file name of the hdf5 gaussian random field file
  char fieldFName[200] = "";
  parser_get_param_string(params, "Boltzmann:field_file_name", fieldFName);

  // Open and load the file with the primordial Gaussian field
  rend_load_primordial_field(rend, fieldFName);

  // Verify that the physical dimensions of the primordial field match the
  // cosmology
  for (int i = 0; i < 3; i++) {
    if (fabs(rend->primordial_dims[i] - e->s->dim[i]) / e->s->dim[i] > 1e-3) {
      error(
          "Dimensions[%i] of primordial field do not agree with engine->space.",
          i);
    }
  }

  // Verify that the primordial grid has the same grid size as the gravity mesh
  if (rend->primordial_grid_N != (size_t)e->mesh->N) {
    error(
        "Primordial grid is not the same size as the gravity mesh %zu != %zu.",
        rend->primordial_grid_N, (size_t)e->mesh->N);
  }
}

void rend_load_primordial_field(struct renderer *rend, const char *fname) {
  // Open the file containing the primordial fluctuation field
  const hid_t field_file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (field_file < 0) {
    error("Error opening the primordial field file.");
  }

  // Open the header group
  hid_t h_grp = H5Gopen(field_file, "/Header", H5P_DEFAULT);
  if (h_grp < 0) {
    error("Error while opening file header\n");
  }

  // Read the physical dimensions of the box
  const hid_t hid_bsz = H5Aexists(h_grp, "BoxSize");
  if (hid_bsz < 0) {
    error("Error while testing existance of 'BoxSize' attribute");
  }

  double field_dims[3];
  io_read_attribute(h_grp, "BoxSize", DOUBLE, field_dims);
  rend->primordial_dims = (double *)malloc(3 * sizeof(double));
  for (int i = 0; i < 3; i++) {
    rend->primordial_dims[i] = field_dims[i];
  }

  // Now load the actual grid
  h_grp = H5Gopen(field_file, "/Field", H5P_DEFAULT);
  if (h_grp < 0) {
    error("Error while opening field group\n");
  }

  hid_t h_data = H5Dopen(h_grp, "GaussianRandomField", H5P_DEFAULT);
  hid_t space = H5Dget_space(h_data);

  // The number of dimensions in the dataset (expected 3)
  const int rank = H5Sget_simple_extent_ndims(space);
  if (rank != 3) {
    error("Incorrect dataset dimensions for primordial field.");
  }

  // Find the extent of each dimension (the grid size; not physical size)
  hsize_t grid_dims[rank];
  H5Sget_simple_extent_dims(space, grid_dims, NULL);

  // The grid must be cubic
  if (grid_dims[0] != grid_dims[1] || grid_dims[0] != grid_dims[2]) {
    error("Primordial grid is not cubic.");
  }

  size_t N = grid_dims[0];
  rend->primordial_grid_N = N;

  // Create a temporary array to read the data
  float grf[N][N][N];
  H5Dread(h_data, H5T_NATIVE_FLOAT, space, space, H5P_DEFAULT, grf);
  H5Dclose(h_data);

  // Allocate memory in the main program
  rend->primordial_grid = malloc(N * N * N * sizeof(double));

  // Transfer the data to rend->primordial_grid
  for (size_t i = 0; i < N; i++) {
    for (size_t j = 0; j < N; j++) {
      for (size_t k = 0; k < N; k++) {
        rend->primordial_grid[k + j * N + i * N * N] = grf[i][j][k];
      }
    }
  }

  message("The primordial field has dimensions %fx%fx%f U_L^3", field_dims[0],
          field_dims[1], field_dims[2]);
  message("The primordial grid has dimensions %zu x %zu x %zu",
          (size_t)grid_dims[0], (size_t)grid_dims[1], (size_t)grid_dims[2]);

  // Close the file
  H5Fclose(field_file);
}

void rend_interp_init(struct renderer *rend) {
#ifdef HAVE_LIBGSL

  /* The memory for the transfer functions is located here */
  struct transfer *tr = &rend->transfer;

  /* We will use bilinear interpolation in (tau, k) space */
  interp_type = gsl_interp2d_bilinear;

  /* Initialize the spline */
  spline = gsl_spline2d_alloc(interp_type, tr->k_size, tr->tau_size);
  gsl_spline2d_init(spline, tr->k, tr->log_tau, tr->delta, tr->k_size,
                    tr->tau_size);

  /* Initialize accelerator objects */
  k_acc = gsl_interp_accel_alloc();
  tau_acc = gsl_interp_accel_alloc();

#else
  error("No GSL library found. Cannot perform cosmological interpolation.");
#endif
}

void rend_interp_free(struct renderer *rend) {
#ifdef HAVE_LIBGSL
  /* Done with the GSL interpolation */
  gsl_spline2d_free(spline);
  gsl_interp_accel_free(k_acc);
  gsl_interp_accel_free(tau_acc);
#else
  error("No GSL library found. Cannot perform cosmological interpolation.");
#endif
}

void rend_clean(struct renderer *rend) {
  /* Free the Gaussian field */
  free(rend->primordial_grid);
  free(rend->primordial_dims);

  /* Free interpolation tables */
  free(rend->transfer.delta);
  free(rend->transfer.k);
  free(rend->transfer.log_tau);

  /* Clean up the interpolation spline */
  rend_interp_free(rend);
}

void rend_add_to_mesh(struct renderer *rend, const struct engine *e) {
#ifdef HAVE_FFTW
#ifdef HAVE_LIBGSL
  /* Grid size */
  const int N = rend->primordial_grid_N;
  const double box_len = e->s->dim[0];
  const double box_volume = pow(box_len, 3);
  const double delta_k = 2 * M_PI / box_len;  // U_L^-1

  /* Current conformal time */
  double tau = e->cosmology->conformal_time;

  // cosmology_get_grav_kick_factor

  /* Boxes in configuration and momentum space */
  double *restrict prime;
  fftw_complex *restrict fp;

  /* Allocate memory for the rendered field */
  prime = (double *)fftw_malloc(sizeof(double) * N * N * N);
  fp = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N * N * (N / 2 + 1));

  if (prime == NULL || fp == NULL) {
    error("Error allocating memory for density mesh or its Fourier transform.");
  }

  memuse_log_allocation("prime", prime, 1, sizeof(fftw_complex) * N * N * N);
  memuse_log_allocation("f", fp, 1, sizeof(fftw_complex) * N * N * (N / 2 + 1));

  /* Prepare the FFTW plans */
  fftw_plan r2c = fftw_plan_dft_r2c_3d(N, N, N, prime, fp, FFTW_ESTIMATE);
  fftw_plan c2r = fftw_plan_dft_c2r_3d(N, N, N, fp, prime, FFTW_ESTIMATE);

  // Transfer the data to prime
  for (int i = 0; i < N * N * N; i++) {
    prime[i] = rend->primordial_grid[i];
  }

  // Transform to momentum space
  fftw_execute(r2c);

  // Normalization
  for (int i = 0; i < N * N * (N / 2 + 1); i++) {
    fp[i][0] *= box_volume / (N * N * N);
    fp[i][1] *= box_volume / (N * N * N);
  }

  // Apply the transfer function
  for (int x = 0; x < N; x++) {
    for (int y = 0; y < N; y++) {
      for (int z = 0; z <= N / 2; z++) {
        double k_x = (x > N / 2) ? (x - N) * delta_k : x * delta_k;  // U_L^-1
        double k_y = (y > N / 2) ? (y - N) * delta_k : y * delta_k;  // U_L^-1
        double k_z = (z > N / 2) ? (z - N) * delta_k : z * delta_k;  // U_L^-1

        double k = sqrt(k_x * k_x + k_y * k_y + k_z * k_z);

        /* Ignore the DC mode */
        if (k > 0) {
          double Tr = gsl_spline2d_eval(spline, k, log(tau), k_acc, tau_acc);

          fp[half_box_idx(N, x, y, z)][0] *= Tr * Tr;
          fp[half_box_idx(N, x, y, z)][1] *= Tr * Tr;
        }
      }
    }
  }

  // Transform back
  fftw_execute(c2r);

  // Normalization
  for (int i = 0; i < N * N * N; i++) {
    prime[i] /= box_volume;
  }

  write_doubles_as_floats("haha.box", prime, N * N * N);

  // Add the contribution to the gravity mesh
  for (int i = 0; i < N * N * N; i++) {
    e->mesh->potential[i] += prime[i];
  }

  message("Adding contributions to mesh.");

#else
  error("No GSL library found. Cannot perform cosmological interpolation.");
#endif

#else
  error("No FFTW library found. Cannot compute periodic long-range forces.");
#endif
}
