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

/* We also use the interpolation methods from the Eagle cooling codebase */
#include "../cooling/EAGLE/interpolate.h"

#include "class.h"

// Array index (this is the row major format)
inline int box_idx(int N, int x, int y, int z) { return z + N * (y + N * x); }

// Row major index for a half-complex array with N*N*(N/2+1) complex entries
inline int half_box_idx(int N, int x, int y, int z) {
  return z + (N / 2 + 1) * (y + N * x);
}

double rend_transfer(double k) { return exp(-k); }

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

  message("The primordial field has dimensions %f x %f x %f", field_dims[0],
          field_dims[1], field_dims[2]);
  message("The primordial grid has dimensions %zu x %zu x %zu",
          (size_t)grid_dims[0], (size_t)grid_dims[1], (size_t)grid_dims[2]);

  // Close the file
  H5Fclose(field_file);
}

void rend_add_to_mesh(struct renderer *rend, const struct engine *e) {
  // Grid size
  const size_t N = rend->primordial_grid_N;
  const double box_len = e->s->dim[0];
  const double box_volume = pow(box_len, 3);
  const double delta_k = 2 * M_PI / box_len;  // Mpc^-1

  double *restrict prime = (double *)fftw_malloc(sizeof(double) * N * N * N);
  if (prime == NULL) error("Error allocating memory for density mesh");
  // bzero(prime, N * N * N * sizeof(double));

  /* Allocates some memory for the mesh in Fourier space */
  fftw_complex *restrict fprime =
      (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N * N * (N / 2 + 1));
  if (fprime == NULL)
    error("Error allocating memory for transform of density mesh");
  memuse_log_allocation("fftw_fprime", fprime, 1,
                        sizeof(fftw_complex) * N * N * (N / 2 + 1));

  /* Prepare the FFT library */
  fftw_plan forward_plan = fftw_plan_dft_r2c_3d(
      N, N, N, prime, fprime, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
  fftw_plan inverse_plan = fftw_plan_dft_c2r_3d(
      N, N, N, fprime, prime, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);

  // Transfer the data to prime
  for (size_t i = 0; i < N; i++) {
    for (size_t j = 0; j < N; j++) {
      for (size_t k = 0; k < N; k++) {
        prime[box_idx(N, i, j, k)] = rend->primordial_grid[box_idx(N, i, j, k)];
      }
    }
  }

  // Transform to momentum space
  fftw_execute(forward_plan);

  // Normalization
  for (size_t x = 0; x < N; x++) {
    for (size_t y = 0; y < N; y++) {
      for (size_t z = 0; z <= N / 2; z++) {
        fprime[half_box_idx(N, x, y, z)][0] *= box_volume / (N * N * N);
        fprime[half_box_idx(N, x, y, z)][1] *= box_volume / (N * N * N);
      }
    }
  }

  // Apply the transfer function
  for (size_t x = 0; x < N; x++) {
    for (size_t y = 0; y < N; y++) {
      for (size_t z = 0; z <= N / 2; z++) {
        double k_x = (x > N / 2) ? (x - N) * delta_k : x * delta_k;  // Mpc^-1
        double k_y = (y > N / 2) ? (y - N) * delta_k : y * delta_k;  // Mpc^-1
        double k_z = (z > N / 2) ? (z - N) * delta_k : z * delta_k;  // Mpc^-1

        double k = sqrt(k_x * k_x + k_y * k_y + k_z * k_z);

        double T = rend_transfer(k);
        double TT = T * T;

        fprime[half_box_idx(N, x, y, z)][0] *= TT;
        fprime[half_box_idx(N, x, y, z)][1] *= TT;
      }
    }
  }

  // Transform back
  fftw_execute(inverse_plan);

  // Normalization
  for (size_t x = 0; x < N; x++) {
    for (size_t y = 0; y < N; y++) {
      for (size_t z = 0; z < N; z++) {
        prime[z + y * N + x * N * N] /= box_volume;
      }
    }
  }

  // Add the contribution to the gravity mesh
  for (size_t i = 0; i < N; i++) {
    for (size_t j = 0; j < N; j++) {
      for (size_t k = 0; k < N; k++) {
        e->mesh->potential[box_idx(N, i, j, k)] += prime[box_idx(N, i, j, k)];
      }
    }
  }

  message("Adding contributions to mesh.");
}

void rend_compute_perturbations(struct renderer *rend) {
    struct precision pr;  /* for precision parameters */
    struct background ba; /* for cosmological background */
    struct thermo th;     /* for thermodynamics */
    struct perturbs pt;   /* for source functions */
    struct transfers tr;  /* for transfer functions */
    struct primordial pm; /* for primordial spectra */
    struct spectra sp;    /* for output spectra */
    struct nonlinear nl;  /* for non-linear spectra */
    struct lensing le;    /* for lensed spectra */
    struct output op;     /* for output files */
    ErrorMsg errmsg;      /* for error messages */

    int class_argc = 2;
    char *class_argv[] = {"", "class_file.ini", NULL};

    if (input_init_from_arguments(class_argc, class_argv, &pr, &ba, &th, &pt, &tr,
                                  &pm, &sp, &nl, &le, &op, errmsg) == _FAILURE_) {
      error("Error running input_init_from_arguments \n=>%s\n", errmsg);
    }

    if (background_init(&pr, &ba) == _FAILURE_) {
      error("Error running background_init \n%s\n", ba.error_message);
    }

    if (thermodynamics_init(&pr, &ba, &th) == _FAILURE_) {
      error("Error in thermodynamics_init \n%s\n", th.error_message);
    }

    if (perturb_init(&pr, &ba, &th, &pt) == _FAILURE_) {
      error("Error in perturb_init \n%s\n", pt.error_message);
    }

    /* Try getting a source */
    int index_md = pt.index_md_scalars; //scalar mode
    int index_ic = 0; //index of the initial condition
    int index_tp = pt.index_tp_delta_ncdm1; //type of source function

    /* Size of the perturbations */
    int k_size = pt.k_size[index_md];
    int tau_size = pt.tau_size;

    /* Vector of the wavenumbers */
    rend->transfer.k_size = k_size;
    rend->transfer.k = (double*) calloc(k_size, sizeof(double));

    /* Vector of the conformal times at which the perturbation is sampled */
    rend->transfer.tau_size = tau_size;
    rend->transfer.tau = (double*) calloc(tau_size, sizeof(double));

    /* Vector with the transfer functions T(tau, k) */
    rend->transfer.delta = (double*) calloc(k_size*tau_size, sizeof(double));

    /* Read out the perturbation */
    for (int index_tau=0; index_tau<tau_size; index_tau++) {
        for (int index_k=0; index_k<k_size; index_k++) {
            double k = pt.k[index_md][index_k];
            double p = pt.sources[index_md][index_ic * pt.tp_size[index_md] + index_tp][index_tau * k_size + index_k];

            rend->transfer.k[index_k] = k;
            rend->transfer.delta[index_tau * k_size + index_k] = p;
        }
        rend->transfer.tau[index_tau] = pt.tau_sampling[index_tau];
    }

    message("The sizes are %i * %i", k_size, tau_size);

    /* Pre-empt segfault in CLASS if there is no interacting dark radiation */
    if (ba.has_idr == _FALSE_) {
        pt.alpha_idm_dr = (double*) malloc(0);
        pt.beta_idr = (double*) malloc(0);
    }

    /* Close CLASS again */
    if (perturb_free(&pt) == _FAILURE_) {
      error("Error in freeing class memory \n%s\n", pt.error_message);
    }

    if (thermodynamics_free(&th) == _FAILURE_) {
      error("Error in thermodynamics_free \n%s\n", th.error_message);
    }

    if (background_free(&ba) == _FAILURE_) {
      error("Error in background_free \n%s\n", ba.error_message);
    }
}
