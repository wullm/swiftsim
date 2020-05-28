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

/* Array index (this is the row major format) */
inline int box_idx(int N, int x, int y, int z) { return z + N * (y + N * x); }

/* Row major index for a half-complex array with N*N*(N/2+1) complex entries */
inline int half_box_idx(int N, int x, int y, int z) {
  return z + (N / 2 + 1) * (y + N * x);
}

/* Sinc function */
inline double sinc(double x) { return x == 0 ? 0. : sin(x) / x; }

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

  /* Read the file name of the hdf5 gaussian random field file */
  char fieldFName[200] = "";
  parser_get_param_string(params, "Boltzmann:field_file_name", fieldFName);

  /* The file names of the perturbation data (either for reading or writing) */
  rend->in_perturb_fname = (char *)malloc(200 * sizeof(char));
  rend->out_perturb_fname = (char *)malloc(200 * sizeof(char));
  parser_get_opt_param_string(params, "Boltzmann:in_perturb_file_name",
                              rend->in_perturb_fname, "");
  parser_get_opt_param_string(params, "Boltzmann:out_perturb_file_name",
                              rend->out_perturb_fname, "perturb.hdf5");

  /* The file names of the CLASS parameter files */
  rend->class_ini_fname = (char *)malloc(200 * sizeof(char));
  rend->class_pre_fname = (char *)malloc(200 * sizeof(char));
  parser_get_opt_param_string(params, "Boltzmann:class_ini_file",
                              rend->class_ini_fname, "");
  parser_get_opt_param_string(params, "Boltzmann:class_pre_file",
                              rend->class_pre_fname, "");

  /* Open and load the file with the primordial Gaussian field */
  rend_load_primordial_field(rend, fieldFName);

  /* Print the loaded field dimensions */
  message(
      "Loaded %zu^3 primordial grid with dimensions: (%.1f, %.1f, %.1f) U_L "
      "on this node.",
      (size_t)rend->primordial_grid_N, rend->primordial_dims[0],
      rend->primordial_dims[1], rend->primordial_dims[2]);

  /* Verify that the physical dimensions of the primordial field match the
     cosmology */
  for (int i = 0; i < 3; i++) {
    if (fabs(rend->primordial_dims[i] - e->s->dim[i]) / e->s->dim[i] > 1e-3) {
      error("Dimensions[%i] of primordial field do not agree with space.", i);
    }
  }

  /* Verify that the primordial grid has the same size as the gravity mesh */
  if (rend->primordial_grid_N != (size_t)e->mesh->N) {
    error("Primordial grid is not the same size as the gravity mesh %zu!=%zu.",
          rend->primordial_grid_N, (size_t)e->mesh->N);
  }

  // /* Allocate memory for the rendered density grid */
  // const int N = e->mesh->N;
  // const int bytes = sizeof(double) * N * N * N;
  // rend->density_grid = (double *)swift_malloc("density_grid", bytes);
  //
  // if (rend->density_grid == NULL) {
  //   error("Error allocating memory for density grid.");
  // }
}

void rend_grids_alloc(struct renderer *rend) {
  /* Allocate memory for the perturbation theory grids */
  const int N = rend->primordial_grid_N;
  const size_t Nf = rend->transfer.n_functions;
  const size_t bytes = N * N * N * Nf * sizeof(double);
  rend->the_grids = (double *)swift_malloc("the_grids", bytes);

  if (rend->the_grids == NULL) {
    error("Error allocating memory for perturbation theory grids.");
  }

  /* Create pointer to the density grid, which is the first field in
     the array */
   rend->density_grid = rend->the_grids;
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
  const int bytes = N * N * N * sizeof(double);
  rend->primordial_grid = (double *)swift_malloc("primordial_grid", bytes);

  // Transfer the data to rend->primordial_grid
  for (size_t i = 0; i < N; i++) {
    for (size_t j = 0; j < N; j++) {
      for (size_t k = 0; k < N; k++) {
        rend->primordial_grid[k + j * N + i * N * N] = grf[i][j][k];
      }
    }
  }

  // Close the file
  H5Fclose(field_file);
}

void rend_interp_init(struct renderer *rend) {
#ifdef HAVE_LIBGSL

  /* The memory for the transfer functions is located here */
  struct transfer *tr = &rend->transfer;

  /* We will use bilinear interpolation in (tau, k) space */
  interp_type = gsl_interp2d_bilinear;

  /* Allocate memory for the spline */
  spline = gsl_spline2d_alloc(interp_type, tr->k_size, tr->tau_size);
  /* Note: this only copies the first transfer function from tr->delta */
  gsl_spline2d_init(spline, tr->k, tr->log_tau, tr->delta, tr->k_size,
                    tr->tau_size);

  /* Allocate memory for the accelerator objects */
  k_acc = gsl_interp_accel_alloc();
  tau_acc = gsl_interp_accel_alloc();
#else
  error("No GSL library found. Cannot perform cosmological interpolation.");
#endif
}

/* index_src is the index of the transfer function type */
void rend_interp_switch_source(struct renderer *rend, int index_src) {
#ifdef HAVE_LIBGSL

  /* The memory for the transfer functions is located here */
  struct transfer *tr = &rend->transfer;

  /* The array tr->delta contains a sequence of all transfer functions T(k,tau),
   * each of size tr->k_size * tr->tau_size doubles */
  int chunk_size = tr->k_size * tr->tau_size;

  /* Copy the desired transfer function to the spline */
  double *destination = spline->zarr;
  double *source_address = tr->delta + index_src * chunk_size;
  memcpy(destination, source_address, chunk_size * sizeof(double));

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

  /* Free function title strings */
  for (size_t i=0; i<rend->transfer.n_functions; i++) {
    free(rend->transfer.titles[i]);
  }
  free(rend->transfer.titles);

  /* Free density & perturbation theory grids */
  free(rend->the_grids);

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
  const struct cosmology *cosmo = e->cosmology;
  const double tau = cosmo->conformal_time;
  const double H_conformal = cosmo->H * cosmo->a;

  message("H-conformal = %f", H_conformal);

  /* Prevent out of interpolation range error */
  const int tau_size = rend->transfer.tau_size;
  const double final_log_tau = rend->transfer.log_tau[tau_size - 1];
  const double log_tau = min(log(tau), final_log_tau);

  /* Boxes in configuration and momentum space */
  double *restrict potential;
  fftw_complex *restrict fp;

  /* Allocate memory for the rendered field */
  potential = (double *)fftw_malloc(sizeof(double) * N * N * N);
  fp = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N * N * (N / 2 + 1));

  if (potential == NULL || fp == NULL) {
    error("Error allocating memory for rendering.");
  }

  memuse_log_allocation("potential", potential, 1, sizeof(double) * N * N * N);
  memuse_log_allocation("f", fp, 1, sizeof(fftw_complex) * N * N * (N / 2 + 1));

  /* Prepare the FFTW plans */
  fftw_plan pr2c = fftw_plan_dft_r2c_3d(N, N, N, potential, fp, FFTW_ESTIMATE);
  fftw_plan pc2r = fftw_plan_dft_c2r_3d(N, N, N, fp, potential, FFTW_ESTIMATE);

  /* Realize all the perturbation theory grids */
  for (size_t index_f = 0; index_f < rend->transfer.n_functions; index_f++) {
    /* Switch the interpolation spline to the desired transfer function */
    rend_interp_switch_source(rend, index_f);

    /* Use memory that has already been allocated */
    double *grid = rend->the_grids + index_f * N * N * N;

    /* First, copy the primordial field into the array */
    double *source_address = rend->primordial_grid;
    double *destination = grid;
    memcpy(destination, source_address, N * N * N * sizeof(double));

    /* Create plans */
    fftw_plan r2c_grid = fftw_plan_dft_r2c_3d(N, N, N, grid, fp, FFTW_ESTIMATE);
    fftw_plan c2r_grid = fftw_plan_dft_c2r_3d(N, N, N, fp, grid, FFTW_ESTIMATE);

    /* Transform to momentum space */
    fftw_execute(r2c_grid);

    /* Normalization */
    for (int i = 0; i < N * N * (N / 2 + 1); i++) {
      fp[i][0] *= box_volume / (N * N * N);
      fp[i][1] *= box_volume / (N * N * N);
    }

    /* Apply the transfer function */
    for (int x = 0; x < N; x++) {
      for (int y = 0; y < N; y++) {
        for (int z = 0; z <= N / 2; z++) {
          double k_x = (x > N / 2) ? (x - N) * delta_k : x * delta_k;  // U_L^-1
          double k_y = (y > N / 2) ? (y - N) * delta_k : y * delta_k;  // U_L^-1
          double k_z = (z > N / 2) ? (z - N) * delta_k : z * delta_k;  // U_L^-1

          double k = sqrt(k_x * k_x + k_y * k_y + k_z * k_z);

          /* Ignore the DC mode */
          if (k > 0) {
            double Tr = gsl_spline2d_eval(spline, k, log_tau, k_acc, tau_acc);

            /* The CIC Window function in Fourier space */
            double W_x = (k_x == 0) ? 1 : pow(sinc(0.5 * k_x * box_len / N), 2);
            double W_y = (k_y == 0) ? 1 : pow(sinc(0.5 * k_y * box_len / N), 2);
            double W_z = (k_z == 0) ? 1 : pow(sinc(0.5 * k_z * box_len / N), 2);
            double W = W_x * W_y * W_z;

            Tr /= W * W;

            fp[half_box_idx(N, x, y, z)][0] *= Tr;
            fp[half_box_idx(N, x, y, z)][1] *= Tr;
          } else {
            fp[half_box_idx(N, x, y, z)][0] = 0;
            fp[half_box_idx(N, x, y, z)][1] = 0;
          }
        }
      }
    }

    /* Transform back */
    fftw_execute(c2r_grid);

    /* Normalization */
    for (int i = 0; i < N * N * N; i++) {
      grid[i] /= box_volume;
    }

    /* Export the data block (only on master node) */
    if (e->nodeID == 0) {
      char boxname[40];
      sprintf(boxname, "grid_%zu.box", index_f);
      write_doubles_as_floats(boxname, grid, N * N * N);
    }
  }

  /* Next, compute the potential due to neutrinos (modulo a factor G_newt) */

   /* Calculate the background neutrino density at the present time */
   const double Omega_nu = cosmology_get_neutrino_density_param(cosmo, cosmo->a);
   const double rho_crit0 = cosmo->critical_density_0;
   const double neutrino_density = Omega_nu * rho_crit0;
   const double photon_density = cosmo->Omega_g * rho_crit0;
   const double ultra_relativistic_density = cosmo->Omega_ur * rho_crit0;

  /* The starting indices of the respective grids */
  double *ncdm_grid = rend->density_grid + rend->index_transfer_delta_ncdm * N * N * N;
  double *g_grid = rend->density_grid + rend->index_transfer_delta_g * N * N * N;
  double *ur_grid = rend->density_grid + rend->index_transfer_delta_ur * N * N * N;
  double *HT_prime_grid = rend->density_grid + rend->index_transfer_H_T_Nb_prime * N * N * N;
  double *HT_prime_prime_grid = rend->density_grid + rend->index_transfer_H_T_Nb_pprime * N * N * N;
  double *phi_grid = rend->density_grid + rend->index_transfer_phi * N * N * N;
  double *psi_grid = rend->density_grid + rend->index_transfer_psi * N * N * N;

  /* The potential is multiplied by G_newton later, so for phi, psi & H_T_Nb,
   * which should not be multiplied by G_newton, we need to divide now.
   * Also, we multiply by the scale factor a to get peculiar potentials.
   */
  const float G_newton = e->physical_constants->const_newton_G;
  const double potential_factor = cosmo->a / G_newton;

  /* Apply this factor to phi, psi, and the H_T_Nb derivatives */
  for (int i = 0; i < N * N * N; i++) {
    phi_grid[i] *= potential_factor;
    psi_grid[i] *= potential_factor;
    HT_prime_grid[i] *= potential_factor;
    HT_prime_prime_grid[i] *= potential_factor;
  }

  /* Compute RHS of Poisson's equation (modulo G_newton) */
  for (int i = 0; i < N * N * N; i++) {
    /* Neutrino contribution */
    double rho_ncdm = (1.0 + ncdm_grid[i]) * neutrino_density;
    /* Ultra-relativistic fluid contribution */
    double rho_ur = (1.0 + ur_grid[i]) * ultra_relativistic_density;
    /* Photon contribution (gamma) */
    double rho_g = (1.0 + g_grid[i]) * photon_density;
    /* H_T_Nb term = (H*a + d/dtau) * (d/dtau) * H_T_Nb */
    double H_T_term = HT_prime_grid[i] * H_conformal + HT_prime_prime_grid[i];

    /* We will apply the 1/k^2 kernel to the Fourier transform of this */
    potential[i] = -4 * M_PI * (rho_ncdm + rho_g + rho_ur) + H_T_term;

    /* Note: the phi & psi contribution is added later. */
  }

  /* Export the grids for troubleshooting */
  if (e->nodeID == 0) {
    write_doubles_as_floats("grid_ncdm.box", ncdm_grid, N * N * N);
    write_doubles_as_floats("grid_g.box", g_grid, N * N * N);
    write_doubles_as_floats("grid_ur.box", ur_grid, N * N * N);
    write_doubles_as_floats("grid_HT_prime.box", HT_prime_grid, N * N * N);
    write_doubles_as_floats("grid_HT_prime_prime.box", HT_prime_prime_grid, N * N * N);
    write_doubles_as_floats("gr_dens.box", potential, N * N * N);
    write_doubles_as_floats("grid_phi.box", phi_grid, N * N * N);
    write_doubles_as_floats("grid_psi.box", psi_grid, N * N * N);
  }

  /* Transform to momentum space */
  fftw_execute(pr2c);

  /* Normalization */
  for (int i = 0; i < N * N * (N / 2 + 1); i++) {
    fp[i][0] *= box_volume / (N * N * N);
    fp[i][1] *= box_volume / (N * N * N);
  }

  /* Multiply by the 1/k^2 kernel */
  for (int x = 0; x < N; x++) {
    for (int y = 0; y < N; y++) {
      for (int z = 0; z <= N / 2; z++) {
        double k_x = (x > N / 2) ? (x - N) * delta_k : x * delta_k;  // U_L^-1
        double k_y = (y > N / 2) ? (y - N) * delta_k : y * delta_k;  // U_L^-1
        double k_z = (z > N / 2) ? (z - N) * delta_k : z * delta_k;  // U_L^-1

        double k = sqrt(k_x * k_x + k_y * k_y + k_z * k_z);

        double kernel = 1.0 / k / k;

        /* Ignore the DC mode */
        if (k > 0) {
          fp[half_box_idx(N, x, y, z)][0] *= kernel;
          fp[half_box_idx(N, x, y, z)][1] *= kernel;
        } else {
          fp[half_box_idx(N, x, y, z)][0] = 0;
          fp[half_box_idx(N, x, y, z)][1] = 0;
        }
      }
    }
  }

  /* Transform back */
  fftw_execute(pc2r);

  /* Normalization */
  for (int i = 0; i < N * N * N; i++) {
    potential[i] /= box_volume;
  }

  /* Export the potentials */
  if (e->nodeID == 0) {
    write_doubles_as_floats("m_potential.box", e->mesh->potential, N * N * N);
    write_doubles_as_floats("gr_potential_without_stress.box", potential, N * N * N);
  }

  /* Add the contribution from anisotropic stress = (phi - psi) */
  for (int i = 0; i < N * N * N; i++) {
    potential[i] -= phi_grid[i] - psi_grid[i];
  }

  write_doubles_as_floats("gr_potential.box", potential, N * N * N);

  /* Add the contribution to the gravity mesh */
  for (int i = 0; i < N * N * N; i++) {
    e->mesh->potential[i] += potential[i];
  }

  write_doubles_as_floats("full_potential.box", e->mesh->potential, N * N * N);

  fftw_free(potential);
  fftw_free(fp);

#else
  error("No GSL library found. Cannot perform cosmological interpolation.");
#endif

#else
  error("No FFTW library found. Cannot compute periodic long-range forces.");
#endif
}

/* Read the perturbation data from a file */
void rend_read_perturb(struct renderer *rend, const struct engine *e,
                       char *fname) {
  /* The memory for the transfer functions is located here */
  struct transfer *tr = &rend->transfer;
  const struct unit_system *us = e->internal_units;

  hid_t h_file, h_grp, h_data, h_err, h_attr, h_tp;

  /* Open file */
  h_file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (h_file < 0) error("Error while opening file '%s'.", fname);

  message("Reading the perturbation from '%s'.", fname);

  /* Open header to read simulation properties */
  h_grp = H5Gopen(h_file, "/Header", H5P_DEFAULT);
  if (h_grp < 0) error("Error while opening file header\n");

  // tr->k_size;
  // tr->tau_size;
  // tr->n_functions;

  io_read_attribute(h_grp, "k_size", INT, &tr->k_size);
  io_read_attribute(h_grp, "tau_size", INT, &tr->tau_size);
  io_read_attribute(h_grp, "n_functions", INT, &tr->n_functions);

  /* Read the relevant units (length and time) */
  double file_length_us, file_time_us;
  io_read_attribute(h_grp, "Unit length in cgs (U_L)", DOUBLE, &file_length_us);
  io_read_attribute(h_grp, "Unit time in cgs (U_t)", DOUBLE, &file_time_us);

  message("Converting perturbation file units:");
  message("(perturb) Unit system: U_L = \t %.6e cm", file_length_us);
  message("(perturb) Unit system: U_T = \t %.6e s", file_time_us);
  message("to:");
  message("(internal) Unit system: U_L = \t %.6e cm", us->UnitLength_in_cgs);
  message("(internal) Unit system: U_T = \t %.6e s", us->UnitTime_in_cgs);

  /* Allocate memory for the transfer function titles */
  tr->titles = (char **)swift_calloc("titles", tr->n_functions, sizeof(char*));

  /* Read the titles of the transfer functions */
  h_attr = H5Aopen(h_grp, "FunctionTitles", H5P_DEFAULT);
  h_tp = H5Aget_type(h_attr);
  h_err = H5Aread(h_attr, h_tp, tr->titles);
  H5Aclose(h_attr);
  H5Tclose(h_tp);

  /* Print the titles */
  for (size_t i=0; i<tr->n_functions; i++) {
    message("Loaded perturbation vector '%s'.", tr->titles[i]);
  }

  /* Identify commonly used indices by their titles */
  for (size_t i=0; i<tr->n_functions; i++) {
    if (strcmp(tr->titles[i], "d_ncdm[0]") == 0) {
      rend->index_transfer_delta_ncdm = i;
      message("Identified ncdm density vector '%s'.", tr->titles[i]);
    } else if (strcmp(tr->titles[i], "d_g") == 0) {
      rend->index_transfer_delta_g = i;
      message("Identified photon density vector '%s'.", tr->titles[i]);
    } else if (strcmp(tr->titles[i], "d_ur") == 0) {
      rend->index_transfer_delta_ur = i;
      message("Identified ultra-relatistic fluid density vector '%s'.", tr->titles[i]);
    } else if (strcmp(tr->titles[i], "phi") == 0) {
      rend->index_transfer_phi = i;
      message("Identified scalar potential phi vector '%s'.", tr->titles[i]);
    } else if (strcmp(tr->titles[i], "psi") == 0) {
      rend->index_transfer_psi = i;
      message("Identified scalar potential psi vector '%s'.", tr->titles[i]);
    } else if (strcmp(tr->titles[i], "H_T_Nb_prime") == 0) {
      rend->index_transfer_H_T_Nb_prime = i;
      message("Identified N-body gauge H_T_prime vector '%s'.", tr->titles[i]);
    } else if (strcmp(tr->titles[i], "H_T_Nb_prime_prime") == 0) {
      rend->index_transfer_H_T_Nb_pprime = i;
      message("Identified N-body gauge H_T_prime_prime vector '%s'.", tr->titles[i]);
    }
  }

  /* Close header */
  H5Gclose(h_grp);

  // free(tr->k);
  // free(tr->log_tau);
  // free(tr->delta);

  tr->k = (double *)swift_calloc("k", tr->k_size, sizeof(double));
  tr->log_tau =
      (double *)swift_malloc("log_tau", tr->tau_size * sizeof(double));
  tr->delta = (double *)swift_malloc(
      "delta", tr->n_functions * tr->k_size * tr->tau_size * sizeof(double));

  message("We read the perturbation size %zu * %zu * %zu", tr->n_functions,
          tr->k_size, tr->tau_size);

  /* Open the perturbation data group */
  h_grp = H5Gopen(h_file, "/Perturb", H5P_DEFAULT);
  if (h_grp < 0) error("Error while opening perturbation group\n");

  /* Read the wavenumbers */
  h_data = H5Dopen2(h_grp, "Wavenumbers", H5P_DEFAULT);
  if (h_data < 0) error("Error while opening data space '%s'.", "Wavenumbers");

  h_err = H5Dread(h_data, io_hdf5_type(DOUBLE), H5S_ALL, H5S_ALL, H5P_DEFAULT,
                  tr->k);
  if (h_err < 0) error("Error while reading data array '%s'.", "Wavenumbers");

  /* Close the dataset */
  H5Dclose(h_data);

  /* Read the conformal times */
  h_data = H5Dopen2(h_grp, "Log conformal times", H5P_DEFAULT);
  if (h_data < 0)
    error("Error while opening data space '%s'.", "Log conformal times");

  h_err = H5Dread(h_data, io_hdf5_type(DOUBLE), H5S_ALL, H5S_ALL, H5P_DEFAULT,
                  tr->log_tau);
  if (h_err < 0)
    error("Error while reading data array '%s'.", "Log conformal times");

  /* Close the dataset */
  H5Dclose(h_data);

  /* Read the transfer functions */
  h_data = H5Dopen2(h_grp, "Transfer functions", H5P_DEFAULT);
  if (h_data < 0)
    error("Error while opening data space '%s'.", "Transfer functions");

  h_err = H5Dread(h_data, io_hdf5_type(DOUBLE), H5S_ALL, H5S_ALL, H5P_DEFAULT,
                  tr->delta);
  if (h_err < 0)
    error("Error while reading data array '%s'.", "Transfer functions");

  /* Convert the units of the wavenumbers (inverse length) */
  for (size_t i = 0; i < tr->k_size; i++) {
    tr->k[i] *= us->UnitLength_in_cgs / file_length_us;
  }

  /* Convert the units of the conformal time. This is log(tau) ! */
  for (size_t i = 0; i < tr->tau_size; i++) {
    tr->log_tau[i] += log(file_time_us) - log(us->UnitTime_in_cgs);
  }

  // for (size_t i=0; i<tr->k_size; i++) {
  //     printf("%e\n", tr->k[i]);
  // }
  //
  // for (size_t i=0; i<tr->tau_size; i++) {
  //     printf("%e\n", tr->log_tau[i]);
  // }
  //
  // for (size_t i=0; i<tr->k_size * tr->tau_size; i++) {
  //     printf("%e\n", tr->delta[i]);
  // }

  /* Close the dataset */
  H5Dclose(h_data);

  /* Close the perturbation group */
  H5Gclose(h_grp);

  /* Close file */
  H5Fclose(h_file);
}

/* Save the perturbation data to a file */
void rend_write_perturb(struct renderer *rend, const struct engine *e,
                        char *fname) {
  /* The memory for the transfer functions is located here */
  struct transfer *tr = &rend->transfer;
  const struct unit_system *us = e->internal_units;

  hid_t h_file, h_grp, h_data, h_err;

  /* Open file */
  h_file = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (h_file < 0) error("Error while opening file '%s'.", fname);

  message("Writing the perturbation to '%s'.", fname);

  /* Open header to write simulation properties */
  h_grp = H5Gcreate(h_file, "/Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (h_grp < 0) error("Error while creating file header\n");

  /* Print the relevant information and print status */
  io_write_attribute(h_grp, "k_size", INT, &tr->k_size, 1);
  io_write_attribute(h_grp, "tau_size", INT, &tr->tau_size, 1);
  io_write_attribute(h_grp, "n_functions", INT, &tr->n_functions, 1);

  /* Write unit system for this data set */
  io_write_attribute_d(h_grp, "Unit mass in cgs (U_M)",
                       units_get_base_unit(us, UNIT_MASS));
  io_write_attribute_d(h_grp, "Unit length in cgs (U_L)",
                       units_get_base_unit(us, UNIT_LENGTH));
  io_write_attribute_d(h_grp, "Unit time in cgs (U_t)",
                       units_get_base_unit(us, UNIT_TIME));
  io_write_attribute_d(h_grp, "Unit current in cgs (U_I)",
                       units_get_base_unit(us, UNIT_CURRENT));
  io_write_attribute_d(h_grp, "Unit temperature in cgs (U_T)",
                       units_get_base_unit(us, UNIT_TEMPERATURE));

  /* Close header */
  H5Gclose(h_grp);

  /* Open group to write the perturbation arrays */
  h_grp = H5Gcreate(h_file, "/Perturb", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (h_grp < 0) error("Error while creating perturbation group\n");

  /* Create data space */
  const hid_t h_space = H5Screate(H5S_SIMPLE);
  if (h_space < 0) error("Error while creating data space.");

  /* Set the extent of the data */
  int rank = 1;
  hsize_t shape[1] = {tr->k_size};
  h_err = H5Sset_extent_simple(h_space, rank, shape, shape);
  if (h_err < 0) error("Error while changing data space shape.");

  /* Dataset properties */
  const hid_t h_prop = H5Pcreate(H5P_DATASET_CREATE);

  /* Create dataset */
  h_data = H5Dcreate(h_grp, "Wavenumbers", io_hdf5_type(DOUBLE), h_space,
                     H5P_DEFAULT, h_prop, H5P_DEFAULT);
  if (h_data < 0) error("Error while creating dataspace '%s'.", "Wavenumbers");

  /* Write temporary buffer to HDF5 dataspace */
  h_err = H5Dwrite(h_data, io_hdf5_type(DOUBLE), h_space, H5S_ALL, H5P_DEFAULT,
                   tr->k);
  if (h_err < 0) error("Error while writing data array '%s'.", "tr->k");

  /* Close the dataset */
  H5Dclose(h_data);

  /* Set the extent of the tau data */
  rank = 1;
  hsize_t shape_tau[1] = {tr->tau_size};
  h_err = H5Sset_extent_simple(h_space, rank, shape_tau, shape_tau);
  if (h_err < 0) error("Error while changing data space shape.");

  /* Create dataset */
  h_data = H5Dcreate(h_grp, "Log conformal times", io_hdf5_type(DOUBLE),
                     h_space, H5P_DEFAULT, h_prop, H5P_DEFAULT);
  if (h_data < 0)
    error("Error while creating dataspace '%s'.", "Log conformal times");

  /* Write temporary buffer to HDF5 dataspace */
  h_err = H5Dwrite(h_data, io_hdf5_type(DOUBLE), h_space, H5S_ALL, H5P_DEFAULT,
                   tr->log_tau);
  if (h_err < 0) error("Error while writing data array '%s'.", "tr->log_tau");

  /* Close the dataset */
  H5Dclose(h_data);

  /* Set the extent of the transfer function data */
  rank = 3;
  hsize_t shape_delta[3] = {tr->n_functions, tr->k_size, tr->tau_size};
  h_err = H5Sset_extent_simple(h_space, rank, shape_delta, shape_delta);
  if (h_err < 0) error("Error while changing data space shape.");

  /* Create dataset */
  h_data = H5Dcreate(h_grp, "Transfer functions", io_hdf5_type(DOUBLE), h_space,
                     H5P_DEFAULT, h_prop, H5P_DEFAULT);
  if (h_data < 0)
    error("Error while creating dataspace '%s'.", "Transfer functions");

  /* Write temporary buffer to HDF5 dataspace */
  h_err = H5Dwrite(h_data, io_hdf5_type(DOUBLE), h_space, H5S_ALL, H5P_DEFAULT,
                   tr->delta);
  if (h_err < 0) error("Error while writing data array '%s'.", "tr->delta");

  /* Close the dataset */
  H5Dclose(h_data);

  /* Close the properties */
  H5Pclose(h_prop);

  /* Close the group */
  H5Gclose(h_grp);

  /* Close file */
  H5Fclose(h_file);
}

void rend_init_perturb_vec(struct renderer *rend, struct swift_params *params,
                           const struct engine *e, int myrank) {

    if (myrank == 0) {

      /* If a perturbation file & a CLASS ini file are both specified */
      if (strlen(rend->in_perturb_fname) > 1 && strlen(rend->class_ini_fname) > 1) {
        error("Specified both perturbation file & CLASS .ini file. '%s' '%s' (%ld)", rend->in_perturb_fname, rend->class_ini_fname,  strlen(rend->class_ini_fname));
      } else if (strlen(rend->in_perturb_fname) == '\0') {
#ifdef WITH_CLASS_INTERFACE
        /* Initialize perturbations to the cosmology with CLASS */
        message("We run CLASS to calculate perturbations to the cosmology.");
        rend_perturb_from_class(rend, params, e);
        message("Done with CLASS. The perturbations are now available.");

        /* Save to disk */
        if (strlen(rend->out_perturb_fname) > 1) {
          rend_write_perturb(rend, e, rend->out_perturb_fname);
        }
#else
        error("No CLASS library found. Cannot compute transfer functions.");
#endif
      } else if (strlen(rend->in_perturb_fname) > 1) {
        /* Read from disk */
        rend_read_perturb(rend, e, rend->in_perturb_fname);
      }

      /* Initialize our own interpolation spline */
      rend_interp_init(rend);
      rend_grids_alloc(rend);
    }

      /* Broadcast the cosmological perturbations to the other ranks */
#ifdef WITH_MPI
    /* The memory for the transfer functions is located here */
    struct transfer *tr = &rend->transfer;

    /* First broadcast the size of the perturbation to the other ranks */
    MPI_Bcast(&tr->k_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tr->tau_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tr->n_functions, 1, MPI_INT, 0, MPI_COMM_WORLD);

    /* Allocate memory on the other ranks */
    if (myrank != 0) {
      tr->delta =  (double *)swift_malloc("delta", tr->n_functions * tr->k_size * tr->tau_size * sizeof(double));
      tr->k = (double *)swift_malloc("k", tr->k_size * sizeof(double));
      tr->log_tau = (double *)swift_malloc("log_tau", tr->tau_size * sizeof(double));
    }

    /* Broadcast the perturbation to the other ranks */
    MPI_Bcast(tr->k, tr->k_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(tr->log_tau, tr->tau_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(tr->delta, tr->k_size * tr->tau_size * tr->n_functions, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /* Initialize the interpolation spline on the other ranks */
    if (myrank != 0) {
      rend_interp_init(rend);
      rend_grids_alloc(rend);
    }
#endif

}
