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
//#ifdef HAVE_LIBGSL
//#include <gsl/gsl_spline2d.h>
//#endif

/* Array index (this is the row major format) */
inline int box_idx(int N, int x, int y, int z) { return z + N * (y + N * x); }

/* Row major index for a half-complex array with N*N*(N/2+1) complex entries */
inline int half_box_idx(int N, int x, int y, int z) {
  return z + (N / 2 + 1) * (y + N * x);
}

/* Sinc function */
inline double sinc(double x) { return x == 0 ? 1. : sin(x) / x; }

/* Write binary boxes in HDF5 format */
int writeGRF_H5(const double *box, int N, double boxlen, const char *fname) {
    /* Create the hdf5 file */
    hid_t h_file = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /* Create the Header group */
    hid_t h_grp = H5Gcreate(h_file, "/Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /* Create dataspace for BoxSize attribute */
    const hsize_t arank = 1;
    const hsize_t adims[1] = {3}; //3D space
    hid_t h_aspace = H5Screate_simple(arank, adims, NULL);

    /* Create the BoxSize attribute and write the data */
    hid_t h_attr = H5Acreate1(h_grp, "BoxSize", H5T_NATIVE_DOUBLE, h_aspace, H5P_DEFAULT);
    double boxsize[3] = {boxlen, boxlen, boxlen};
    H5Awrite(h_attr, H5T_NATIVE_DOUBLE, boxsize);

    /* Close the attribute, corresponding dataspace, and the Header group */
    H5Aclose(h_attr);
    H5Sclose(h_aspace);
    H5Gclose(h_grp);

    /* Create the Field group */
    h_grp = H5Gcreate(h_file, "/Field", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /* Create dataspace for the field */
    const hsize_t frank = 3;
    const hsize_t fdims[3] = {N, N, N}; //3D space
    hid_t h_fspace = H5Screate_simple(frank, fdims, NULL);

    /* Create the dataset for the field */
    hid_t h_data = H5Dcreate(h_grp, "Field", H5T_NATIVE_DOUBLE, h_fspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /* Write the data */
    H5Dwrite(h_data, H5T_NATIVE_DOUBLE, h_fspace, h_fspace, H5P_DEFAULT, box);

    /* Close the dataset, corresponding dataspace, and the Field group */
    H5Dclose(h_data);
    H5Sclose(h_fspace);
    H5Gclose(h_grp);

    /* Close the file */
    H5Fclose(h_file);

    return 0;
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
  // rend->class_ini_fname = (char *)malloc(200 * sizeof(char));
  // rend->class_pre_fname = (char *)malloc(200 * sizeof(char));
  // parser_get_opt_param_string(params, "Boltzmann:class_ini_file",
  //                             rend->class_ini_fname, "");
  // parser_get_opt_param_string(params, "Boltzmann:class_pre_file",
  //                             rend->class_pre_fname, "");

  /* Open and load the file with the primordial Gaussian field */
  rend_load_primordial_field(rend, fieldFName);

  /* Print the loaded field dimensions */
  message(
      "Loaded %d^3 primordial grid with dimensions: (%.1f, %.1f, %.1f) U_L "
      "on this node.",
      rend->primordial_grid_N, rend->primordial_dims[0],
      rend->primordial_dims[1], rend->primordial_dims[2]);

  /* Verify that the physical dimensions of the primordial field match the
     cosmology */
  // for (int i = 0; i < 3; i++) {
  //   if (fabs(rend->primordial_dims[i] - e->s->dim[i]) / e->s->dim[i] > 1e-3) {
  //     error("Dimensions[%i] of primordial field do not agree with space.", i);
  //   }
  // }

  /* Verify that the primordial grid has the same size as the gravity mesh */
  if (rend->primordial_grid_N != e->mesh->N) {
    error("Primordial grid is not the same size as the gravity mesh %d!=%d.",
          rend->primordial_grid_N, e->mesh->N);
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
  const int Nf = rend->transfer.n_functions;
  const int bytes = N * N * N * Nf * sizeof(double);
  rend->the_grids = (double *)swift_malloc("the_grids", bytes);

  if (rend->the_grids == NULL) {
    error("Error allocating memory for perturbation theory grids.");
  }

  /* Create pointer to the density grid, which is the first field in
     the array */
  rend->density_grid = rend->the_grids;
}

void rend_load_primordial_field(struct renderer *rend, const char *fname) {
  /* Open the file containing the primordial fluctuation field */
  const hid_t field_file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (field_file < 0) {
    error("Error opening the primordial field file.");
  }

  /* Open the Header group */
  hid_t h_grp = H5Gopen(field_file, "/Header", H5P_DEFAULT);
  if (h_grp < 0) {
    error("Error while opening file header\n");
  }

  /* Determine whether the BoxSize attribute exists */
  const hid_t hid_bsz = H5Aexists(h_grp, "BoxSize");
  if (hid_bsz < 0) {
    error("Error while testing existance of 'BoxSize' attribute");
  }

  /* Load the physical BoxSize attribute */
  rend->primordial_dims = (double *)malloc(3 * sizeof(double));
  io_read_attribute(h_grp, "BoxSize", DOUBLE, rend->primordial_dims);

  /* Open the Field group */
  h_grp = H5Gopen(field_file, "/Field", H5P_DEFAULT);
  if (h_grp < 0) {
    error("Error while opening field group\n");
  }

  /* Open the dataset */
  hid_t h_data = H5Dopen(h_grp, "Field", H5P_DEFAULT);

  /* Open the dataspace */
  hid_t space = H5Dget_space(h_data);

  /* The number of dimensions in the dataset (expected 3) */
  const int rank = H5Sget_simple_extent_ndims(space);
  if (rank != 3) {
    error("Incorrect dataset dimensions for primordial field.");
  }

  /* Find the extent of each dimension (the grid size; not physical size) */
  hsize_t grid_dims[rank];
  H5Sget_simple_extent_dims(space, grid_dims, NULL);

  /* The grid must be cubic */
  if (grid_dims[0] != grid_dims[1] || grid_dims[0] != grid_dims[2]) {
    error("Primordial grid is not cubic.");
  }

  /* Store the grid dimensions */
  int N = grid_dims[0];
  rend->primordial_grid_N = N;

  /* Allocate memory for the Gaussian random field */
  const int bytes = N * N * N * sizeof(double);
  rend->primordial_grid = (double *)swift_malloc("primordial_grid", bytes);

  /* Read the data */
  H5Dread(h_data, H5T_NATIVE_DOUBLE, space, space, H5P_DEFAULT, rend->primordial_grid);

  /* Close the dataset */
  H5Dclose(h_data);

  /* Close the file */
  H5Fclose(field_file);
}

// void rend_interp_init(struct renderer *rend) {
// #ifdef HAVE_LIBGSL
//
//   /* The memory for the transfer functions is located here */
//   struct transfer *tr = &rend->transfer;
//
//   if (rend->spline == NULL) {
//     /* We will use bilinear interpolation in (tau, k) space */
//     rend->interp_type = gsl_interp2d_bilinear;
//
//     /* Allocate memory for the spline */
//     rend->spline =
//         gsl_spline2d_alloc(rend->interp_type, tr->k_size, tr->tau_size);
//     /* Note: this only copies the first transfer function from tr->delta */
//     gsl_spline2d_init(rend->spline, tr->k, tr->log_tau, tr->delta, tr->k_size,
//                       tr->tau_size);
//
//
//     /* Allocate memory for the accelerator objects */
//     rend->k_acc = gsl_interp_accel_alloc();
//     rend->tau_acc = gsl_interp_accel_alloc();
//   } else {
//    message("Warning, spline already allocated.");
//   }
// #else
//   error("No GSL library found. Cannot perform cosmological interpolation.");
// #endif
// }

/* index_src is the index of the transfer function type */
// void rend_interp_switch_source(struct renderer *rend, int index_src) {
// #ifdef HAVE_LIBGSL
//
//   /* The memory for the transfer functions is located here */
//   struct transfer *tr = &rend->transfer;
//
//   /* The array tr->delta contains a sequence of all transfer functions T(k,tau),
//    * each of size tr->k_size * tr->tau_size doubles */
//   int chunk_size = tr->k_size * tr->tau_size;
//
//   /* Copy the desired transfer function to the spline */
//   double *destination = rend->spline->zarr;
//   double *source_address = tr->delta + index_src * chunk_size;
//   memcpy(destination, source_address, chunk_size * sizeof(double));
//
// #else
//   error("No GSL library found. Cannot perform cosmological interpolation.");
// #endif
// }
//
// void rend_interp_free(struct renderer *rend) {
// #ifdef HAVE_LIBGSL
//   /* Done with the GSL interpolation */
//   gsl_spline2d_free(rend->spline);
//   gsl_interp_accel_free(rend->k_acc);
//   gsl_interp_accel_free(rend->tau_acc);
// #else
//   error("No GSL library found. Cannot perform cosmological interpolation.");
// #endif
// }

void rend_clean(struct renderer *rend) {
  /* Free the Gaussian field */
  free(rend->primordial_grid);
  free(rend->primordial_dims);

  /* Free strings */
  free(rend->in_perturb_fname);
  free(rend->out_perturb_fname);
  // free(rend->class_ini_fname);
  // free(rend->class_pre_fname);

  /* Free interpolation tables */
  free(rend->transfer.delta);
  free(rend->transfer.k);
  free(rend->transfer.log_tau);

  /* Free function title strings */
  for (int i = 0; i < rend->transfer.n_functions; i++) {
    free(rend->transfer.titles[i]);
  }
  free(rend->transfer.titles);

  /* Free density & perturbation theory grids */
// #ifdef RENDERER_FULL_GR
//   free(rend->the_grids);
// #endif

  /* Free the index table */
  free(rend->k_acc_table);

  /* Clean up the interpolation spline */
  // rend_interp_free(rend);
}

/* Add neutrinos using the Bird & Ali-Haïmoud method */
void rend_add_rescaled_nu_mesh(struct renderer *rend, const struct engine *e) {
#ifdef HAVE_FFTW
  /* Grid size */
  const int N = rend->primordial_grid_N;
  const double box_len = e->s->dim[0];
  const double box_volume = pow(box_len, 3);
  const double delta_k = 2 * M_PI / box_len;  // U_L^-1

  /* Current conformal time */
  const struct cosmology *cosmo = e->cosmology;
  const double tau = cosmo->conformal_time;

  /* Prevent out of interpolation range error */
  const int tau_size = rend->transfer.tau_size;
  const double final_log_tau = rend->transfer.log_tau[tau_size - 1];
  const double log_tau = min(log(tau), final_log_tau);

  /* Bilinear interpolation indices in (log_tau, k) space */
  int tau_index = 0, k_index = 0;
  double u_tau = 0.f, u_k = 0.f;

  /* Find the time index */
  rend_interp_locate_tau(rend, log_tau, &tau_index, &u_tau);

  /* Calculate the background neutrino density at the present time */
  const double Omega_nu = cosmology_get_neutrino_density_param(cosmo, cosmo->a);
  const double Omega_m = cosmo->Omega_m;  // does not include neutrinos
  /* The comoving density is (Omega_nu * a^-4) * a^3  = Omega_nu / a */
  const double bg_density_ratio = (Omega_nu / cosmo->a) / Omega_m;

  /* Long-range potential smoothing length */
  const double r_s = e->gravity_properties->a_smooth * box_len / N;

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

  /* Start a timer */
  ticks tic = getticks();

  /* First, copy the matter potential from SWIFT into the array */
  double *source_address = e->mesh->potential;
  double *destination = potential;
  memcpy(destination, source_address, N * N * N * sizeof(double));

  /* Transform to momentum space */
  fftw_execute(pr2c);

  /* Normalization */
  for (int i = 0; i < N * N * (N / 2 + 1); i++) {
    fp[i][0] *= box_volume / (N * N * N);
    fp[i][1] *= box_volume / (N * N * N);
  }

  if (e->verbose)
    message("Forward Fourier transform took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  /* Time the next stage */
  tic = getticks();

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
          /* Find the k-space interpolation index */
          rend_interp_locate_k(rend, k, &k_index, &u_k);

          /* Bilinear interpolation of the ncdm transfer function */
          double Trn = rend_custom_interp(rend, k_index, tau_index, u_tau, u_k,
                                          rend->index_transfer_delta_ncdm);

          /* Bilinear interpolation of the cdm transfer function */
          double Trc = rend_custom_interp(rend, k_index, tau_index, u_tau, u_k,
                                        rend->index_transfer_delta_cdm);

          /* The long-range kernel */
          double K = 1.;
          fourier_kernel_long_grav_eval(k * k * r_s * r_s, &K);

          fp[half_box_idx(N, x, y, z)][0] *= (Trn / Trc) * bg_density_ratio / K;
          fp[half_box_idx(N, x, y, z)][1] *= (Trn / Trc) * bg_density_ratio / K;
        } else {
          fp[half_box_idx(N, x, y, z)][0] = 0;
          fp[half_box_idx(N, x, y, z)][1] = 0;
        }
      }
    }
  }

  if (e->verbose)
    message("Applying transfer function took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  /* Time the next stage */
  tic = getticks();

  /* Transform back */
  fftw_execute(pc2r);

  /* Normalization */
  for (int i = 0; i < N * N * N; i++) {
    potential[i] /= box_volume;
  }

  if (e->verbose)
    message("Backward Fourier transform took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  /* Time the next stage */
  tic = getticks();

  /* Export the potentials if necessary */
  if (e->nodeID == 0) {
    /* Are we also exporting snapshots? */
    if (e->step % 50 == 0) {
      char one[40];
      char two[40];
      double z = e->cosmology->z;
      sprintf(one, "m_potential_z_%.2f.hdf5", z);
      sprintf(two, "scaled_nu_potential_z_%.2f.hdf5", z);
      writeGRF_H5(e->mesh->potential, N, box_len, one);
      writeGRF_H5(potential, N, box_len, two);

      if (e->verbose) {
        /* Print some statistics */
        double rms_matter = 0.f;
        double rms_nu = 0.f;

        for (int i = 0; i < N * N * N; i++) {
          rms_matter += e->mesh->potential[i] * e->mesh->potential[i];
          rms_nu += potential[i] * potential[i];
        }

        message("Dumping render boxes took %.3f %s.",
                clocks_from_ticks(getticks() - tic), clocks_getunit());
        message("[Phi_m, Phi_nu] = [%e, %e].", rms_matter, rms_nu);
      }
    }
  }

  /* Add the contribution to the gravity mesh */
  for (int i = 0; i < N * N * N; i++) {
    e->mesh->potential[i] += potential[i];
  }

  /* Free memory */
  fftw_free(potential);
  fftw_free(fp);
  fftw_destroy_plan(pr2c);
  fftw_destroy_plan(pc2r);

#else
  error("No FFTW library found. Cannot compute periodic long-range forces.");
#endif
}

/* Add neutrinos using the linear theory transfer function */
void rend_add_linear_nu_mesh(struct renderer *rend, const struct engine *e) {
#ifdef HAVE_FFTW
  /* Grid size */
  const int N = rend->primordial_grid_N;
  const double box_len = e->s->dim[0];
  const double box_volume = pow(box_len, 3);
  const double delta_k = 2 * M_PI / box_len;  // U_L^-1

  /* Current conformal time */
  const struct cosmology *cosmo = e->cosmology;
  const double tau = cosmo->conformal_time;

  /* Prevent out of interpolation range error */
  const int tau_size = rend->transfer.tau_size;
  const double final_log_tau = rend->transfer.log_tau[tau_size - 1];
  const double log_tau = min(log(tau), final_log_tau);

  /* Throw an error if the neutrino transfer function was not found */
  if (rend->index_transfer_delta_ncdm < 0) {
    error("Running the neutrino renderer: transfer function d_ncdm is absent.");
  }

  /* Bilinear interpolation indices in (log_tau, k) space */
  int tau_index = 0, k_index = 0;
  double u_tau = 0.f, u_k = 0.f;

  /* Find the time index */
  rend_interp_locate_tau(rend, log_tau, &tau_index, &u_tau);

  /* Calculate the background neutrino density at the current time step */
  const double Omega_nu = cosmology_get_neutrino_density_param(cosmo, cosmo->a);
  const double rho_crit0 = cosmo->critical_density_0;
  /* The comoving density is (Omega_nu * a^-4) * a^3  = Omega_nu / a */
  const double neutrino_density = Omega_nu * rho_crit0 / cosmo->a;

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

  /* Start a timer */
  ticks tic = getticks();

  /* First, copy the primordial field into the array */
  double *source_address = rend->primordial_grid;
  double *destination = potential;
  memcpy(destination, source_address, N * N * N * sizeof(double));

  /* Transform to momentum space */
  fftw_execute(pr2c);

  /* Normalization */
  for (int i = 0; i < N * N * (N / 2 + 1); i++) {
    fp[i][0] *= box_volume / (N * N * N);
    fp[i][1] *= box_volume / (N * N * N);
  }

  if (e->verbose)
    message("Forward Fourier transform took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  /* Time the next stage */
  tic = getticks();

  /* Apply the neutrino transfer function */
  for (int x = 0; x < N; x++) {
    for (int y = 0; y < N; y++) {
      for (int z = 0; z <= N / 2; z++) {
        double k_x = (x > N / 2) ? (x - N) * delta_k : x * delta_k;  // U_L^-1
        double k_y = (y > N / 2) ? (y - N) * delta_k : y * delta_k;  // U_L^-1
        double k_z = (z > N / 2) ? (z - N) * delta_k : z * delta_k;  // U_L^-1

        double k = sqrt(k_x * k_x + k_y * k_y + k_z * k_z);

        /* Ignore the DC mode */
        if (k > 0) {
          /* Find the k-space interpolation index */
          rend_interp_locate_k(rend, k, &k_index, &u_k);

          /* Bilinear interpolation of the ncdm transfer function */
          double Tr = rend_custom_interp(rend, k_index, tau_index, u_tau, u_k,
                                         rend->index_transfer_delta_ncdm);


          /* The CIC Window function in Fourier space */
          const double sqrt_W_x = sinc(0.5 * k_x * box_len / N);
          const double sqrt_W_y = sinc(0.5 * k_y * box_len / N);
          const double sqrt_W_z = sinc(0.5 * k_z * box_len / N);
          const double sqrt_W = sqrt_W_x * sqrt_W_y * sqrt_W_z;
          const double W = sqrt_W * sqrt_W;

          /* Only deconvolve once for the interpolation (no gpart assignment) */
          Tr /= W;

          /* Convert from overdensity to density (we can ignore the k=0 mode) */
          Tr *= neutrino_density;

          /* Convert from density to potential by applying the -4 pi/k^2 kernel */
          Tr *= -4 * M_PI / k / k;

          fp[half_box_idx(N, x, y, z)][0] *= Tr;
          fp[half_box_idx(N, x, y, z)][1] *= Tr;
        } else {
          fp[half_box_idx(N, x, y, z)][0] = 0;
          fp[half_box_idx(N, x, y, z)][1] = 0;
        }
      }
    }
  }

  if (e->verbose)
    message("Applying transfer function took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  /* Time the next stage */
  tic = getticks();

  /* Transform back */
  fftw_execute(pc2r);

  /* Normalization */
  for (int i = 0; i < N * N * N; i++) {
    potential[i] /= box_volume;
  }

  if (e->verbose)
    message("Backward Fourier transform took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  /* Time the next stage */
  tic = getticks();

  /* Export the potentials if necessary */
  if (e->nodeID == 0) {
    /* Are we also exporting snapshots? */
    if (e->step % 50 == 0) {
      char one[40];
      char two[40];
      double z = e->cosmology->z;
      sprintf(one, "m_potential_z_%.2f.hdf5", z);
      sprintf(two, "linear_nu_potential_z_%.2f.hdf5", z);
      writeGRF_H5(e->mesh->potential, N, box_len, one);
      writeGRF_H5(potential, N, box_len, two);

      if (e->verbose) {
        /* Print some statistics */
        double rms_matter = 0.f;
        double rms_nu = 0.f;

        for (int i = 0; i < N * N * N; i++) {
          rms_matter += e->mesh->potential[i] * e->mesh->potential[i];
          rms_nu += potential[i] * potential[i];
        }

        message("Dumping render boxes took %.3f %s.",
                clocks_from_ticks(getticks() - tic), clocks_getunit());
        message("[Phi_m, Phi_nu] = [%e, %e].", rms_matter, rms_nu);
      }
    }
  }

  /* Add the contribution to the gravity mesh */
  for (int i = 0; i < N * N * N; i++) {
    e->mesh->potential[i] += potential[i];
  }

  /* Free memory */
  fftw_free(potential);
  fftw_free(fp);
  fftw_destroy_plan(pr2c);
  fftw_destroy_plan(pc2r);

#else
  error("No FFTW library found. Cannot compute periodic long-range forces.");
#endif
}



/* Add all the linear theory transfer functions for relativistic components */
void rend_add_gr_potential_mesh(struct renderer *rend, const struct engine *e) {
#ifdef HAVE_FFTW
  /* Grid size */
  const int N = rend->primordial_grid_N;
  const double box_len = e->s->dim[0];
  const double box_volume = pow(box_len, 3);
  const double delta_k = 2 * M_PI / box_len;  // U_L^-1

  /* Current conformal time */
  const struct cosmology *cosmo = e->cosmology;
  const double tau = cosmo->conformal_time;
  const double H_conformal = cosmo->H * cosmo->a;

  /* Prevent out of interpolation range error */
  const int tau_size = rend->transfer.tau_size;
  const double final_log_tau = rend->transfer.log_tau[tau_size - 1];
  const double log_tau = min(log(tau), final_log_tau);

  /* Throw an error if essential transfer functions are missing */
  if (rend->index_transfer_H_T_Nb_prime < 0) {
    error("Running the GR renderer: transfer function HT_Nb_prime is absent.");
  }
  if (rend->index_transfer_phi < 0) {
    error("Running the GR renderer: transfer function phi is absent.");
  }
  if (rend->index_transfer_psi < 0) {
    error("Running the GR renderer: transfer function psi is absent.");
  }

  /* Bilinear interpolation indices in (log_tau, k) space */
  int tau_index = 0, k_index = 0;
  double u_tau = 0.f, u_k = 0.f;

  /* Find the time index */
  rend_interp_locate_tau(rend, log_tau, &tau_index, &u_tau);

  /* Calculate the background neutrino density at the present time */
  // const double Omega_nu = cosmology_get_neutrino_density_param(cosmo, cosmo->a);
  const double Omega_g = cosmo->Omega_g;
  const double Omega_ur = cosmo->Omega_ur;
  const double rho_crit0 = cosmo->critical_density_0;
  /* The comoving density is (Omega_nu * a^-4) * a^3  = Omega_nu / a */
  // const double neutrino_density = Omega_nu * rho_crit0 / cosmo->a;
  const double photon_density = Omega_g * rho_crit0 / cosmo->a;
  const double ultra_relativistic_density = Omega_ur * rho_crit0 / cosmo->a;

  /* The potential is multiplied by G_newton later, so for phi, psi & H_T_Nb,
   * which should not be multiplied by G_newton, we need to divide here.
   * Also, we multiply by the scale factor a to get peculiar potentials.
   */
  const float G_newton = e->physical_constants->const_newton_G;
  const double potential_factor = cosmo->a / G_newton;

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

  /* Start a timer */
  ticks tic = getticks();

  /* First, copy the primordial field into the array */
  double *source_address = rend->primordial_grid;
  double *destination = potential;
  memcpy(destination, source_address, N * N * N * sizeof(double));

  /* Transform to momentum space */
  fftw_execute(pr2c);

  /* Normalization */
  for (int i = 0; i < N * N * (N / 2 + 1); i++) {
    fp[i][0] *= box_volume / (N * N * N);
    fp[i][1] *= box_volume / (N * N * N);
  }

  if (e->verbose)
    message("Forward Fourier transform took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  /* Time the next stage */
  tic = getticks();

  /* Apply the neutrino transfer function */
  for (int x = 0; x < N; x++) {
    for (int y = 0; y < N; y++) {
      for (int z = 0; z <= N / 2; z++) {
        double k_x = (x > N / 2) ? (x - N) * delta_k : x * delta_k;  // U_L^-1
        double k_y = (y > N / 2) ? (y - N) * delta_k : y * delta_k;  // U_L^-1
        double k_z = (z > N / 2) ? (z - N) * delta_k : z * delta_k;  // U_L^-1

        double k = sqrt(k_x * k_x + k_y * k_y + k_z * k_z);

        /* Ignore the DC mode */
        if (k > 0) {
          /* Find the k-space interpolation index */
          rend_interp_locate_k(rend, k, &k_index, &u_k);

          /* Interpolate relativistic species density transfer functions */
          // double Tr_nu = 0;
          double Tr_g = 0;
          double Tr_ur = 0;

          /* Bilinear interpolation of the ncdm transfer function */
          // if (rend->index_transfer_delta_ncdm > -1) {
          //   Tr_nu = rend_custom_interp(rend, k_index, tau_index, u_tau, u_k,
          //                              rend->index_transfer_delta_ncdm);
          // }

          /* Bilinear interpolation of the photon transfer function */
          if (rend->index_transfer_delta_g > -1) {
            Tr_g = rend_custom_interp(rend, k_index, tau_index, u_tau, u_k,
                                      rend->index_transfer_delta_g);
          }

          /* Bilinear interpolation of the ur transfer function */
          if (rend->index_transfer_delta_ur > -1) {
            Tr_ur = rend_custom_interp(rend, k_index, tau_index, u_tau, u_k,
                                       rend->index_transfer_delta_ur);
          }

          /* Convert from overdensity to density (we ignore the k=0 mode) */
          // Tr_nu *= neutrino_density;
          Tr_g *= photon_density;
          Tr_ur *= ultra_relativistic_density;

          /* Collect the relativistic fluid contributions to the potential */
          double Tr_pot = -4 * M_PI * (Tr_g + Tr_ur);

          /* Interpolate metric derivative transfer functions */
          double Tr_HT_p = 0;
          double Tr_HT_pp = 0;

          /* Bilinear interpolation of the metric derivative functions */
          if (rend->index_transfer_H_T_Nb_prime > -1) {
            Tr_HT_p = rend_custom_interp(rend, k_index, tau_index, u_tau, u_k,
                                         rend->index_transfer_H_T_Nb_prime);
          }
          if (rend->index_transfer_H_T_Nb_pprime > -1) {
            Tr_HT_pp = rend_custom_interp(rend, k_index, tau_index, u_tau, u_k,
                                          rend->index_transfer_H_T_Nb_pprime);
          }

          /* Compute the contributiom from the transverse metric term H_T */
          double Tr_HT_term = Tr_HT_p * H_conformal + Tr_HT_pp;

          /* Interpolate scalar metric transfer functions */
          double Tr_phi = 0;
          double Tr_psi = 0;

          /* Bilinear interpolation of the scalar metric transfer functions */
          if (rend->index_transfer_phi > -1) {
            Tr_phi = rend_custom_interp(rend, k_index, tau_index, u_tau, u_k,
                                        rend->index_transfer_phi);
          }
          if (rend->index_transfer_psi > -1) {
            Tr_psi = rend_custom_interp(rend, k_index, tau_index, u_tau, u_k,
                                        rend->index_transfer_psi);
          }

          /* The anisotropic stress term */
          double Tr_as_term = Tr_psi - Tr_phi;

          /* Divide Newton's constant out of the potential terms */
          Tr_HT_term *= potential_factor;
          Tr_as_term *= potential_factor;

          /* Add all the contributions */
          double Tr = (Tr_pot + Tr_HT_term) / (k * k) + Tr_as_term;

          /* The CIC Window function in Fourier space */
          const double sqrt_W_x = sinc(0.5 * k_x * box_len / N);
          const double sqrt_W_y = sinc(0.5 * k_y * box_len / N);
          const double sqrt_W_z = sinc(0.5 * k_z * box_len / N);
          const double sqrt_W = sqrt_W_x * sqrt_W_y * sqrt_W_z;
          const double W = sqrt_W * sqrt_W;

          /* Only deconvolve once for the interpolation (no gpart assignment) */
          Tr /= W;

          fp[half_box_idx(N, x, y, z)][0] *= Tr;
          fp[half_box_idx(N, x, y, z)][1] *= Tr;
        } else {
          fp[half_box_idx(N, x, y, z)][0] = 0;
          fp[half_box_idx(N, x, y, z)][1] = 0;
        }
      }
    }
  }

  if (e->verbose)
    message("Applying transfer function took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  /* Time the next stage */
  tic = getticks();

  /* Transform back */
  fftw_execute(pc2r);

  /* Normalization */
  for (int i = 0; i < N * N * N; i++) {
    potential[i] /= box_volume;
  }

  if (e->verbose)
    message("Backward Fourier transform took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  /* Time the next stage */
  tic = getticks();

  /* Export the potentials if necessary */
  if (e->nodeID == 0) {
    /* Are we also exporting snapshots? */
    if (e->step % 500 == 0) {
      char one[40];
      char two[40];
      double z = e->cosmology->z;
      sprintf(one, "m_potential_z_%.2f.hdf5", z);
      sprintf(two, "linear_gr_potential_z_%.2f.hdf5", z);
      writeGRF_H5(e->mesh->potential, N, box_len, one);
      writeGRF_H5(potential, N, box_len, two);

      if (e->verbose) {
        /* Print some statistics */
        double rms_matter = 0.f;
        double rms_nu = 0.f;

        for (int i = 0; i < N * N * N; i++) {
          rms_matter += e->mesh->potential[i] * e->mesh->potential[i];
          rms_nu += potential[i] * potential[i];
        }

        message("Dumping render boxes took %.3f %s.",
                clocks_from_ticks(getticks() - tic), clocks_getunit());
        message("[Phi_m, Phi_gr] = [%e, %e].", rms_matter, rms_nu);
      }
    }
  }

  /* Add the contribution to the gravity mesh */
  for (int i = 0; i < N * N * N; i++) {
    e->mesh->potential[i] += potential[i];
  }

  /* Free memory */
  fftw_free(potential);
  fftw_free(fp);
  fftw_destroy_plan(pr2c);
  fftw_destroy_plan(pc2r);

#else
  error("No FFTW library found. Cannot compute periodic long-range forces.");
#endif
}

// void rend_add_gr_potential_mesh(struct renderer *rend, const struct engine *e) {
// #ifdef HAVE_FFTW
// #ifdef HAVE_LIBGSL
//   /* Grid size */
//   const int N = rend->primordial_grid_N;
//   const double box_len = e->s->dim[0];
//   const double box_volume = pow(box_len, 3);
//   const double delta_k = 2 * M_PI / box_len;  // U_L^-1
//
//   /* Current conformal time */
//   const struct cosmology *cosmo = e->cosmology;
//   const double tau = cosmo->conformal_time;
//   const double H_conformal = cosmo->H * cosmo->a;
//
//   message("H-conformal = %f", H_conformal);
//
//   /* Prevent out of interpolation range error */
//   const int tau_size = rend->transfer.tau_size;
//   const double final_log_tau = rend->transfer.log_tau[tau_size - 1];
//   const double log_tau = min(log(tau), final_log_tau);
//
//   /* Bilinear interpolation indices in (log_tau, k) space */
//   int tau_index = 0, k_index = 0;
//   double u_tau = 0.f, u_k = 0.f;
//
//   /* Find the time index */
//   rend_interp_locate_tau(rend, log_tau, &tau_index, &u_tau);
//
//   /* Boxes in configuration and momentum space */
//   double *restrict potential;
//   fftw_complex *restrict fp;
//
//   /* Allocate memory for the rendered field */
//   potential = (double *)fftw_malloc(sizeof(double) * N * N * N);
//   fp = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N * N * (N / 2 + 1));
//
//   if (potential == NULL || fp == NULL) {
//     error("Error allocating memory for rendering.");
//   }
//
//   memuse_log_allocation("potential", potential, 1, sizeof(double) * N * N * N);
//   memuse_log_allocation("f", fp, 1, sizeof(fftw_complex) * N * N * (N / 2 + 1));
//
//   /* Prepare the FFTW plans */
//   fftw_plan pr2c = fftw_plan_dft_r2c_3d(N, N, N, potential, fp, FFTW_ESTIMATE);
//   fftw_plan pc2r = fftw_plan_dft_c2r_3d(N, N, N, fp, potential, FFTW_ESTIMATE);
//
//   /* Realize all the perturbation theory grids */
//   for (int index_f = 0; index_f < rend->transfer.n_functions; index_f++) {
//     /* Switch the interpolation spline to the desired transfer function */
//     // rend_interp_switch_source(rend, index_f);
//
//     /* Use memory that has already been allocated */
//     double *grid = rend->the_grids + index_f * N * N * N;
//
//     /* First, copy the primordial field into the array */
//     double *source_address = rend->primordial_grid;
//     double *destination = grid;
//     memcpy(destination, source_address, N * N * N * sizeof(double));
//
//     /* Create plans */
//     fftw_plan r2c_grid = fftw_plan_dft_r2c_3d(N, N, N, grid, fp, FFTW_ESTIMATE);
//     fftw_plan c2r_grid = fftw_plan_dft_c2r_3d(N, N, N, fp, grid, FFTW_ESTIMATE);
//
//     /* Transform to momentum space */
//     fftw_execute(r2c_grid);
//
//     /* Normalization */
//     for (int i = 0; i < N * N * (N / 2 + 1); i++) {
//       fp[i][0] *= box_volume / (N * N * N);
//       fp[i][1] *= box_volume / (N * N * N);
//     }
//
//     /* Apply the transfer function */
//     for (int x = 0; x < N; x++) {
//       for (int y = 0; y < N; y++) {
//         for (int z = 0; z <= N / 2; z++) {
//           double k_x = (x > N / 2) ? (x - N) * delta_k : x * delta_k;  // U_L^-1
//           double k_y = (y > N / 2) ? (y - N) * delta_k : y * delta_k;  // U_L^-1
//           double k_z = (z > N / 2) ? (z - N) * delta_k : z * delta_k;  // U_L^-1
//
//           double k = sqrt(k_x * k_x + k_y * k_y + k_z * k_z);
//
//           /* Ignore the DC mode */
//           if (k > 0) {
//             // double Tr = gsl_spline2d_eval(rend->spline, k, log_tau, rend->k_acc,
//             //                               rend->tau_acc);
//
//             /* Find the k-space interpolation index */
//             rend_interp_locate_k(rend, k, &k_index, &u_k);
//
//             /* Bilinear interpolation of the ncdm transfer function */
//             double Tr = rend_custom_interp(rend, k_index, tau_index, u_tau, u_k,
//                                          index_f);
//
//             /* The CIC Window function in Fourier space */
//             const double sqrt_W_x = sinc(0.5 * k_x * box_len / N);
//             const double sqrt_W_y = sinc(0.5 * k_y * box_len / N);
//             const double sqrt_W_z = sinc(0.5 * k_z * box_len / N);
//             const double sqrt_W = sqrt_W_x * sqrt_W_y * sqrt_W_z;
//             const double W = sqrt_W * sqrt_W;
//
//             /* Only deconvolve once for the interpolation (no gpart assignment) */
//             Tr /= W;
//
//             fp[half_box_idx(N, x, y, z)][0] *= Tr;
//             fp[half_box_idx(N, x, y, z)][1] *= Tr;
//           } else {
//             fp[half_box_idx(N, x, y, z)][0] = 0;
//             fp[half_box_idx(N, x, y, z)][1] = 0;
//           }
//         }
//       }
//     }
//
//     /* Transform back */
//     fftw_execute(c2r_grid);
//
//     /* Free memory */
//     fftw_destroy_plan(c2r_grid);
//     fftw_destroy_plan(r2c_grid);
//
//     /* Normalization */
//     for (int i = 0; i < N * N * N; i++) {
//       grid[i] /= box_volume;
//     }
//
//     /* Export the data block (only on master node) */
//     if (e->nodeID == 0) {
//       char boxname[40];
//       sprintf(boxname, "grid_%d.hdf5", index_f);
//       writeGRF_H5(grid, N, box_len, boxname);
//     }
//
//   }
//
//   /* Next, compute the potential due to neutrinos (modulo a factor G_newt) */
//
//   /* Calculate the background neutrino density at the present time */
//   const double Omega_nu = cosmology_get_neutrino_density_param(cosmo, cosmo->a);
//   const double Omega_g = cosmo->Omega_g;
//   const double Omega_ur = cosmo->Omega_ur;
//   const double rho_crit0 = cosmo->critical_density_0;
//   /* The comoving density is (Omega_nu * a^-4) * a^3  = Omega_nu / a */
//   const double neutrino_density = Omega_nu * rho_crit0 / cosmo->a;
//   const double photon_density = Omega_g * rho_crit0 / cosmo->a;
//   const double ultra_relativistic_density = Omega_ur * rho_crit0 / cosmo->a;
//
//   /* The starting indices of the respective grids */
//   double *ncdm_grid =
//       rend->density_grid + rend->index_transfer_delta_ncdm * N * N * N;
//   double *g_grid =
//       rend->density_grid + rend->index_transfer_delta_g * N * N * N;
//   double *ur_grid =
//       rend->density_grid + rend->index_transfer_delta_ur * N * N * N;
//   double *HT_prime_grid =
//       rend->density_grid + rend->index_transfer_H_T_Nb_prime * N * N * N;
//   double *HT_prime_prime_grid =
//       rend->density_grid + rend->index_transfer_H_T_Nb_pprime * N * N * N;
//   double *phi_grid = rend->density_grid + rend->index_transfer_phi * N * N * N;
//   double *psi_grid = rend->density_grid + rend->index_transfer_psi * N * N * N;
//
//   /* The potential is multiplied by G_newton later, so for phi, psi & H_T_Nb,
//    * which should not be multiplied by G_newton, we need to divide now.
//    * Also, we multiply by the scale factor a to get peculiar potentials.
//    */
//   const float G_newton = e->physical_constants->const_newton_G;
//   const double potential_factor = cosmo->a / G_newton;
//
//   /* Apply this factor to phi, psi, and the H_T_Nb derivatives */
//   for (int i = 0; i < N * N * N; i++) {
//     phi_grid[i] *= potential_factor;
//     psi_grid[i] *= potential_factor;
//     HT_prime_grid[i] *= potential_factor;
//     HT_prime_prime_grid[i] *= potential_factor;
//   }
//
//   /* Compute RHS of Poisson's equation (modulo G_newton) */
//   for (int i = 0; i < N * N * N; i++) {
//     /* Neutrino contribution */
//     double rho_ncdm = (1.0 + ncdm_grid[i]) * neutrino_density;
//     /* Ultra-relativistic fluid contribution */
//     double rho_ur = (1.0 + ur_grid[i]) * ultra_relativistic_density;
//     /* Photon contribution (gamma) */
//     double rho_g = (1.0 + g_grid[i]) * photon_density;
//     /* H_T_Nb term = (H*a + d/dtau) * (d/dtau) * H_T_Nb */
//     double H_T_term = HT_prime_grid[i] * H_conformal + HT_prime_prime_grid[i];
//
//     /* We will apply the 1/k^2 kernel to the Fourier transform of this */
//     potential[i] = -4 * M_PI * (rho_ncdm + rho_g + rho_ur) + H_T_term;
//
//     /* Note: the phi & psi contribution is added later. */
//   }
//
//   /* Export the grids for troubleshooting */
//   if (e->nodeID == 0) {
//     writeGRF_H5(ncdm_grid, N, box_len, "grid_ncdm.hdf5");
//     writeGRF_H5(g_grid, N, box_len, "grid_g.hdf5");
//     writeGRF_H5(ur_grid, N, box_len, "grid_ur.hdf5");
//     writeGRF_H5(HT_prime_grid, N, box_len, "grid_HT_prime.hdf5");
//     writeGRF_H5(HT_prime_prime_grid, N, box_len, "grid_HT_prime_prime.hdf5");
//     writeGRF_H5(potential, N, box_len, "gr_dens.hdf5");
//     writeGRF_H5(phi_grid, N, box_len, "grid_phi.hdf5");
//     writeGRF_H5(psi_grid, N, box_len, "grid_psi.hdf5");
//   }
//
//   /* Transform to momentum space */
//   fftw_execute(pr2c);
//
//   /* Normalization */
//   for (int i = 0; i < N * N * (N / 2 + 1); i++) {
//     fp[i][0] *= box_volume / (N * N * N);
//     fp[i][1] *= box_volume / (N * N * N);
//   }
//
//   /* Multiply by the 1/k^2 kernel */
//   for (int x = 0; x < N; x++) {
//     for (int y = 0; y < N; y++) {
//       for (int z = 0; z <= N / 2; z++) {
//         double k_x = (x > N / 2) ? (x - N) * delta_k : x * delta_k;  // U_L^-1
//         double k_y = (y > N / 2) ? (y - N) * delta_k : y * delta_k;  // U_L^-1
//         double k_z = (z > N / 2) ? (z - N) * delta_k : z * delta_k;  // U_L^-1
//
//         double k = sqrt(k_x * k_x + k_y * k_y + k_z * k_z);
//
//         double kernel = 1.0 / k / k;
//
//         /* Ignore the DC mode */
//         if (k > 0) {
//           fp[half_box_idx(N, x, y, z)][0] *= kernel;
//           fp[half_box_idx(N, x, y, z)][1] *= kernel;
//         } else {
//           fp[half_box_idx(N, x, y, z)][0] = 0;
//           fp[half_box_idx(N, x, y, z)][1] = 0;
//         }
//       }
//     }
//   }
//
//   /* Transform back */
//   fftw_execute(pc2r);
//
//   /* Normalization */
//   for (int i = 0; i < N * N * N; i++) {
//     potential[i] /= box_volume;
//   }
//
//   /* Export the potentials */
//   if (e->nodeID == 0) {
//     writeGRF_H5(e->mesh->potential, N, box_len, "m_potential.hdf5");
//     writeGRF_H5(potential, N, box_len, "gr_potential_without_stress.hdf5");
//   }
//
//   /* Add the contribution from anisotropic stress = (phi - psi) */
//   for (int i = 0; i < N * N * N; i++) {
//     potential[i] -= (phi_grid[i] - psi_grid[i]);
//   }
//
//   if (e->nodeID == 0) {
//     writeGRF_H5(potential, N, box_len, "gr_potential.hdf5");
//   }
//
//   /* Add the contribution to the gravity mesh */
//   for (int i = 0; i < N * N * N; i++) {
//     e->mesh->potential[i] += potential[i];
//   }
//
//   if (e->nodeID == 0) {
//     writeGRF_H5(e->mesh->potential, N, box_len, "full_potential.hdf5");
//   }
//
//   /* Free memory */
//   fftw_free(potential);
//   fftw_free(fp);
//   fftw_destroy_plan(pr2c);
//   fftw_destroy_plan(pc2r);
//
//
// #else
//   error("No GSL library found. Cannot perform cosmological interpolation.");
// #endif
//
// #else
//   error("No FFTW library found. Cannot compute periodic long-range forces.");
// #endif
// }

/* Depending on the compilation option, add linear theory potentials
 * to the long-range potential mesh. */
void rend_add_to_mesh(struct renderer *rend, const struct engine *e) {

#ifdef RENDERER_NEUTRINOS
  /* Compute the neutrino contribution by applying the linear transfer function
   * to the primordial phases. */
  rend_add_linear_nu_mesh(rend, e);
#endif

#ifdef RENDERER_RESCALED_NEUTRINOS
  /* Compute the neutrino contribution by applying the ratio of linear transfer
   * functions (d_ncdm / d_cdm) to the non-linear cdm phases. */
  rend_add_rescaled_nu_mesh(rend, e);
#endif

#ifdef RENDERER_FULL_GR
  /* Compute all general relativistic potentials (neutrinos, radiation, etc.)
   * by applying the linear transfer function to the primordial phases. */
  rend_add_gr_potential_mesh(rend, e);
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

  int k_size, tau_size, n_functions;

  io_read_attribute(h_grp, "k_size", INT, &k_size);
  io_read_attribute(h_grp, "tau_size", INT, &tau_size);
  io_read_attribute(h_grp, "n_functions", INT, &n_functions);

  tr->k_size = k_size;
  tr->tau_size = tau_size;
  tr->n_functions = n_functions;

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
  tr->titles = (char **)swift_calloc("titles", tr->n_functions, sizeof(char *));

  /* Read the titles of the transfer functions */
  h_attr = H5Aopen(h_grp, "FunctionTitles", H5P_DEFAULT);
  h_tp = H5Aget_type(h_attr);
  h_err = H5Aread(h_attr, h_tp, tr->titles);
  H5Aclose(h_attr);
  H5Tclose(h_tp);

  /* Print the titles */
  for (int i = 0; i < tr->n_functions; i++) {
    message("Loaded perturbation vector '%s'.", tr->titles[i]);
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

  message("We read the perturbation size %d * %d * %d", tr->n_functions,
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
  for (int i = 0; i < tr->k_size; i++) {
    tr->k[i] *= us->UnitLength_in_cgs / file_length_us;
  }

  /* Convert the units of the conformal time. This is log(tau) ! */
  for (int i = 0; i < tr->tau_size; i++) {
    tr->log_tau[i] += log(file_time_us) - log(us->UnitTime_in_cgs);
  }

  // for (int i=0; i<tr->k_size; i++) {
  //     printf("%e\n", tr->k[i]);
  // }
  //
  // for (int i=0; i<tr->tau_size; i++) {
  //     printf("%e\n", tr->log_tau[i]);
  // }
  //
  // for (int i=0; i<tr->k_size * tr->tau_size; i++) {
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

  /* Dimensions of the primordial field */
  const int N = rend->primordial_grid_N;
  const double box_len = rend->primordial_dims[0];

  /* Determine the maximum and minimum wavenumbers we will ever need */
  const double margin = 1.2; //ensure a safe error margin
  const double dk = 2*M_PI/box_len;
  const double lookup_k_max = sqrt(3)*dk*N/2 * margin;
  const double lookup_k_min = dk / margin;

  if (myrank == 0) {

    /* Check if a perturbation file was specified */
    if (strlen(rend->in_perturb_fname) < 1) {
      error("No perturbation file specified.");
    }

    /* Read from disk */
    rend_read_perturb(rend, e, rend->in_perturb_fname);

    /* Initialize our own interpolation spline */
    rend_custom_interp_init(rend, LOOKUP_TABLE_LENGTH, lookup_k_min,
                            lookup_k_max);

// #ifdef RENDERER_FULL_GR
//     rend_grids_alloc(rend);
// #endif
  }

  /* The memory for the transfer functions is located here */
  struct transfer *tr = &rend->transfer;

  /* Broadcast the cosmological perturbations to the other ranks */
#ifdef WITH_MPI

  /* First broadcast the size of the perturbation to the other ranks */
  MPI_Bcast(&tr->k_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&tr->tau_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&tr->n_functions, 1, MPI_INT, 0, MPI_COMM_WORLD);

  /* Allocate memory on the other ranks */
  if (myrank != 0) {
    tr->delta = (double *)swift_malloc(
        "delta", tr->n_functions * tr->k_size * tr->tau_size * sizeof(double));
    tr->k = (double *)swift_malloc("k", tr->k_size * sizeof(double));
    tr->log_tau =
        (double *)swift_malloc("log_tau", tr->tau_size * sizeof(double));
  }

  /* Broadcast the perturbation to the other ranks */
  MPI_Bcast(tr->k, tr->k_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(tr->log_tau, tr->tau_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(tr->delta, tr->k_size * tr->tau_size * tr->n_functions, MPI_DOUBLE,
            0, MPI_COMM_WORLD);

#endif

  /* Identify commonly used indices by their titles */
  if (myrank == 0) {
    /* Initialize the transfer function indices */
    rend->index_transfer_delta_cdm = -1;
    rend->index_transfer_delta_ncdm = -1;
    rend->index_transfer_delta_g = -1;
    rend->index_transfer_delta_ur = -1;
    rend->index_transfer_phi = -1;
    rend->index_transfer_psi = -1;
    rend->index_transfer_H_T_Nb_prime = -1;
    rend->index_transfer_H_T_Nb_pprime = -1;

    /* Find the indices */
    for (int i = 0; i < tr->n_functions; i++) {
      if (strcmp(tr->titles[i], "d_ncdm[0]") == 0) {
        rend->index_transfer_delta_ncdm = i;
        message("Identified ncdm density vector '%s'.", tr->titles[i]);
      } else if (strcmp(tr->titles[i], "d_cdm") == 0) {
        rend->index_transfer_delta_cdm = i;
        message("Identified cdm density vector '%s'.", tr->titles[i]);
      } else if (strcmp(tr->titles[i], "d_g") == 0) {
        rend->index_transfer_delta_g = i;
        message("Identified photon density vector '%s'.", tr->titles[i]);
      } else if (strcmp(tr->titles[i], "d_ur") == 0) {
        rend->index_transfer_delta_ur = i;
        message("Identified ultra-relatistic fluid density vector '%s'.",
                tr->titles[i]);
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
        message("Identified N-body gauge H_T_prime_prime vector '%s'.",
                tr->titles[i]);
      }
    }
  }

#ifdef WITH_MPI
  /* Broadcast the indices to the other ranks */
  MPI_Bcast(&rend->index_transfer_delta_ncdm, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&rend->index_transfer_delta_cdm, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&rend->index_transfer_delta_g, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&rend->index_transfer_delta_ur, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&rend->index_transfer_phi, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&rend->index_transfer_psi, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&rend->index_transfer_H_T_Nb_prime, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&rend->index_transfer_H_T_Nb_pprime, 1, MPI_INT, 0, MPI_COMM_WORLD);

  /* Initialize the interpolation spline on the other ranks */
  if (myrank != 0) {
    // rend_interp_init(rend);
    rend_custom_interp_init(rend, LOOKUP_TABLE_LENGTH, lookup_k_min,
                            lookup_k_max);

// #ifdef RENDERER_FULL_GR
//     rend_grids_alloc(rend);
// #endif
  }
#endif
}

/* Locate the index such that log_tau[index] <= log_tau < log_tau[index+1] and
 *  the fractional distance w = (y - y[index]) / (y[index+1] - y[index]). */
void rend_interp_locate_tau(struct renderer *rend, double log_tau,
                            int *index, double *w) {
  /* The memory for the transfer functions is located here */
  struct transfer *tr = &rend->transfer;

  /* Number of bins */
  int tau_size = tr->tau_size;

  /* Quickly return if we are in the first or last bin */
  if (log_tau < tr->log_tau[0]) {
    *index = 0;
    *w = 0.f;
  } else if (log_tau >= tr->log_tau[tau_size - 1]) {
    *index = tau_size - 1;
    *w = 1.f;
  }

  /* Run through the array in strides of 10 for initial scan */
  for (int i=1; i<tau_size; i++) {
    if (tr->log_tau[i] >= log_tau) {
      *index = i-1;
      break;
    }
  }

  /* Find the bounding values */
  double left = tr->log_tau[*index];
  double right = tr->log_tau[*index + 1];

  /* Calculate the ratio (X - X_left) / (X_right - X_left) */
  *w = (log_tau - left) / (right - left);
}


void rend_custom_interp_init(struct renderer *rend, int table_size,
                             double lookup_k_min, double lookup_k_max) {

    /* Allocate the search table */
    rend->k_acc_table = malloc(table_size * sizeof(double));
    rend->k_acc_table_size = table_size;
    rend->k_acc_prev_index = 0; //updated later

    /* Bounding values for the look-up table */
    rend->lookup_k_min = lookup_k_min;
    rend->lookup_k_max = lookup_k_max;

    /* The perturbation structure containing the vector that is to be searched */
    struct transfer *tr = &rend->transfer;
    const int k_size = tr->k_size;

    /* Make the index table */
    for (int i=0; i<table_size; i++) {
        double u = (double) i/table_size;
        double v = lookup_k_min + u * (lookup_k_max - lookup_k_min);

        /* Find the largest bin such that w > k */
        double maxJ = 0;
        for(int j=0; j<k_size; j++) {
            if (tr->k[j] < v) {
                maxJ = j;
            }
        }
        rend->k_acc_table[i] = maxJ;
    }
}

/* Locate the index such that k[index] <= k < k[index+1] and the fractional
 * distance w = (k - k[index]) / (k[index+1] - k[index]).
 * To limit the number of branches, we do not check if the value is in bounds!
 */
void rend_interp_locate_k(struct renderer *rend, double k,
                          int *index, double *w) {

  /* The transfer functions structure */
  struct transfer *tr = &rend->transfer;

  /* Before doing anything else, check if the last search is still valid */
  int last_index = rend->k_acc_prev_index;
  if (k >= tr->k[last_index] && k < tr->k[last_index + 1]) {
    /* We found the index */
    *index = last_index;
  } else {
    /* We will use the fast look-up table */
    const int k_acc_table_size = rend->k_acc_table_size;
    const int k_size = tr->k_size;

    /* Bounding values in the lookup table */
    const double lookup_k_min = rend->lookup_k_min;
    const double lookup_k_max = rend->lookup_k_max;

    /* Quickly find a starting index using the look-up table */
    const double u = (k - lookup_k_min) / (lookup_k_max - lookup_k_min);
    const int I = floor(u * k_acc_table_size);
    const int start = rend->k_acc_table[I];

    /* Search in the k vector, starting from the looked up index */
    int i;
    for (i = start; i < k_size; i++) {
      if (k >= tr->k[i] && k <= tr->k[i + 1]) break;
    }

    /* We found the index */
    *index = i;

    /* Store for later searches */
    rend->k_acc_prev_index = i;
  }

  /* Find the bounding values */
  const double left = tr->k[*index];
  const double right = tr->k[*index + 1];

  /* Calculate the ratio (X - X_left) / (X_right - X_left) */
  *w = (k - left) / (right - left);
}

/* Bilinear interpolation */
double rend_custom_interp(struct renderer *rend, int k_index, int tau_index,
                          double u_tau, double u_k, int index_src) {

  /* Bounding values for the larger table */
  struct transfer *tr = &rend->transfer;
  int k_size = tr->k_size;
  int tau_size = tr->tau_size;

  /* Select the desired transfer function */
  double *arr = tr->delta + index_src * k_size * tau_size;

  /* Retrieve the bounding values */
  double T11 = arr[k_size * tau_index + k_index];
  double T21 = arr[k_size * tau_index + k_index + 1];
  double T12 = arr[k_size * (tau_index + 1) + k_index];
  double T22 = arr[k_size * (tau_index + 1) + k_index + 1];

  return (1 - u_tau) * ((1 - u_k) * T11 + u_k * T21)
             + u_tau * ((1 - u_k) * T12 + u_k * T22);
}
