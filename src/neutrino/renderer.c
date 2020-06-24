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

/* We use Firebolt for linear theory neutrino calculations */
#ifdef WITH_FIREBOLT_INTERFACE
#include "firebolt_interface.h"
#endif

/* Array index (this is the row major format) */
inline int box_idx(int N, int x, int y, int z) { return z + N * (y + N * x); }

/* Row major index for a half-complex array with N*N*(N/2+1) complex entries */
inline int half_box_idx(int N, int x, int y, int z) {
  return z + (N / 2 + 1) * (y + N * x);
}

/* Sinc function */
inline double sinc(double x) { return x == 0 ? 0. : sin(x) / x; }

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

/* Read binary boxes in HDF5 format */
int readGRF_H5(double **box, int *N, double *box_len, const char *fname) {
  /* Open the hdf5 file */
  hid_t h_file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  /* Open the Header group */
  hid_t h_grp = H5Gopen(h_file, "Header", H5P_DEFAULT);

  /* Read the size of the field */
  hid_t h_attr, h_err;
  double boxsize[3];

  /* Open and read out the attribute */
  h_attr = H5Aopen(h_grp, "BoxSize", H5P_DEFAULT);
  h_err = H5Aread(h_attr, H5T_NATIVE_DOUBLE, &boxsize);
  assert(h_err >= 0);

  /* It should be a cube */
  assert(boxsize[0] == boxsize[1]);
  assert(boxsize[1] == boxsize[2]);
  *box_len = boxsize[0];

  /* Close the attribute, and the Header group */
  H5Aclose(h_attr);
  H5Gclose(h_grp);

  /* Open the Field group */
  h_grp = H5Gopen(h_file, "Field", H5P_DEFAULT);

  /* Open the Field dataset */
  hid_t h_data = H5Dopen2(h_grp, "GaussianRandomField", H5P_DEFAULT);

  /* Open the dataspace and fetch the grid dimensions */
  hid_t h_space = H5Dget_space(h_data);
  int ndims = H5Sget_simple_extent_ndims(h_space);
  hsize_t *dims = malloc(ndims * sizeof(hsize_t));
  H5Sget_simple_extent_dims(h_space, dims, NULL);

  /* We should be in 3D */
  assert(ndims == 3);
  /* It should be a cube */
  assert(dims[0] == dims[1]);
  assert(dims[1] == dims[2]);
  *N = dims[0];

  /* Allocate the array */
  *box = swift_malloc("box", dims[0] * dims[1] * dims[2] * sizeof(double));

  /* Read out the data */
  h_err = H5Dread(h_data, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, *box);

  /* Close the dataspace and dataset */
  H5Sclose(h_space);
  H5Dclose(h_data);
  free(dims);

  /* Close the Field group */
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
  rend->class_ini_fname = (char *)malloc(200 * sizeof(char));
  rend->class_pre_fname = (char *)malloc(200 * sizeof(char));
  parser_get_opt_param_string(params, "Boltzmann:class_ini_file",
                              rend->class_ini_fname, "");
  parser_get_opt_param_string(params, "Boltzmann:class_pre_file",
                              rend->class_pre_fname, "");

  /* Open and load the file with the primordial Gaussian field */
  double *primordial_field;
  double box_len, box_volume;
  int N;
  readGRF_H5(&primordial_field, &N, &box_len, fieldFName);

  /* Store the grid dimensions */
  rend->primordial_box_N = N;
  rend->primordial_box_len = box_len;
  box_volume = box_len * box_len * box_len;

  /* Print the loaded field dimensions */
  message("Loaded %d^3 primordial grid with physical side length: %.1f U_L "
          "on this node.", N, box_len);

  /* Verify that the dimensions of the primordial field match the gravity mesh */
  for (int i = 0; i < 3; i++) {
    if (fabs(box_len - e->s->dim[i]) / e->s->dim[i] > 1e-3) {
      error("Dimensions[%i] of primordial field do not agree with space.", i);
    }
  }

  /* Verify that the primordial grid has the same size as the gravity mesh */
  if (N != e->mesh->N) {
    error("Primordial grid is not the same size as the gravity mesh %d != %d.",
          N, e->mesh->N);
  }

  /* Allocate memory for the complex phases */
  rend->primordial_phases = (fftw_complex *)swift_malloc("primordial_phases",
                                    sizeof(fftw_complex) * N * N * (N / 2 + 1));

  /* Compute the Fourier transform of the primordial field to get the phases */
  fftw_plan r2c = fftw_plan_dft_r2c_3d(N, N, N, primordial_field,
                                       rend->primordial_phases, FFTW_ESTIMATE);
  fftw_execute(r2c);

  /* Normalization */
  for (int i = 0; i < N * N * (N / 2 + 1); i++) {
    rend->primordial_phases[i][0] *= box_volume / (N * N * N);
    rend->primordial_phases[i][1] *= box_volume / (N * N * N);
  }

  /* Destroy the plan */
  fftw_destroy_plan(r2c);

  /* The size of the smaller (down-sampled) Gaussian random field grid */
  int M = parser_get_opt_param_int(params,
      "Boltzmann:small_mesh_side_length", rend->primordial_box_N);
  rend->primordial_box_small_N = M;

  /* If needed, compute Fourier transform of a smaller copy of the grid */
  if (M != N) {

    /* Down-sampling factor */
    int f = N / M;
    message("Computing %d^3 downsampled version of the primordial grid.", M);

    /* Make sure that the smaller grid divides into the larger grid */
    if (f * M != N) {
      error("Small grid does not divide: [small, large] = [%d, %d].", N, M);
    }

    /* Allocate memory for the down-sampled grid */
    double *small = swift_calloc("small", M * M * M, sizeof(double));

    /* Compute the down-sampled version of the primordial grid by averaging */
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        for (int k = 0; k < N; k++) {
          int small_index = box_idx(M, i / f, j / f, k / f);
          int large_index = box_idx(N, i, j, k);
          small[small_index] += primordial_field[large_index] / (f * f * f);
        }
      }
    }

    /* Allocate memory for the complex phases of the down-sampled box */
    rend->primordial_phases_small = (fftw_complex *)swift_malloc("small_phases",
                                    sizeof(fftw_complex) * M * M * (M / 2 + 1));

    /* Compute the Fourier transform of the primordial field to get the phases */
    fftw_plan r2c_small = fftw_plan_dft_r2c_3d(M, M, M, small,
                                  rend->primordial_phases_small, FFTW_ESTIMATE);
    fftw_execute(r2c_small);

    /* Normalization */
    for (int i = 0; i < M * M * (M / 2 + 1); i++) {
      rend->primordial_phases_small[i][0] *= box_volume / (M * M * M);
      rend->primordial_phases_small[i][1] *= box_volume / (M * M * M);
    }

    /* Destroy the plan */
    fftw_destroy_plan(r2c_small);
    message("Computed Fourier transform of small primordial grid.");

    /* We are done with the down-sampled version */
    swift_free("small", small);
  } else {
    /* No size difference, so let the small grid point to the same array */
    rend->primordial_phases_small = rend->primordial_phases;
  }

  /* We are done with the configuration space primordial field. */
  free(primordial_field);

#ifdef WITH_FIREBOLT_INTERFACE
  /* Initialize the Firebolt interface */
  firebolt_init(params, rend, e);
#endif
}

void rend_grids_alloc(struct renderer *rend) {
  /* Allocate memory for the perturbation theory grids */
  const int N = rend->primordial_box_N;
  const int Nf = rend->transfer.n_functions;
  const int bytes = N * N * N * Nf * sizeof(double);
  rend->the_grids = (double *)swift_malloc("the_grids", bytes);

  if (rend->the_grids == NULL) {
    error("Error allocating memory for perturbation theory grids.");
  }

}

void rend_clean(struct renderer *rend) {
  /* Free the primordial phases */
  swift_free("primordial_phases", rend->primordial_phases);

  /* Free the down-sampled version if it exists */
  if (rend->primordial_box_N != rend->primordial_box_small_N) {
    swift_free("small_phases", rend->primordial_phases_small);
  }

  /* Free strings */
  free(rend->in_perturb_fname);
  free(rend->out_perturb_fname);
  free(rend->class_ini_fname);
  free(rend->class_pre_fname);

  /* Free interpolation tables */
  swift_free("delta", rend->transfer.delta);
  swift_free("k", rend->transfer.k);
  swift_free("log_tau", rend->transfer.log_tau);

  /* Free function title strings */
  for (int i = 0; i < rend->transfer.n_functions; i++) {
    free(rend->transfer.titles[i]);
  }
  swift_free("titles", rend->transfer.titles);

  /* Free density & perturbation theory grids */
#ifdef RENDERER_FULL_GR
  free(rend->the_grids);
#endif

  /* Free the index table */
  swift_free("k_acc_table", rend->k_acc_table);

#ifdef WITH_FIREBOLT_INTERFACE
  /* Clean up the Firebolt interface */
  firebolt_free();
#endif
}

/* Add neutrinos using the Bird & Ali-Haïmoud method */
void rend_add_rescaled_nu_mesh(struct renderer *rend, const struct engine *e) {
#ifdef HAVE_FFTW
#ifdef HAVE_LIBGSL
  /* Grid size */
  const int N = rend->primordial_box_N;
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

  /* Apply the ncdm transfer function, undo the cdm transfer function, and undo
   * the long-range kernel. */
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

          /* Bilinear interpolation of the ncdm transfer function delta_ncdm */
          double d_ncdm = rend_custom_interp(rend, k_index, tau_index, u_tau, u_k,
                                             rend->index_transfer_delta_ncdm);
          /* Bilinear interpolation of the cdm transfer function delta_cdm */
          double d_cdm = rend_custom_interp(rend, k_index, tau_index, u_tau, u_k,
                                            rend->index_transfer_delta_cdm);

          /* The long-range kernel */
          double K = 1.;
          fourier_kernel_long_grav_eval(k * k * r_s * r_s, &K);

          /* The overall factor to be applied to the cdm long-range potential */
          double factor = (d_ncdm / d_cdm) * bg_density_ratio / K;

          fp[half_box_idx(N, x, y, z)][0] *= factor;
          fp[half_box_idx(N, x, y, z)][1] *= factor;
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
    writeGRF_H5(e->mesh->potential, N, box_len, "m_potential.hdf5");
    writeGRF_H5(potential, N, box_len, "scaled_nu_potential.hdf5");
  }

  double Q = 0.f;
  double R = 0.f;

  /* Add the contribution to the gravity mesh */
  for (int i = 0; i < N * N * N; i++) {
    R += e->mesh->potential[i] * e->mesh->potential[i];
    Q += potential[i] * potential[i];
    e->mesh->potential[i] += potential[i];
  }

  message("[Q, R] = [%e, %e]", sqrt(Q/(N*N*N)), sqrt(R/(N*N*N)));

  writeGRF_H5(e->mesh->potential, N, box_len, "full_potential.hdf5");

  /* Free memory */
  fftw_free(potential);
  fftw_free(fp);
  fftw_destroy_plan(pr2c);
  fftw_destroy_plan(pc2r);

#else
  error("No GSL library found. Cannot perform cosmological interpolation.");
#endif

#else
  error("No FFTW library found. Cannot compute periodic long-range forces.");
#endif
}

/* Add neutrinos using the linear theory transfer function */
void rend_add_linear_nu_mesh(struct renderer *rend, const struct engine *e) {
#ifdef HAVE_FFTW
#ifdef HAVE_LIBGSL
  /* Grid size */
  const int N = rend->primordial_box_N;
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

  /* Find the time index, used for interpolation later on */
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

  /* Copy the Fourier transform of the primordial field into the complex array */
  fftw_complex *source_address = rend->primordial_phases;
  memcpy(fp, source_address, sizeof(fftw_complex) * N * N * (N / 2 + 1));

  /* Prepare the FFTW plans */
  fftw_plan pr2c = fftw_plan_dft_r2c_3d(N, N, N, potential, fp, FFTW_ESTIMATE);
  fftw_plan pc2r = fftw_plan_dft_c2r_3d(N, N, N, fp, potential, FFTW_ESTIMATE);

  /* Apply the neutrino transfer function & inverse Poisson kernel */
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

          /* Bilinear interpolation of the ncdm transfer function delta_ncdm */
          double d_ncdm = rend_custom_interp(rend, k_index, tau_index, u_tau, u_k,
                                            rend->index_transfer_delta_ncdm);

          /* The CIC Window function in Fourier space */
          double W_x = (k_x == 0) ? 1 : pow(sinc(0.5 * k_x * box_len / N), 2);
          double W_y = (k_y == 0) ? 1 : pow(sinc(0.5 * k_y * box_len / N), 2);
          double W_z = (k_z == 0) ? 1 : pow(sinc(0.5 * k_z * box_len / N), 2);
          double W = W_x * W_y * W_z;
          double WW_CIC = W * W;

          /* The inverse Poisson kernel */
          double kernel = -4 * M_PI / k / k;

          /* The overall factor to be applied to the primordial phases */
          double factor = d_ncdm * neutrino_density * kernel / WW_CIC;

          fp[half_box_idx(N, x, y, z)][0] *= factor;
          fp[half_box_idx(N, x, y, z)][1] *= factor;
        } else {
          fp[half_box_idx(N, x, y, z)][0] = 0;
          fp[half_box_idx(N, x, y, z)][1] = 0;
        }
      }
    }
  }

  /* Transform from complex to real */
  fftw_execute(pc2r);

  /* Normalization */
  for (int i = 0; i < N * N * N; i++) {
    potential[i] /= box_volume;
  }

  /* Export the potentials */
  if (e->nodeID == 0) {
    // writeGRF_H5(e->mesh->potential, N, box_len, "m_potential.hdf5");
    // writeGRF_H5(potential, N, box_len, "linear_nu_potential.hdf5");

    /* Store the box as a separate box for the timestep */
    if (e->step % 10 == 0) {
      char one[40];
      char two[40];
      double z = e->cosmology->z;
      sprintf(one, "m_potential_z_%.2f.hdf5", z);
      sprintf(two, "linear_nu_potential_z_%.2f.hdf5", z);
      writeGRF_H5(e->mesh->potential, N, box_len, one);
      writeGRF_H5(potential, N, box_len, two);
    }
  }

  double Q = 0.f;
  double R = 0.f;

  /* Add the contribution to the gravity mesh */
  for (int i = 0; i < N * N * N; i++) {
    R += e->mesh->potential[i] * e->mesh->potential[i];
    Q += potential[i] * potential[i];
    e->mesh->potential[i] += potential[i];
  }

  message("[Q, R] = [%e, %e]", sqrt(Q/(N*N*N)), sqrt(R/(N*N*N)));

  if (e->nodeID == 0) {
    writeGRF_H5(e->mesh->potential, N, box_len, "full_potential.hdf5");
  }

  /* Free memory */
  fftw_free(potential);
  fftw_free(fp);
  fftw_destroy_plan(pr2c);
  fftw_destroy_plan(pc2r);

#else
  error("No GSL library found. Cannot perform cosmological interpolation.");
#endif

#else
  error("No FFTW library found. Cannot compute periodic long-range forces.");
#endif
}

void rend_add_gr_potential_mesh(struct renderer *rend, const struct engine *e) {
#ifdef HAVE_FFTW
#ifdef HAVE_LIBGSL
  /* Grid size */
  const int N = rend->primordial_box_N;
  const double box_len = e->s->dim[0];
  const double box_volume = pow(box_len, 3);
  const double delta_k = 2 * M_PI / box_len;  // U_L^-1

  /* Current conformal time */
  const struct cosmology *cosmo = e->cosmology;
  const double tau = cosmo->conformal_time;
  const double H_conformal = cosmo->H * cosmo->a;

  // message("H-conformal = %f", H_conformal);

  /* Prevent out of interpolation range error */
  const int tau_size = rend->transfer.tau_size;
  const double final_log_tau = rend->transfer.log_tau[tau_size - 1];
  const double log_tau = min(log(tau), final_log_tau);

  /* Bilinear interpolation indices in (log_tau, k) space */
  int tau_index = 0, k_index = 0;
  double u_tau = 0.f, u_k = 0.f;

  /* Find the time index */
  rend_interp_locate_tau(rend, log_tau, &tau_index, &u_tau);

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
  for (int index_f = 0; index_f < rend->transfer.n_functions; index_f++) {
    /* Switch the interpolation spline to the desired transfer function */
    // rend_interp_switch_source(rend, index_f);

    /* Use memory that has already been allocated */
    double *grid = rend->the_grids + index_f * N * N * N;

    /* First, copy the primordial phases into the complex array */
    fftw_complex *source_address = rend->primordial_phases;
    fftw_complex *destination = fp;
    memcpy(destination, source_address, sizeof(fftw_complex) * N * N * (N / 2 + 1));

    /* Create plans */
    fftw_plan r2c_grid = fftw_plan_dft_r2c_3d(N, N, N, grid, fp, FFTW_ESTIMATE);
    fftw_plan c2r_grid = fftw_plan_dft_c2r_3d(N, N, N, fp, grid, FFTW_ESTIMATE);

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
            // double Tr = gsl_spline2d_eval(rend->spline, k, log_tau, rend->k_acc,
            //                               rend->tau_acc);

            /* Find the k-space interpolation index */
            rend_interp_locate_k(rend, k, &k_index, &u_k);

            /* Bilinear interpolation of the ncdm transfer function */
            double Tr = rend_custom_interp(rend, k_index, tau_index, u_tau, u_k,
                                         index_f);

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

    /* Free memory */
    fftw_destroy_plan(c2r_grid);
    fftw_destroy_plan(r2c_grid);

    /* Normalization */
    for (int i = 0; i < N * N * N; i++) {
      grid[i] /= box_volume;
    }

    /* Export the data block (only on master node) */
    if (e->nodeID == 0) {
      char boxname[40];
      sprintf(boxname, "grid_%d.hdf5", index_f);
      writeGRF_H5(grid, N, box_len, boxname);
    }

  }

  /* Next, compute the potential due to neutrinos (modulo a factor G_newt) */

  /* Calculate the background neutrino density at the present time */
  const double Omega_nu = cosmology_get_neutrino_density_param(cosmo, cosmo->a);
  const double Omega_g = cosmo->Omega_g;
  const double Omega_ur = cosmo->Omega_ur;
  const double rho_crit0 = cosmo->critical_density_0;
  /* The comoving density is (Omega_nu * a^-4) * a^3  = Omega_nu / a */
  const double neutrino_density = Omega_nu * rho_crit0 / cosmo->a;
  const double photon_density = Omega_g * rho_crit0 / cosmo->a;
  const double ultra_relativistic_density = Omega_ur * rho_crit0 / cosmo->a;

  /* The starting indices of the respective grids */
  double *ncdm_grid =
      rend->the_grids + rend->index_transfer_delta_ncdm * N * N * N;
  double *g_grid =
      rend->the_grids + rend->index_transfer_delta_g * N * N * N;
  double *ur_grid =
      rend->the_grids + rend->index_transfer_delta_ur * N * N * N;
  double *HT_prime_grid =
      rend->the_grids + rend->index_transfer_H_T_Nb_prime * N * N * N;
  double *HT_prime_prime_grid =
      rend->the_grids + rend->index_transfer_H_T_Nb_pprime * N * N * N;
  double *phi_grid = rend->the_grids + rend->index_transfer_phi * N * N * N;
  double *psi_grid = rend->the_grids + rend->index_transfer_psi * N * N * N;

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
    writeGRF_H5(ncdm_grid, N, box_len, "grid_ncdm.hdf5");
    writeGRF_H5(g_grid, N, box_len, "grid_g.hdf5");
    writeGRF_H5(ur_grid, N, box_len, "grid_ur.hdf5");
    writeGRF_H5(HT_prime_grid, N, box_len, "grid_HT_prime.hdf5");
    writeGRF_H5(HT_prime_prime_grid, N, box_len, "grid_HT_prime_prime.hdf5");
    writeGRF_H5(potential, N, box_len, "gr_dens.hdf5");
    writeGRF_H5(phi_grid, N, box_len, "grid_phi.hdf5");
    writeGRF_H5(psi_grid, N, box_len, "grid_psi.hdf5");
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
    writeGRF_H5(e->mesh->potential, N, box_len, "m_potential.hdf5");
    writeGRF_H5(potential, N, box_len, "gr_potential_without_stress.hdf5");
  }

  /* Add the contribution from anisotropic stress = (phi - psi) */
  for (int i = 0; i < N * N * N; i++) {
    potential[i] -= (phi_grid[i] - psi_grid[i]);
  }

  if (e->nodeID == 0) {
    writeGRF_H5(potential, N, box_len, "gr_potential.hdf5");
  }

  /* Add the contribution to the gravity mesh */
  for (int i = 0; i < N * N * N; i++) {
    e->mesh->potential[i] += potential[i];
  }

  if (e->nodeID == 0) {
    writeGRF_H5(e->mesh->potential, N, box_len, "full_potential.hdf5");
  }

  /* Free memory */
  fftw_free(potential);
  fftw_free(fp);
  fftw_destroy_plan(pr2c);
  fftw_destroy_plan(pc2r);

#else
  error("No GSL library found. Cannot perform cosmological interpolation.");
#endif

#else
  error("No FFTW library found. Cannot compute periodic long-range forces.");
#endif
}

/* Depending on the compilation option, add linear theory potentials
 * to the long-range potential mesh. */
void rend_add_to_mesh(struct renderer *rend, const struct engine *e) {

#ifdef WITH_FIREBOLT_INTERFACE
  /* Make sure that Firebolt has updated all the multipoles */
  firebolt_update(rend, e);
#endif


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

  if (myrank == 0) {

    /* If a perturbation file & a CLASS ini file are both specified */
    if (strlen(rend->in_perturb_fname) > 1 &&
        strlen(rend->class_ini_fname) > 1) {
      error(
          "Specified both perturbation file & CLASS .ini file. '%s' '%s' (%ld)",
          rend->in_perturb_fname, rend->class_ini_fname,
          strlen(rend->class_ini_fname));
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
    // rend_interp_init(rend);
    rend_custom_interp_init(rend, 100);

#ifdef RENDERER_FULL_GR
    rend_grids_alloc(rend);
#endif
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
    rend_custom_interp_init(rend, 20);

#ifdef RENDERER_FULL_GR
    rend_grids_alloc(rend);
#endif
  }
#endif
}

/* Locate the greatest lower bounding index in the log_tau table */
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


void rend_custom_interp_init(struct renderer *rend, int table_size) {
    /* Allocate the search table */
    rend->k_acc_table = swift_malloc("k_acc_table", table_size * sizeof(double));
    rend->k_acc_table_size = table_size;

    /* Bounding values for the larger table */
    struct transfer *tr = &rend->transfer;
    int k_size = tr->k_size;
    double k_min = tr->k[0];
    double k_max = tr->k[k_size-1];

    /* Make the index table */
    for (int i=0; i<table_size; i++) {
        double u = (double) i/table_size;
        double v = k_min + u * (k_max - k_min);

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

/* Locate the greatest lower bounding index in the log_tau table */
void rend_interp_locate_k(struct renderer *rend, double k,
                          int *index, double *w) {

  /* Bounding values for the larger table */
  struct transfer *tr = &rend->transfer;
  int k_acc_table_size = rend->k_acc_table_size;
  int k_size = tr->k_size;
  double k_min = tr->k[0];
  double k_max = tr->k[k_size-1];

  if (k > k_max) {
      *index = k_size - 1;
      *w = 1.0;
      return;
  }

  /* Quickly find a starting index using the indexed seach */
  double v = log(k);
  double u = (v - k_min) / (k_max - k_min);
  int J = floor(u * k_acc_table_size);
  int idx = rend->k_acc_table[J < k_acc_table_size ? J : k_acc_table_size - 1];

  /* Search in the k vector */
  int i;
  for (i = idx; i < k_size; i++) {
    if (k >= tr->k[i] && k <= tr->k[i + 1]) break;
  }

  /* We found the index */
  *index = i;

  /* Find the bounding values */
  double left = tr->k[*index];
  double right = tr->k[*index + 1];

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
