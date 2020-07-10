/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Willem Elbers (whe@willemelbers.com).
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

/**
 * renderer.h  -  linear response implementation of the relic neutrinos
 */

#ifndef SWIFT_RENDERER_H
#define SWIFT_RENDERER_H

#include "../common_io.h"
#include "../engine.h"

/* We use CLASS for the transfer functions */
#ifdef WITH_CLASS_INTERFACE
#include "class.h"
#endif

/* We use GSL for accelerated 2D interpolation */
// #ifdef HAVE_LIBGSL
// #include <gsl/gsl_spline2d.h>
// #endif

/* We use FFTW for Fourier transforming the primordial Gaussian field */
#ifdef HAVE_FFTW
#include <fftw3.h>
#endif

/* We use HDF5 to load the primordial Gaussian field */
#ifdef HAVE_HDF5
#include <hdf5.h>
#endif

#define LOOKUP_TABLE_LENGTH 1000

/**
 * @brief Data structure for the renderer, generating perturbation theory grids
 */
struct renderer {

  /*! Array containing the primordial fluctuations on a grid */
  double *primordial_grid;
  double *primordial_dims;
  int primordial_grid_N;

  /*! Array containing sequency of perturbation theory grids */
  double *the_grids;

  /*! Pointer to the neutrino over-density field, corresponding to the first
      N * N * N doubles in the_grids */
  double *density_grid;

  /* File name of the perturbation data */
  char *in_perturb_fname;
  char *out_perturb_fname;

  /* Names of the CLASS parameter & precision files */
  char *class_ini_fname;
  char *class_pre_fname;

  /* Commonly used indices of transfer functions */
  int index_transfer_delta_cdm;
  int index_transfer_delta_ncdm;
  int index_transfer_delta_g;
  int index_transfer_delta_ur;
  int index_transfer_phi;
  int index_transfer_psi;
  int index_transfer_H_T_Nb_prime;
  int index_transfer_H_T_Nb_pprime;

  /*! Neutrino transfer functions */
  struct transfer {
    int k_size;
    int tau_size;
    int n_functions;
    double *delta;
    double *k;
    double *log_tau;
    char **titles;
  } transfer;

  /*! Desired length of the neutrino perturbation along the k dimension */
  int num_of_k_bins;  // user-defined

  /* Search table for interpolation acceleration in the k direction */
  int k_acc_table_size;
  double *k_acc_table;
  int k_acc_last_index;

// #ifdef HAVE_LIBGSL
//   /* GSL interpolation objects */
//   const gsl_interp2d_type *interp_type;
//   gsl_interp_accel *k_acc;
//   gsl_interp_accel *tau_acc;
//   gsl_spline2d *spline;
// #endif
};

/* The renderer object renders transfer functions onto the grid */
void rend_init(struct renderer *rend, struct swift_params *params,
               const struct engine *e);
void rend_clean(struct renderer *rend);

/* Input of the primordial Gaussian field data */
void rend_load_primordial_field(struct renderer *rend, const char *fname);

/* Rendering the transfer functions onto the primordial field */
void rend_add_to_mesh(struct renderer *rend, const struct engine *e);
void rend_add_rescaled_nu_mesh(struct renderer *rend, const struct engine *e);
void rend_add_linear_nu_mesh(struct renderer *rend, const struct engine *e);
void rend_add_gr_potential_mesh(struct renderer *rend, const struct engine *e);

/* The GSL interpolation structures */
// void rend_interp_init(struct renderer *rend);
// void rend_interp_switch_source(struct renderer *rend, int index_src);
// void rend_interp_free(struct renderer *rend);
void rend_grids_alloc(struct renderer *rend);

/* Custom interpolation functions */
void rend_custom_interp_init(struct renderer *rend, int table_size);
void rend_interp_locate_tau(struct renderer *rend, double log_tau,
                            int *index, double *w);
void rend_interp_locate_k(struct renderer *rend, double k, int *index,
                          double *w);
double rend_custom_interp(struct renderer *rend, int k_index, int tau_index,
                          double u_tau, double u_k, int index_src);

/* Input and output for the perturbation data */
void rend_read_perturb(struct renderer *rend, const struct engine *e,
                       char *fname);
void rend_write_perturb(struct renderer *rend, const struct engine *e,
                        char *fname);

void rend_init_perturb_vec(struct renderer *rend, struct swift_params *params,
                           const struct engine *e, int myrank);

/* Additional functions related to CLASS */
#ifdef WITH_CLASS_INTERFACE
#include "class_interface.h"
#endif

#endif /* SWIFT_RENDERER_H */
