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

#define BOLTZ_DEFAULT_BINS 20

#include "../common_io.h"
#include "../engine.h"

/* We use CLASS for the transfer functions */
#ifdef WITH_CLASS_INTERFACE
#include "class.h"
#endif

/* We use FFTW for Fourier transforming the primordial Gaussian field */
#ifdef HAVE_FFTW
#include <fftw3.h>
#endif

/* We use HDF5 to load the primordial Gaussian field */
#ifdef HAVE_HDF5
#include <hdf5.h>
#endif

/**
 * @brief Data structure for the renderer, generating perturbation theory grids
 */
struct renderer {

  /*! Array containing the primordial fluctuations on a grid */
  double *primordial_grid;
  double *primordial_dims;
  size_t primordial_grid_N;

  /* File name of the perturbation data */
  char *in_perturb_fname;
  char *out_perturb_fname;

  /* Names of the CLASS parameter & precision files */
  char *class_ini_fname;
  char *class_pre_fname;

  /*! Neutrino transfer functions */
  struct transfer {
    size_t k_size;
    size_t tau_size;
    double *delta;
    double *k;
    double *log_tau;
  } transfer;

  /*! Desired length of the neutrino perturbation along the k dimension */
  size_t num_of_k_bins;  // user-defined
};

/* The renderer object renders transfer functions onto the grid */
void rend_init(struct renderer *rend, struct swift_params *params,
               const struct engine *e);
void rend_clean(struct renderer *rend);

/* Input of the primordial Gaussian field data */
void rend_load_primordial_field(struct renderer *rend, const char *fname);

/* Rendering the transfer functions onto the primordial field */
void rend_add_to_mesh(struct renderer *rend, const struct engine *e);

/* The GSL interpolation structures */
void rend_interp_init(struct renderer *rend);
void rend_interp_free(struct renderer *rend);

/* Input and output for the perturbation data */
void rend_read_perturb(struct renderer *rend, const struct engine *e,
                       char *fname);
void rend_write_perturb(struct renderer *rend, const struct engine *e,
                        char *fname);

/* Additional functions related to CLASS */
#ifdef WITH_CLASS_INTERFACE
#include "class_interface.h"
#endif

#endif /* SWIFT_RENDERER_H */
