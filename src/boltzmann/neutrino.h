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
 * neutrino.h  -  linear response implementation of the relic neutrinos
 */

 #ifndef SWIFT_NEUTRINO_H
 #define SWIFT_NEUTRINO_H

 #define BOLTZ_DEFAULT_BINS 20

 #include "../engine.h"
 #include "../mesh_gravity.h"
 #include "../common_io.h"
 #include "powerspec.h"

 #ifdef HAVE_HDF5
 #include <hdf5.h>
 #endif

/**
 * @brief Data structure for the Boltzmann solver
 */
struct boltz {

    /*! The gravity mesh */
    struct pm_mesh *mesh;

    /*! Array containing the primordial fluctuations on a grid */
    float* primordial_grid;

    /*! Array containing the neutrino perturbation \Psi */
    float* Psi;
    /*! Array containing the time derivative of \Psi */
    float* Psi_dot;

    /*! Array containing the momentum bins */
    float* q_bins;
    /*! Array containing the wavenumber bins */
    float* k_bins;

    /*! Dimensions of the neutrino perturbation */
    size_t Nl; //number of multipoles
    size_t Nk; //number of wavenumbers
    size_t Nq; //number of momentum bins
    size_t Nn; //number of species (not used yet)

    /* Here be the power spectrum measurements */
    size_t bins; //user-defined
  	double *k_in_bins;
  	double *power_in_bins;
  	int *obs_in_bins;
    bool refactored;

    /* Logirithmic derivative of the 0th order distribution function */
    double* dlogf0_dlogq;
    double f0_prefactor; //this is g_s/h^3

    /* Here be Boltzmann related quantities */
    double* phi;
    double* phi_dot;

    /* Neutrino related constants (should be user-defined) */
    double m_nu; //neutrino mass (just one for now)
};

/**
 * @brief Inline 0th order neutrino distribution function without g/h^3 prefactor
 *
 * @param q Energy or momentum divided by k_b*T_nu
 */
__attribute__((always_inline)) INLINE static double f0_nu(double q) {
    return 1./(exp(q)+1.);
}

void boltz_init(struct boltz *bolt, struct swift_params *params, struct engine *e);
void boltz_export_phi(struct boltz *bolt, const char *fname);
void boltz_update_phi(struct boltz *bolt);
void boltz_refactor_bins(struct boltz *bolt, size_t bins);
void boltz_step(struct boltz *bolt, struct engine *e);

#endif /* SWIFT_NEUTRINO_H */
