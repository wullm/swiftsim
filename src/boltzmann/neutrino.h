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

 #include "../mesh_gravity.h"
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
};

void boltz_init(struct boltz *bolt, struct pm_mesh *mesh);



#endif /* SWIFT_NEUTRINO_H */
