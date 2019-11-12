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

#ifdef HAVE_FFTW
#include <fftw3.h>
#endif

/* This object's header. */
#include "neutrino.h"

void boltz_init(struct boltz *bolt, struct pm_mesh *mesh) {
    //Initialize mesh reference
    bolt->mesh = mesh;

    //Open the file containing the neutrino perturbation (e.g. from CLASS)
    const hid_t h_file = H5Fopen("multipoles_0.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
    if (h_file < 0) error("Error opening the multipoles file.");

    //First, load \Psi(k,q,l) from the main dataset
    hid_t h_data = H5Dopen(h_file, "Psi(k,q,l)", H5P_DEFAULT);
    hid_t space = H5Dget_space(h_data);

    //The number of dimensions in the dataset (expected 3)
    const int rank = H5Sget_simple_extent_ndims(space);
    if (rank != 3) {
        error("Incorrect dataset dimensions for the neutrino perturbation.");
    }

    //Find the extent of each dimension
    hsize_t dims[rank];
    H5Sget_simple_extent_dims(space, dims, NULL);
    //The number of wavenumbers, momentum bins, and multipoles
    bolt->Nk = dims[0];
    bolt->Nq = dims[1];
    bolt->Nl = dims[2];

    message("Reading ncdm perturbation (%zu wavenumbers, %zu momentum bins, %zu multipoles)",
    bolt->Nk, bolt->Nq, bolt->Nl);

    //Create a temporary array to read the data
    float Psi[bolt->Nk][bolt->Nq][bolt->Nl];
    H5Dread(h_data, H5T_NATIVE_FLOAT, space, space, H5P_DEFAULT, Psi);
    H5Dclose(h_data);

    //Allocate memory in the main program
    bolt->Psi = malloc(bolt->Nk * bolt->Nq * bolt->Nl * sizeof(float));
    bolt->Psi_dot = malloc(bolt->Nk * bolt->Nq * bolt->Nl * sizeof(float));
    bolt->k_bins = malloc(bolt->Nk * sizeof(float));
    bolt->q_bins = malloc(bolt->Nq * sizeof(float));

    //Transfer the data to bolt->Psi
    for (size_t ik=0; ik<bolt->Nk; ik++) {
        for (size_t iq=0; iq<bolt->Nq; iq++) {
            for (size_t il=0; il<bolt->Nl; il++) {
                bolt->Psi[ik + iq*bolt->Nk + il*bolt->Nk*bolt->Nq]
                                = Psi[ik][iq][il];
            }
        }
    }

    //Next, load the wavenumber bins (1 column of length Nk)
    h_data = H5Dopen(h_file, "wavenumber_bins (k)", H5P_DEFAULT);
    space = H5Dget_space(h_data);

    //Find the number of bins
    hsize_t k_dims[1];
    H5Sget_simple_extent_dims(space, k_dims, NULL);
    if (k_dims[0] != bolt->Nk) {
        error("Incorrect dataset dimensions for the wavenumber bins.");
    }

    //Read the data
    H5Dread(h_data, H5T_NATIVE_FLOAT, space, space, H5P_DEFAULT, bolt->k_bins);
    H5Dclose(h_data);

    //Next, load the momentum bins (1 column of length Nq)
    h_data = H5Dopen(h_file, "momentum_bins (q)", H5P_DEFAULT);
    space = H5Dget_space(h_data);

    //Find the number of bins
    hsize_t q_dims[1];
    H5Sget_simple_extent_dims(space, q_dims, NULL);
    if (q_dims[0] != bolt->Nq) {
        error("Incorrect dataset dimensions for the momentum bins.");
    }

    //Read the data
    H5Dread(h_data, H5T_NATIVE_FLOAT, space, space, H5P_DEFAULT, bolt->q_bins);
    H5Dclose(h_data);

    //Close the file
    H5Fclose(h_file);
}
