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

/* We also use the interpolation methods from the Eagle cooling codebase */
#include "../cooling/EAGLE/interpolate.h"

void boltz_init(struct boltz *bolt, struct pm_mesh *mesh) {
    //Initialize mesh reference
    bolt->mesh = mesh;
    //Initialize the power spectrum data
    bolt->bins = 20;
    //The perturbations has not yet been refactored
    bolt->refactored = 0;

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


void boltz_update_phi(struct boltz *bolt) {

    /* Some useful constants */
    const double box_size = bolt->mesh->dim[0];
    const int N = bolt->mesh->N;
    const int N_half = N / 2;
    const double cell_fac = N / box_size;

    /* Use the memory allocated for the potential to temporarily store rho */
    double* restrict rho = bolt->mesh->potential;
    if (rho == NULL) error("Error allocating memory for density mesh");

    fftw_complex* restrict frho =
        (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N * N * (N_half + 1));
    if (frho == NULL)
      error("Error allocating memory for transform of density mesh");
    memuse_log_allocation("fftw_frho", frho, 1,
                          sizeof(fftw_complex) * N * N * (N_half + 1));

  /* Prepare the FFT library */
  fftw_plan forward_plan = fftw_plan_dft_r2c_3d(
      N, N, N, rho, frho, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);


      /* Fourier transform to go to magic-land */
      fftw_execute(forward_plan);

      //Normalization
      for (int x=0; x<N; x++) {
          for (int y=0; y<N; y++) {
              for (int z=0; z<=N_half; z++) {
                  frho[row_major_id_half_periodic(x, y, z, N)][0] /= cell_fac;
                  frho[row_major_id_half_periodic(x, y, z, N)][1] /= cell_fac;
              }
          }
      }

      fftw_destroy_plan(forward_plan);
      memuse_log_allocation("fftw_frho", frho, 0, 0);
      fftw_free(frho);

      //Calculate the power spectrum
  	bolt->k_in_bins = (double*) malloc(sizeof(double)*bolt->bins);
  	bolt->power_in_bins = (double*) malloc(sizeof(double)*bolt->bins);
  	bolt->obs_in_bins = (int*) malloc(sizeof(int)*bolt->bins);

  	for (size_t i=0; i<bolt->bins; i++) {
  		bolt->k_in_bins[i] = 0;
  		bolt->power_in_bins[i] = 0;
  		bolt->obs_in_bins[i] = 0;
  	}

    calc_cross_powerspec(N,box_size,frho,frho,bolt->bins,bolt->k_in_bins,
                        bolt->power_in_bins,bolt->obs_in_bins,CLOUD_IN_CELL);

    if (!bolt->refactored) {
        boltz_refactor_bins(bolt,bolt->bins);
    }
}

//Export the power spectrum table
void boltz_export_phi(struct boltz *bolt, const char *fname) {
    FILE *of = fopen(fname,"w");
    fprintf(of,"k(Mpc) P(k) observations\n");
    for (size_t i=0;i<bolt->bins;i++) {
        if (bolt->obs_in_bins[i]>0) {
            fprintf(of,"%f %f %i\n",bolt->k_in_bins[i],bolt->power_in_bins[i],bolt->obs_in_bins[i]);
        }
    }
    fclose(of);
}


//Restructure the neutrino pertubration
void boltz_refactor_bins(struct boltz *bolt, size_t bins) {
    //After this, we have refactored!
    bolt->refactored = 1;

    //Define some constants
    size_t old_Nk = bolt->Nk; //this will be changed
    size_t Nl = bolt->Nl;
    size_t Nq = bolt->Nq;

    //Recompute the power spectrum
    bolt->bins = bins;
    boltz_update_phi(bolt);

    double k_min = bolt->k_in_bins[0];
    double k_max = bolt->k_in_bins[bins-1];

    message("k_min = %f, k_max = %f", k_min, k_max);

    float new_Psi[bins*Nl*Nq];

    double k = 0;
    double k_glb = 0; //the greatest lower bound of (k,inf) in k_bin
    size_t ik_glb = 0; //the corresponding index in k_bin

    //Find k_glb for each k in k_in_bins
    for (size_t bin=0; bin<bins; bin++) {
        //Handle bins without any obervations
        if (bolt->obs_in_bins[bin] > 0) {
            k = bolt->k_in_bins[bin];
        }

        //Find the greatest lower bound of (k,inf) in k_bin q
        for (size_t ik=ik_glb; ik<bolt->Nk; ik++) {
            if (ik==ik_glb || bolt->k_bins[ik]<k) {
                k_glb = bolt->k_bins[ik];
                ik_glb = ik;
            }
        }

        //Now interpolate Psi using the new k's (defined by k_in_bins)
        for (size_t iq=0; iq<Nq; iq++) {
            for (size_t il=0; il<Nl; il++) {
                new_Psi[bin + iq*bins + il*Nq*bins] =
                    interpolation_3d(bolt->Psi, il, iq, ik_glb, 0, 0, k-k_glb,
                                        Nl, Nq, old_Nk);
            }
        }

        /**
        message("%zu %f %f %f %zu %f %f %f", bin, k, k_glb, bolt->k_bins[ik_glb+1], ik_glb, new_Psi[bin + 3*bins + 2*bins*Nq],
                                                                bolt->Psi[ik_glb + 3*old_Nk + 2*old_Nk*Nq],
                                                                bolt->Psi[ik_glb+1 + 3*old_Nk + 2*old_Nk*Nq]);
        */
    }

    //Update the length of k_bins
    size_t Nk = bins;
    bolt->Nk = Nk;

    //Update the data tables with the newly refactored tables
    free(bolt->Psi);
    free(bolt->Psi_dot);
    free(bolt->k_bins);
    bolt->Psi = malloc(bolt->Nk * bolt->Nq * bolt->Nl * sizeof(float));
    bolt->Psi_dot = malloc(bolt->Nk * bolt->Nq * bolt->Nl * sizeof(float));
    bolt->k_bins = malloc(bolt->Nk * sizeof(float));

    for (size_t ik=0; ik<Nk; ik++) {
        for (size_t iq=0; iq<Nq; iq++) {
            for (size_t il=0; il<Nl; il++) {
                bolt->Psi[ik + iq*Nk + il*Nk*Nq] = new_Psi[ik + iq*Nk + il*Nk*Nq];
            }
        }
    }

    message("Refactored the neutrino perturbation (%zu wavenumbers, %zu momentum bins, %zu multipoles)",
    bolt->Nk, bolt->Nq, bolt->Nl);

}
