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

void boltz_init(struct boltz *bolt, struct swift_params *params, struct engine *e) {
    //Initialize mesh reference
    bolt->mesh = e->mesh;
    //Initialize parameters
    bolt->m_nu = e->cosmology->m_nu;
    //The user-specified number of k-bins used in the power spectrum calculation
    bolt->num_of_k_bins = parser_get_opt_param_int(params, "Boltzmann:k_bins", BOLTZ_DEFAULT_BINS);
    //Allocate memory for the power spectrum data
    bolt->powerSpec.k_in_bins = (double*) malloc(sizeof(double)*bolt->num_of_k_bins);
    bolt->powerSpec.power_in_bins = (double*) malloc(sizeof(double)*bolt->num_of_k_bins);
    bolt->powerSpec.obs_in_bins = (int*) malloc(sizeof(int)*bolt->num_of_k_bins);
    //Allocate memory for the gravitational potential vector and its derivative
    bolt->phi = (double*) malloc(sizeof(double)*bolt->num_of_k_bins);
    bolt->phi_dot = (double*) malloc(sizeof(double)*bolt->num_of_k_bins);

    //Read the file name of the hdf5 neutrino perturbation file
    char perturbFName[200] = "";
    parser_get_param_string(params, "Boltzmann:file_name", perturbFName);

    //Open and load the file with the neutrino perturbation (e.g. from CLASS)
    boltz_load_nu_perturb(bolt, perturbFName);

    //Now compute the logarithmic derivative of the distribution function f0
    bolt->f0_prefactor = 2.0 / e->physical_constants->const_planck_h;
    message("Distribution function pre-factor %.10e", bolt->f0_prefactor);

    //Allocate memory for the derivative
    bolt->dlogf0_dlogq = malloc(bolt->Nq * sizeof(double));

    //Numerically differentiate the 0th order distribution function
    double precision = 1e-8;
    for (size_t iq=0; iq<bolt->Nq; iq++) {
        double q = bolt->q_bins[iq];
        double fl = f0_nu(q-precision);
        double fr = f0_nu(q+precision);
        double df_dq = (fr-fl)/(2*precision);
        bolt->dlogf0_dlogq[iq] = df_dq*q/fl;
        // message("q = %.10e, dln(f)/dln(q) = %.10e %.10e", q, df_dq, dlogf_dlogq);
    }

    //Read the file name of the hdf5 gaussian random field file
    char fieldFName[200] = "";
    parser_get_param_string(params, "Boltzmann:field_file_name", fieldFName);

    //Open and load the file with the primordial Gaussian field
    boltz_load_primordial_field(bolt, fieldFName);

    //Verify that the physical dimensions of the primordial field match the cosmology
    for (int i=0; i<3; i++) {
        if (fabs(bolt->primordial_dims[i] - e->s->dim[i]) > 1e-2) {
            error("Dimensions[%i] of primordial field do not agree with engine->space.", i);
        }
    }
}

void boltz_load_nu_perturb(struct boltz *bolt, const char *fname) {
    //Open the file containing the neutrino perturbation (e.g. from CLASS)
    const hid_t h_file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (h_file < 0) {
        error("Error opening the multipoles file.");
    }

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

    //Are there enough multipoles?
    if (bolt->Nl < 4) {
        error("Not enough neutrino perturbation multipoles.");
    }

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

void boltz_load_primordial_field(struct boltz *bolt, const char *fname) {
    //Open the file containing the primordial fluctuation field
    const hid_t field_file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (field_file < 0) {
        error("Error opening the primordial field file.");
    }

    //Open the header group
    hid_t h_grp = H5Gopen(field_file, "/Header", H5P_DEFAULT);
    if (h_grp < 0) {
        error("Error while opening file header\n");
    }

    //Read the physical dimensions of the box
    const hid_t hid_bsz = H5Aexists(h_grp, "BoxSize");
    if (hid_bsz < 0) {
        error("Error while testing existance of 'BoxSize' attribute");
    }

    double field_dims[3];
    io_read_attribute(h_grp, "BoxSize", DOUBLE, field_dims);
    bolt->primordial_dims = (double*) malloc(3 * sizeof(double));
    for (int i=0; i<3; i++) {
        bolt->primordial_dims[i] = field_dims[i];
    }

    //Now load the actual grid
    h_grp = H5Gopen(field_file, "/Field", H5P_DEFAULT);
    if (h_grp < 0) {
        error("Error while opening field group\n");
    }

    hid_t h_data = H5Dopen(h_grp, "GaussianRandomField", H5P_DEFAULT);
    hid_t space = H5Dget_space(h_data);

    //The number of dimensions in the dataset (expected 3)
    const int rank = H5Sget_simple_extent_ndims(space);
    if (rank != 3) {
        error("Incorrect dataset dimensions for primordial field.");
    }

    //Find the extent of each dimension (the grid size; not physical size)
    hsize_t grid_dims[rank];
    H5Sget_simple_extent_dims(space, grid_dims, NULL);

    //The grid must be cubic and have the same dimensions as the gravity mesh
    if (grid_dims[0] != grid_dims[1] || grid_dims[0] != grid_dims[2]) {
        error("Primordial grid is not cubic.");
    } else if (grid_dims[0] != (size_t) bolt->mesh->N) {
        error("Primordial grid is not the same size as the gravity mesh %zu != %zu.", (size_t) grid_dims[0], (size_t) bolt->mesh->N);
    }

    size_t N = bolt->mesh->N;

    //Create a temporary array to read the data
    float grf[N][N][N];
    H5Dread(h_data, H5T_NATIVE_FLOAT, space, space, H5P_DEFAULT, grf);
    H5Dclose(h_data);

    //Allocate memory in the main program
    bolt->primordial_grid = malloc(N*N*N * sizeof(float));

    //Transfer the data to bolt->primordial_grid
    for (size_t i=0; i<N; i++) {
        for (size_t j=0; j<N; j++) {
            for (size_t k=0; k<N; k++) {
                bolt->primordial_grid[k + j*N + i*N*N] = grf[i][j][k];
            }
        }
    }

    message("The primordial field has dimensions %fx%fx%f", field_dims[0], field_dims[1], field_dims[2]);
    message("The primordial grid has dimensions %zux%zux%zu", (size_t) grid_dims[0], (size_t) grid_dims[1], (size_t) grid_dims[2]);

    //Close the file
    H5Fclose(field_file);
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

    fftw_complex* restrict frho = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N * N * (N_half + 1));
    if (frho == NULL)  error("Error allocating memory for transform of density mesh");
    memuse_log_allocation("fftw_frho", frho, 1, sizeof(fftw_complex) * N * N * (N_half + 1));

    /* Prepare the FFT library */
    fftw_plan forward_plan = fftw_plan_dft_r2c_3d(N, N, N, rho, frho, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);


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

    //Zero the power spectrum arraus
    for (size_t i=0; i<bolt->num_of_k_bins; i++) {
        bolt->powerSpec.k_in_bins[i] = 0;
        bolt->powerSpec.power_in_bins[i] = 0;
        bolt->powerSpec.obs_in_bins[i] = 0;
    }

    //Calculate the power spectrum
    calc_cross_powerspec(N,box_size,frho,frho,bolt->num_of_k_bins,bolt->powerSpec.k_in_bins, bolt->powerSpec.power_in_bins,bolt->powerSpec.obs_in_bins,CLOUD_IN_CELL);

    //Check if the neutrino perturbation needs to be refactored
    if (bolt->Nk != bolt->num_of_k_bins) {
        boltz_refactor_bins(bolt,bolt->num_of_k_bins);
    }
}

//Export the power spectrum table
void boltz_export_phi(struct boltz *bolt, const char *fname) {
    FILE *of = fopen(fname,"w");
    fprintf(of,"k(Mpc) P(k) observations\n");
    for (size_t i=0;i<bolt->num_of_k_bins;i++) {
        if (bolt->powerSpec.obs_in_bins[i]>0) {
            fprintf(of,"%f %f %i\n",bolt->powerSpec.k_in_bins[i],bolt->powerSpec.power_in_bins[i],bolt->powerSpec.obs_in_bins[i]);
        }
    }
    fclose(of);
}


//Restructure the neutrino pertubration by chaning the number of k-bins
void boltz_refactor_bins(struct boltz *bolt, size_t num_of_k_bins) {
    //Define some constants
    size_t old_Nk = bolt->Nk; //this will be changed
    size_t Nl = bolt->Nl;
    size_t Nq = bolt->Nq;

    //Recompute the power spectrum if necessary
    if (num_of_k_bins != bolt->num_of_k_bins) {
        bolt->num_of_k_bins = num_of_k_bins;
        boltz_update_phi(bolt);
    }

    //The range of k values found when computing the power spectrum with "num_of_k_bins" bins
    double k_min = bolt->powerSpec.k_in_bins[0];
    double k_max = bolt->powerSpec.k_in_bins[num_of_k_bins-1];

    message("k_min = %f, k_max = %f", k_min, k_max);

    //We will interpolate the multi-dimensional array Psi along the k-dimension
    float new_Psi[num_of_k_bins*Nl*Nq];

    //For each k in the new table, find the largest k' in the old table, s.t. k'<k
    double k_new = 0;
    double k_old = 0; //the greatest lower bound of (k,inf) in k_bin
    size_t ik_old = 0; //the corresponding index in k_bin

    //Find k_old for each k_new in k_in_bins
    for (size_t bin=0; bin<num_of_k_bins; bin++) {
        if (bolt->powerSpec.obs_in_bins[bin] > 0) {
            k_new = bolt->powerSpec.k_in_bins[bin];

            //Find the greatest lower bound of (k,inf) in k_bin q
            for (size_t ik=ik_old; ik<bolt->Nk; ik++) {
                if (ik==ik_old || bolt->k_bins[ik]<k_new) {
                    k_old = bolt->k_bins[ik];
                    ik_old = ik;
                }
            }
        }

        //Now interpolate Psi using the new k's (defined by k_in_bins)
        for (size_t iq=0; iq<Nq; iq++) {
            for (size_t il=0; il<Nl; il++) {
                new_Psi[bin + iq*num_of_k_bins + il*Nq*num_of_k_bins] =
                    interpolation_3d(bolt->Psi, il, iq, ik_old, 0, 0, k_new-k_old,
                                        Nl, Nq, old_Nk);
            }
        }

        /**
        message("%zu %f %f %f %zu %f %f %f", bin, k, k_old, bolt->k_bins[ik_old+1], ik_old, new_Psi[bin + 3*num_of_k_bins + 2*num_of_k_bins*Nq],
                                                                bolt->Psi[ik_old + 3*old_Nk + 2*old_Nk*Nq],
                                                                bolt->Psi[ik_old+1 + 3*old_Nk + 2*old_Nk*Nq]);
        */
    }

    //Update the length of the neutrino perturbation along the k dimension
    size_t Nk = num_of_k_bins;
    bolt->Nk = Nk;

    //Update the data tables with the newly refactored tables
    free(bolt->Psi);
    free(bolt->Psi_dot);
    free(bolt->k_bins);
    bolt->Psi = malloc(bolt->Nk * bolt->Nq * bolt->Nl * sizeof(float));
    bolt->Psi_dot = malloc(bolt->Nk * bolt->Nq * bolt->Nl * sizeof(float));
    bolt->k_bins = malloc(bolt->Nk * sizeof *bolt->k_bins);

    for (size_t ik=0; ik<Nk; ik++) {
        for (size_t iq=0; iq<Nq; iq++) {
            for (size_t il=0; il<Nl; il++) {
                bolt->Psi[ik + iq*Nk + il*Nk*Nq] = new_Psi[ik + iq*Nk + il*Nk*Nq];
            }
        }
        bolt->k_bins[ik] = bolt->powerSpec.k_in_bins[ik];
    }

    message("Refactored the neutrino perturbation (%zu wavenumbers, %zu momentum bins, %zu multipoles)",
    bolt->Nk, bolt->Nq, bolt->Nl);

}

void boltz_step(struct boltz *bolt, struct engine *e) {
    //Useful constants
    double a = e->cosmology->a;
    size_t Nk = bolt->Nk;
    size_t Nq = bolt->Nq;
    int max_l = (int) bolt->Nl; //maximum multipole
    double m_nu = bolt->m_nu; //mass of the neutrino


    /**
    * The evolution of each (wavenumber k, momentum bin q) pair is completely
    * independent. Therefore, we could parallelize this.
    */
    for (size_t ik=0; ik<Nk; ik++) {
        for (size_t iq=0; iq<Nq; iq++) {
            double k = bolt->k_bins[ik]; //Mpc
            double q = bolt->q_bins[iq];
            double eps = sqrt(q*q + a*a*m_nu*m_nu);

            //For each multipole
            for (int l=0; l<max_l; l++) {
                int bix = ik + iq*Nk;
                int ix = bix + l*Nq*Nk;

                //The first three multipoles have a unique evolution equation
                if (l == 0) {
                    bolt->Psi_dot[ix] =
                        -q*k/eps * bolt->Psi[bix + 1*Nk*Nq]
                        +(1./6.) * bolt->phi_dot[ik] * bolt->dlogf0_dlogq[iq];
                } else if (l == 1) {
                    bolt->Psi_dot[ix] =
                        +q*k/(3.*eps) *(bolt->Psi[bix + 0*Nk*Nq] - 2*bolt->Psi[bix + 2*Nk*Nq]);
                } else if (l == 2) {
                    bolt->Psi_dot[ix] =
                        +q*k/(5.*eps) *(2*bolt->Psi[bix + 1*Nk*Nq] - 3*bolt->Psi[bix + 3*Nk*Nq])
                        -(1./15. + 2./5.) * bolt->phi_dot[ik] * bolt->dlogf0_dlogq[iq];
                } else if (l < max_l-1) {
                    bolt->Psi_dot[ix] = q*k/(eps * (2*l+1)) * (l*bolt->Psi[bix + (l-1)*Nk*Nq] - (l+1)*bolt->Psi[bix + (l+1)*Nk*Nq]);
                } else {
                    //The last multipole is truncated (easy in flat cosmology)
                    bolt->Psi_dot[ix] = 0;

                    if (e->cosmology->Omega_k > 0) {
                        error("Neutrino evolution not yet implemented for non-flat cosmologies.");
                    }
                }

                //Update the multipole
                bolt->Psi[ix] += bolt->Psi_dot[ix] * e->time_step;
            }
        }
    }
    // if (ik>3)
    // message("Psi_0 %.10e", bolt->Psi[2] );
}
