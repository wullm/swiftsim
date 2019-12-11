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
#include "neutrino.h"

/* We also use the interpolation methods from the Eagle cooling codebase */
#include "../cooling/EAGLE/interpolate.h"

void boltz_init(struct boltz *bolt, struct swift_params *params, const struct engine *e) {
    //Initialize parameters
    bolt->m_nu = e->cosmology->m_nu;
    //Initialize parameters
    bolt->T_nu = parser_get_param_double(params, "Cosmology:T_nu");
    //Time unit used in the Boltzmann initial conditions
    bolt->boltz_time_unit = parser_get_param_double(params, "Boltzmann:boltz_time_unit");
    //The user-specified number of k-bins used in the power spectrum calculation
    bolt->num_of_k_bins = parser_get_opt_param_int(params, "Boltzmann:k_bins", BOLTZ_DEFAULT_BINS);
    //Number of time steps over which to average the power spectrum measurements
    bolt->PS_lags = 50;
    //Allocate memory for the power spectrum data
    bolt->powerSpec.k_in_bins = calloc(bolt->num_of_k_bins*bolt->PS_lags, sizeof(double));
    bolt->powerSpec.power_in_bins = calloc(bolt->num_of_k_bins*bolt->PS_lags, sizeof(double));
    bolt->powerSpec.obs_in_bins = calloc(bolt->num_of_k_bins*bolt->PS_lags, sizeof(int));
    //Allocate memory for the cosmological times at which the power spectrum was recorded
    bolt->record_times = calloc(bolt->PS_lags, sizeof(double));
    //Allocate memory for the gravitational potential vector and its derivative
    bolt->d_cdm = calloc(bolt->num_of_k_bins, sizeof(double));
    bolt->d_cdm_prime = calloc(bolt->num_of_k_bins, sizeof(double));
    bolt->d_cdm_prime_error = calloc(bolt->num_of_k_bins, sizeof(double));

    // //Zero the power spectrum arrays
    // for (size_t i=0; i<bolt->num_of_k_bins; i++) {
    //     bolt->powerSpec.k_in_bins[i] = 0;
    //     bolt->powerSpec.power_in_bins[i] = 0;
    //     bolt->powerSpec.obs_in_bins[i] = 0;
    // }


    //Read the file name of the hdf5 neutrino perturbation file
    char perturbFName[200] = "";
    parser_get_param_string(params, "Boltzmann:file_name", perturbFName);

    //Open and load the file with the neutrino perturbation (e.g. from CLASS)
    boltz_load_nu_perturb(bolt, perturbFName);

    message("Psi_0 %.10e", bolt->Psi[0] );

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

    //Verify that the primordial grid has the same grid size as the gravity mesh
    if (bolt->primordial_grid_N != (size_t) e->mesh->N) {
       error("Primordial grid is not the same size as the gravity mesh %zu != %zu.", bolt->primordial_grid_N, (size_t) e->mesh->N);
   }

message("Psi_0 %.10e", bolt->Psi[0] );
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
    bolt->Psi_prime = malloc(bolt->Nk * bolt->Nq * bolt->Nl * sizeof(float));
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


        for (size_t ik=0; ik<bolt->Nk; ik++) {
            for (size_t iq=0; iq<bolt->Nq; iq++) {
                double k = bolt->k_bins[ik]; //Mpc
                double q = bolt->q_bins[iq];

                //For each multipole
                for (size_t l=0; l<bolt->Nl; l++) {
                    int bix = ik + iq*bolt->Nk;
                    int ix = bix + l*bolt->Nq*bolt->Nk;

                    if (iq==0) {
                        message("%f %f %zu %.10e", q, k, l, bolt->Psi[ix]);
                    }
                }
            }
        }
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

    //The grid must be cubic
    if (grid_dims[0] != grid_dims[1] || grid_dims[0] != grid_dims[2]) {
        error("Primordial grid is not cubic.");
    }

    size_t N = grid_dims[0];
    bolt->primordial_grid_N = N;

    //Create a temporary array to read the data
    float grf[N][N][N];
    H5Dread(h_data, H5T_NATIVE_FLOAT, space, space, H5P_DEFAULT, grf);
    H5Dclose(h_data);

    //Allocate memory in the main program
    bolt->primordial_grid = malloc(N*N*N * sizeof(double));

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


void boltz_update_powerspec(struct boltz *bolt, const struct engine *e, fftw_complex* restrict frho) {
message("Psi_0 %.10e", bolt->Psi[0] );

    /* Some useful constants */
    const double box_size = e->mesh->dim[0];
    const int N = e->mesh->N;
    const int N_half = N / 2;
    const double cell_fac = N / box_size;

    // /* Use the memory allocated for the potential to temporarily store rho */
    // double* restrict rho = e->mesh->potential;
    // if (rho == NULL) error("Error allocating memory for density mesh");
    //
    // fftw_complex* restrict frho = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N * N * (N_half + 1));
    // if (frho == NULL)  error("Error allocating memory for transform of density mesh");
    // memuse_log_allocation("fftw_frho", frho, 1, sizeof(fftw_complex) * N * N * (N_half + 1));
    //
    // /* Prepare the FFT library */
    // fftw_plan forward_plan = fftw_plan_dft_r2c_3d(N, N, N, rho, frho, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
    //
    //
    // /* Fourier transform to go to magic-land */
    // fftw_execute(forward_plan);
    //

    //Normalization (this is the ordinary Fourier transform normalization)
    for (int x=0; x<N; x++) {
        for (int y=0; y<N; y++) {
            for (int z=0; z<=N_half; z++) {
                frho[row_major_id_half_periodic(x, y, z, N)][0] /= cell_fac*cell_fac*cell_fac;
                frho[row_major_id_half_periodic(x, y, z, N)][1] /= cell_fac*cell_fac*cell_fac;
            }
        }
    }

    //Divide by average density to get the overdensity field
    double rho_0 = e->total_mass/(box_size*box_size*box_size);
    for (int x=0; x<N; x++) {
        for (int y=0; y<N; y++) {
            for (int z=0; z<=N_half; z++) {
                frho[row_major_id_half_periodic(x, y, z, N)][0] /= rho_0;
                frho[row_major_id_half_periodic(x, y, z, N)][1] /= rho_0;
            }
        }
    }



    // /* Use the memory allocated for the potential to temporarily store rho */
    // double* restrict rho = (double*) malloc(N*N*N*sizeof(double));
    // if (rho == NULL) error("Error allocating memory for density mesh");
    //
    //     /* Prepare the FFT library */
    // fftw_plan forward_plan = fftw_plan_dft_c2r_3d(N, N, N, frho, rho, FFTW_ESTIMATE);
    //
    // /* Fourier transform to go to magic-land */
    // fftw_execute(forward_plan);
    //
    // //Normalization
    // for (int x=0; x<N; x++) {
    //     for (int y=0; y<N; y++) {
    //         for (int z=0; z<N; z++) {
    //             rho[x + y*N + z*N*N] /= box_size*box_size*box_size;
    //         }
    //     }
    // }
    //
    //
    // //Compute the total
    // double tot =0;
    // for (int x=0; x<N; x++) {
    //     for (int y=0; y<N; y++) {
    //         for (int z=0; z<N; z++) {
    //             // double a = frho[row_major_id_half_periodic(x, y, z, N)][0];
    //             // double b = frho[row_major_id_half_periodic(x, y, z, N)][1];
    //             // tot += a*a + b*b;
    //             tot += rho[x + y*N + z*N*N];
    //         }
    //     }
    // }
    //
    // message("Total %f", sqrt(tot));

    //
    // fftw_destroy_plan(forward_plan);
    // memuse_log_allocation("fftw_frho", frho, 0, 0);
    // fftw_free(frho);

    message("Box size %f Cell frac %f", box_size, cell_fac);

    // //Zero the power spectrum arrays
    // for (size_t i=0; i<bolt->num_of_k_bins * bolt->PS_lags; i++) {
    //     bolt->powerSpec.k_in_bins[i] = 0;
    //     bolt->powerSpec.power_in_bins[i] = 0;
    //     bolt->powerSpec.obs_in_bins[i] = 0;
    // }

    //Calculate d_cdm as the average over all the power spectrum columns
    for (size_t i=0; i<bolt->num_of_k_bins; i++) {
        double avg = 0;
        for (size_t j=0; j<bolt->PS_lags; j++) {
            avg += bolt->powerSpec.power_in_bins[i + j*bolt->num_of_k_bins];
        }
        avg *= 1. / bolt->PS_lags;

        //Convert the average from power spectrum to transfer function
        bolt->d_cdm[i] = sqrt(avg / pow(bolt->powerSpec.k_in_bins[i], 0.9667));

        // //Calculate the variance
        // double var = 0;
        // for (size_t j=0; j<bolt->PS_lags; j++) {
        //     var += pow(bolt->powerSpec.power_in_bins[i + j*bolt->num_of_k_bins] - avg,2);
        // }
        // var *= 1. / (bolt->PS_lags-1);
        //
        // bolt->d_cdm_prime[i] = sqrt(var);
    }

    // //Calculate the derivative from the historical power spectrum table
    // //The most recent calculations are on the left
    // for (size_t i=0; i<bolt->num_of_k_bins; i++) {
    //     double d_cdm_left = 0;
    //     double d_cdm_right = 0;
    //     for (size_t j=0; j<25; j++) {
    //         d_cdm_left += bolt->powerSpec.power_in_bins[i + j*bolt->num_of_k_bins];
    //         d_cdm_right += bolt->powerSpec.power_in_bins[i + (j+25+1)*bolt->num_of_k_bins];
    //     }
    //
    //     d_cdm_left *= 1. / 25;
    //     d_cdm_right *= 1. / 25;
    //
    //     bolt->d_cdm_prime[i] = (0.5*d_cdm_left - 0.5*d_cdm_right) / (e->time_step);
    //     bolt->d_cdm_prime[i] = (0.5*d_cdm_left - 0.5*d_cdm_right) / (1.);
    // }

    //Calculate the derivative from a simple linear fit
    //The most recent calculations are on the left
    for (size_t i=0; i<bolt->num_of_k_bins; i++) {
        double Sx = 0;
        double Sy = 0;
        double Sxx = 0;
        double Sxy = 0;
        double Syy = 0;
        size_t n = bolt->PS_lags;

        for (size_t j=0; j<bolt->PS_lags; j++) {
            double time = bolt->record_times[j];
            double power = sqrt(bolt->powerSpec.power_in_bins[i + j*bolt->num_of_k_bins] / pow(bolt->powerSpec.k_in_bins[i], 0.9667));

            Sy += power;
            Syy += power*power;
            Sx += time;
            Sxx += time*time;
            Sxy += time*power;

            // if (i==6) {
            //     message("%.10e %.10e", time, power);
            // }
        }


        double beta = (n*Sxy - Sx*Sy)/(n*Sxx-Sx*Sx);
        double sigma_eps_squared = (n*Syy - Sy*Sy - beta*beta*(n*Sxx-Sx*Sx))/(n*(n-2));
        double sigma_beta_squared = n*sigma_eps_squared/(n*Sxx-Sx*Sx);

        if (i==6)
        message("%.10e", beta);

        bolt->d_cdm_prime[i] = beta;
        bolt->d_cdm_prime_error[i] = sqrt(sigma_beta_squared);

        //Convert the time derivative to a conformal time derivative (using dt/dtau = a)
        bolt->d_cdm_prime[i] *= e->cosmology->a;
        bolt->d_cdm_prime_error[i] *= e->cosmology->a;
    }


    //Move all the power spectrum columns to the right
    for (size_t i=0; i<bolt->num_of_k_bins; i++) {
        for (size_t j=bolt->PS_lags-1; j>0; j--) {
            bolt->powerSpec.k_in_bins[i + j*bolt->num_of_k_bins] = bolt->powerSpec.k_in_bins[i + (j-1)*bolt->num_of_k_bins];
            bolt->powerSpec.power_in_bins[i + j*bolt->num_of_k_bins] = bolt->powerSpec.power_in_bins[i + (j-1)*bolt->num_of_k_bins];
            bolt->powerSpec.obs_in_bins[i + j*bolt->num_of_k_bins] = bolt->powerSpec.obs_in_bins[i + (j-1)*bolt->num_of_k_bins];
        }
    }

    for (size_t j=bolt->PS_lags-1; j>0; j--) {
        bolt->record_times[j] = bolt->record_times[j-1];
    }

    //Record the present time
    bolt->record_times[0] = e->time;

    // message("The mores the times %.10e %.10e %.10e %.10e %.10e", bolt->record_times[0], bolt->record_times[1], bolt->record_times[2], bolt->record_times[3], bolt->record_times[4]);

    //Zero the first column of the power spectrum arrays
    for (size_t i=0; i<bolt->num_of_k_bins; i++) {
        bolt->powerSpec.k_in_bins[i] = 0;
        bolt->powerSpec.power_in_bins[i] = 0;
        bolt->powerSpec.obs_in_bins[i] = 0;
    }

    //Calculate the power spectrum into the first column
    calc_cross_powerspec(N,box_size,frho,frho,bolt->num_of_k_bins,bolt->powerSpec.k_in_bins, bolt->powerSpec.power_in_bins,bolt->powerSpec.obs_in_bins,CLOUD_IN_CELL);

    //Check if the neutrino perturbation needs to be refactored
    if (bolt->Nk != bolt->num_of_k_bins) {
        boltz_refactor_bins(bolt, e, bolt->num_of_k_bins);
    }

    for (size_t i=0; i<bolt->Nk; i++) {
        message("%f %.10e", bolt->k_bins[i], bolt->Psi[i + 0*bolt->Nk*bolt->Nq]);
    }

    //Update d_cdm and d_cdm_prime
    // double G_newton = e->physical_constants->const_newton_G;
    // double normalization = 1./71.;
    // double normalization = e->cosmology->D_linear_growth/e->cosmology->D_linear_growth0;
    // message("Growth %f %f OL Om %f %f", e->cosmology->D_linear_growth, e->cosmology->D_linear_growth0, e->cosmology->Omega_lambda, e->cosmology->Omega_m);
    // message("Hubble %f", e->cosmology->H / e->cosmology->H0);
    // message("Or %f", sqrt(e->cosmology->Omega_r * pow(1+40,4)) / (e->cosmology->H / e->cosmology->H0));
    // message("And the step is %d", engine_current_step);
    // for (size_t i=0; i<bolt->Nk; i++) {
    //     // double k = bolt->powerSpec.k_in_bins[i];
    //     double old_d_cdm = bolt->d_cdm[i];
    //     double new_d_cdm = sqrt(bolt->powerSpec.power_in_bins[i]);
    //     // double new_d_cdm = sqrt(bolt->powerSpec.power_in_bins[i]/pow(k, 0.97))*(4*M_PI*G_newton/(k*k))*normalization;
    //
    //     bolt->d_cdm[i] = new_d_cdm;
    //     bolt->d_cdm_prime[i] = (new_d_cdm - old_d_cdm) / e->time_step;
    // }
}

//Export the power spectrum table
void boltz_export_powerspec(struct boltz *bolt, const char *fname) {
    FILE *of = fopen(fname,"w");
    fprintf(of,"k(Mpc) P(k) observations d_cdm(k) d_cdm_prime(k) d_cdm_prime_err(k)\n");
    for (size_t i=0;i<bolt->num_of_k_bins;i++) {
        if (bolt->powerSpec.obs_in_bins[i]>0) {
            fprintf(of,"%f %f %i %f %f %f\n",bolt->powerSpec.k_in_bins[i],bolt->powerSpec.power_in_bins[i],bolt->powerSpec.obs_in_bins[i],bolt->d_cdm[i],bolt->d_cdm_prime[i],bolt->d_cdm_prime_error[i]);
        }
    }
    fclose(of);
}


//Restructure the neutrino pertubration by chaning the number of k-bins
void boltz_refactor_bins(struct boltz *bolt, const struct engine *e, size_t num_of_k_bins) {
    //Define some constants
    size_t old_Nk = bolt->Nk; //this will be changed
    size_t Nl = bolt->Nl;
    size_t Nq = bolt->Nq;

    //Recompute the power spectrum if necessary
    if (num_of_k_bins != bolt->num_of_k_bins) {
        bolt->num_of_k_bins = num_of_k_bins;
        // boltz_update_d_cdm(bolt, e, e->mesh->potential);
        error("Try computing the potential first.\n");
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
    free(bolt->Psi_prime);
    free(bolt->k_bins);
    bolt->Psi = malloc(bolt->Nk * bolt->Nq * bolt->Nl * sizeof(float));
    bolt->Psi_prime = malloc(bolt->Nk * bolt->Nq * bolt->Nl * sizeof(float));
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

void boltz_step(struct boltz *bolt, const struct engine *e) {

    //Convert internal time units to Mpc/c used by CLASS
    // double Mpc_in_cm = 3.086e24;
    // double time_unit = e->physical_constants->const_speed_light_c * e->internal_units->UnitLength_in_cgs / Mpc_in_cm;

    message("%.10e", bolt->boltz_time_unit);

    //Useful constants
    double a = e->cosmology->a;
    size_t Nk = bolt->Nk;
    size_t Nq = bolt->Nq;
    int max_l = (int) bolt->Nl-1; //maximum multipole (accounting for l=0)
    double eV_in_kT_nu = e->physical_constants->const_electron_volt/(bolt->T_nu * e->physical_constants->const_boltzmann_k); //dimensionless quantity
    double m_nu = bolt->m_nu * eV_in_kT_nu/1.; //mass of the neutrino
    message("neutrino mass %.10e", m_nu);


    message("Psi_0 %.10e", bolt->Psi[0] );

    //We need to have performed enough steps to get a reliable estimate of h_dot
    if ((size_t) engine_current_step < bolt->PS_lags) {
        message("Still working on it %d < %zu", engine_current_step, bolt->PS_lags);
        return;
    }


    // if (ik>3)

    size_t substeps = 100000;
    for (size_t it=0; it<substeps; it++) {

        /**
        * The evolution of each (wavenumber k, momentum bin q) pair is completely
        * independent. Therefore, we could parallelize this.
        */
        for (size_t ik=0; ik<Nk; ik++) {
            for (size_t iq=0; iq<Nq; iq++) {
                double k = bolt->k_bins[ik]; //Mpc
                double q = bolt->q_bins[iq];
                double eps = sqrt(q*q + a*a*m_nu*m_nu);

                // message("==========%f %f", k, eps);

                //For each multipole
                for (int l=0; l<=max_l; l++) {
                    int bix = ik + iq*Nk;
                    int ix = bix + l*Nq*Nk;

                    //Use the fact that h' = -2*d_cdm' and eta' \approx 0
                    double h_prime = 0; //-2*bolt->d_cdm_prime[ik]*bolt->boltz_time_unit;
                    double eta_prime = 0;

                    // message("%f %.10e %.10e %.10e", k, q*k/eps, h_prime, bolt->dlogf0_dlogq[iq]);

                    //The first three multipoles have a unique evolution equation
                    if (l == 0) {
                        bolt->Psi_prime[ix] =
                            -q*k/eps * bolt->Psi[bix + 1*Nk*Nq]
                            +(1./6.) * h_prime * bolt->dlogf0_dlogq[iq];

                        if ((size_t) engine_current_step == bolt->PS_lags && iq == 0 && (it==0 || it==substeps-1))
                        message("%f %f %d %.10e %.10e %.10e %.10e", q, k, l, bolt->Psi[ix], bolt->Psi_prime[ix], 1./a* e->time_step / bolt->boltz_time_unit, bolt->Psi[bix + 1*Nk*Nq]);
                    } else if (l == 1) {
                        bolt->Psi_prime[ix] =
                            +q*k/(3.*eps) *(bolt->Psi[bix + 0*Nk*Nq] - 2*bolt->Psi[bix + 2*Nk*Nq]);

                        if ((size_t) engine_current_step == bolt->PS_lags && iq == 0 && (it==0 || it==substeps-1))
                        message("%f %f %d %.10e %.10e %.10e %.10e", q, k, l, bolt->Psi[ix], bolt->Psi_prime[ix], 1./a* e->time_step / bolt->boltz_time_unit, bolt->Psi[bix + 0*Nk*Nq]);
                    } else if (l == 2) {
                        bolt->Psi_prime[ix] =
                            +q*k/(5.*eps) *(2*bolt->Psi[bix + 1*Nk*Nq] - 3*bolt->Psi[bix + 3*Nk*Nq])
                            -(1./15. * h_prime + 2./5. * eta_prime) * bolt->dlogf0_dlogq[iq];

                        if ((size_t) engine_current_step == bolt->PS_lags && iq == 0 && (it==0 || it==substeps-1))
                        message("%f %f %d %.10e %.10e %.10e %.10e %.10e", q, k, l, bolt->Psi[ix], bolt->Psi_prime[ix], 1./a* e->time_step / bolt->boltz_time_unit, bolt->Psi[bix + 3*Nk*Nq], q*k/(5.*eps));
                    } else if (l <= max_l-1) {
                        bolt->Psi_prime[ix] = q*k/(eps * (2*l+1)) * (l*bolt->Psi[bix + (l-1)*Nk*Nq] - (l+1)*bolt->Psi[bix + (l+1)*Nk*Nq]);

                        if ((size_t) engine_current_step == bolt->PS_lags && iq == 0 && (it==0 || it==substeps-1))
                        message("%f %f %d %.10e %.10e %.10e %.10e", q, k, l, bolt->Psi[ix], bolt->Psi_prime[ix], 1./a* e->time_step / bolt->boltz_time_unit, bolt->Psi[bix + (l+1)*Nk*Nq]);
                    } else {
                        //The last multipole is truncated (easy in flat cosmology)
                        bolt->Psi_prime[ix] = 0;

                        if ((size_t) engine_current_step == bolt->PS_lags && iq == 0 && (it==0 || it==substeps-1))
                        message("%f %f %d %.10e %.10e %.10e %.10e", q, k, l, bolt->Psi[ix], bolt->Psi_prime[ix], 1./a* e->time_step / bolt->boltz_time_unit, bolt->Psi[bix + (l-1)*Nk*Nq]);

                        if (e->cosmology->Omega_k > 0) {
                            error("Neutrino evolution not yet implemented for non-flat cosmologies.");
                        }
                    }

                    //Update the multipole
                    // bolt->Psi[ix] += (bolt->Psi_prime[ix] / a) * e->time_step / bolt->boltz_time_unit;
                }
            }
        }

        //Now update the multipoles
        for (size_t ik=0; ik<Nk; ik++) {
            for (size_t iq=0; iq<Nq; iq++) {
                for (int l=0; l<=max_l; l++) {
                    int bix = ik + iq*Nk;
                    int ix = bix + l*Nq*Nk;

                    bolt->Psi[ix] += (bolt->Psi_prime[ix] / a) * (e->time_step/substeps) / bolt->boltz_time_unit;
                }
            }
        }
    }
}
