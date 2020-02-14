/**
 * @icgen.cpp
 * @author  Willem Elbers <whe@willemelbers.com>
 * @version 1.0
 * @license MIT (https://opensource.org/licenses/MIT)
 *
 * @section DESCRIPTION
 *
 * Generate initial conditions for a cosmological N-body simulation.
 */

#include "config.h"
#include "icgen.h"

#include "create_grf.h"
#include "cosmo.h"
#include "read_transfer.h"
#include "sampler.h"
#include "infini_lpt.h"
#include "H5Cpp.h"

int main() {
    //Start the clock
    auto time = std::chrono::seconds(std::time(NULL));
    long int unix_time = std::chrono::milliseconds(time).count();

    std::cout << "Welcome to the initial condition generator." << std::endl;
    std::cout << "Provide a random seed (enter 0 for current time): ";

    long int seed;
    std::cin >> seed;

    if (!seed) {
        seed = unix_time;
        std::cout << "Using timestamp as seed: " << seed << "." << std::endl << std::endl;
    } else {
        std::cout << "Using seed: " << seed << "." << std::endl << std::endl;
    }

    std::default_random_engine oracle;
    oracle.seed(seed);
	std::normal_distribution<double> Gaussian(0.0, 1.0); //standard normal

    const int N = GRID_WIDTH;
    const double box_len = BOX_WIDTH; //Mpc
    const double delta_k = 2*M_PI/box_len; //Mpc^-1
    const double box_volume = box_len*box_len*box_len; //Mpc^3
    const double z_start = Z_START; //starting redshift

    //The primordial Gaussian random field and its Fourier transform
    double *primordial_box = (double*) fftw_malloc(sizeof(double)*N*N*N);
    fftw_complex *k_box = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*N*(N/2+1));
    fftw_plan r2c_plan  = fftw_plan_dft_r2c_3d(N, N, N, primordial_box, k_box, FFTW_ESTIMATE);

    std::cout << "PHASE 0A - Cosmology check" << std::endl;
    std::cout << "1) Starting redshift z = " << z_start << ", a = " << a_scale_factor_of_z(z_start) << "." << std::endl;
    std::cout << "2) Cosmology table written to " << std::string(OUTPUT_DIR) << "cosmology.txt." << std::endl;
    std::cout << "3) Omega_m = " << Omega_m << ", Omega_L = " << Omega_L << ", Omega_b = " << Omega_b << "," << std::endl;
    std::cout << "   Omega_nu = " << Omega_nu << ", Omega_g = " << Omega_r << "." << std::endl;
    std::cout << "4) Note that Omega_nu is not part of Omega_m." << std::endl;
    std::cout << std::endl;

    test_cosmology(std::string(OUTPUT_DIR) + "cosmology.txt");

    //Compute the relative contributions of cdm and baryons at the starting redshift
    std::cout << "PHASE 0B - Computing contributions to cold component" << std::endl;

    double weight_b = Omega_b / Omega_m;
    double weight_cdm = (Omega_m - Omega_b) / Omega_m;
    std::cout << "Weight of CDM in cold component: " << weight_cdm << "." << std::endl;
    std::cout << "Weight of baryons in cold component: " << weight_b << "." << std::endl;
    std::cout << std::endl;


    std::cout << "PHASE 0C - Reading in transfer function files" << std::endl;

    //Prepare an indexed search table for the transfer functions
    TF_index = (float*) malloc(TF_I_max * sizeof(float));

    //Neutrino and CDM density Transfer function data (loaded from CLASS output)
    read_transfer(TF_ks, TF_T_rho_cdm, TF_T_rho_nu, TF_T_rho_b, TF_T_rho_cb, weight_cdm, weight_b,
                    TF_T_theta_cdm, TF_T_theta_nu, TF_T_theta_b, TF_T_theta_cb);
    make_index_table(TF_I_max, TF_index, TF_ks, &log_k_min, &log_k_max);

	//Export transfer functions
	std::ofstream of(std::string(OUTPUT_DIR) + "transfer_functions.txt");
	of << "k(1/Mpc);T_cdm;T_nu;T_b;T_cb;T_theta_cdm;T_theta_nu;T_theta_b;T_theta_cb\n";

	for (int i = 0; i < TF_ks.size(); i++) {
		of << TF_ks[i] << ";" << TF_T_rho_cdm[i] << ";" << TF_T_rho_nu[i] << ";" << TF_T_rho_b[i] << ";" << TF_T_rho_cb[i] << ";" << TF_T_theta_cdm[i] << ";" << TF_T_theta_nu[i] << ";" << TF_T_theta_b[i] << ";" << TF_T_theta_cb[i] << std::endl;
	}

	of.close();

    std::cout << "1) Interpreted transfer functions exported to " << std::string(OUTPUT_DIR) << "transfer_functions.txt." << std::endl;
    std::cout << std::endl;

    std::cout << "PHASE 1A - Generating a Gaussian random field" << std::endl;

    //Generate a Gaussian random field
    generate_grf(oracle, k_box, N, box_len, sigma_func_cdm);

    std::cout << "1) Done with generating the box (" << N << "^3 complex numbers)." << std::endl;
    std::cout << std::endl;


    //We still need to normalize the power spectrum
    const double R_filter = 8/h; //Mpc
    double integrated_sigma_8 = integrate_sigma_R(N, box_len, R_filter, sigma_func_cdm);

    double global_PS_normalization = 1;

    if (NORMALIZATION_METHOD == NORM_CMB) {
        std::cout << "PHASE 1B - Normalizing the random field using A_s*k_pivot^-n_s" << std::endl;
        std::cout << "1a) The pivot scale is k = " << PIVOT_SCALE << " Mpc^-1." << std::endl;
        std::cout << "1b) The spectral index is n_s = " << N_S << "." << std::endl;
        std::cout << "1c) The normalization is A_s = " << A_S << "." << std::endl;

        global_PS_normalization = sqrt(A_S * pow(1.0 / PIVOT_SCALE, N_S));
    } else if (NORMALIZATION_METHOD == NORM_SIGMA){
        std::cout << "PHASE 1B - Normalizing the random field using sigma_8" << std::endl;
        std::cout << "1a) Planck sigma_8 = " << sigma_8 << "." << std::endl;
        std::cout << "1b) Integrated sigma_8 = " << integrated_sigma_8 << " from unnormalized power spectrum." << std::endl;
        std::cout << "1ba) Hubble ratio: " << H_hubble_of_z(z_start)/H_hubble_of_z(0) << "." << std::endl;
        std::cout << "   See " << std::string(OUTPUT_DIR) << "sigma_8_integration.txt." << std::endl;
        std::cout << "2a) Growth factor at z=" << z_start << " is " << D_growth_factor(z_start) << "." << std::endl;
        std::cout << "2b) Growth factor at z=" << 0 << " is " << D_growth_factor(0) << "." << std::endl;

        //Normalize the Gaussian random field by multiplying the Fourier modes by the appropriate factor
        global_PS_normalization = (sigma_8 / integrated_sigma_8) * (D_growth_factor(z_start) / D_growth_factor(0));
    } else {
        std::cout << "No valid normalization method specified. Exiting" << std::endl;
        return 0;
    }


    for (int x=0; x<N; x++) {
        for (int y=0; y<N; y++) {
            for (int z=0; z<=N/2; z++) {
                k_box[half_box_idx(N, x, y, z)][0] *= global_PS_normalization;
                k_box[half_box_idx(N, x, y, z)][1] *= global_PS_normalization;
            }
        }
    }

    std::cout << "=> The overall normalization is " << global_PS_normalization << "." << std::endl;
    std::cout << std::endl;
    std::cout << "Fluctuation at smallest k:" << std::endl;
    std::cout << "1) k_min = " << 2*M_PI/box_len << " Mpc^-1." << std::endl;
    std::cout << "2) P(k) = " << pow(sigma_func_cdm(2*M_PI/box_len)*global_PS_normalization,2) << " Mpc^3." << std::endl;
    std::cout << std::endl;

    std::cout << "PHASE 1C - Fourier transform of the random field to real field" << std::endl;

    //FTT
    fftw_plan c2r_plan  = fftw_plan_dft_c2r_3d(N, N, N, k_box, primordial_box, FFTW_ESTIMATE);
    fftw_execute(c2r_plan);

    //Normalization (from Fourier conventions alone)
    for (int x=0; x<N; x++) {
        for (int y=0; y<N; y++) {
            for (int z=0; z<N; z++) {
                primordial_box[box_idx(N, x, y, z)] /= box_volume;
            }
        }
    }

    //Write GRF to binary file
    write_array_to_disk(std::string(OUTPUT_DIR) + "gaussian.box", primordial_box, N);
    //Write GRF to H5 format
    writeGRF_H5(primordial_box, N, box_len, std::string(OUTPUT_DIR) + "gaussian.hdf5");

    std::cout << "1) The result has been written to " << std::string(OUTPUT_DIR) << "gaussian.box" << std::endl;
    std::cout << "2) The result has been written to " << std::string(OUTPUT_DIR) << "gaussian.hdf5" << std::endl;

    //Find the maximum and minimum values
    double box_min = primordial_box[0];
    double box_max = primordial_box[0];
    for (int i=1; i<N*N*N; i++) {
        box_min = (primordial_box[i] < box_min ? primordial_box[i] : box_min);
        box_max = (primordial_box[i] > box_max ? primordial_box[i] : box_max);
    }

    std::cout << "3) Random field range: " << box_min << " <= x <= " << box_max << "." << std::endl;
    std::cout << std::endl;

    //FFT the primordial box for later
    fftw_execute(r2c_plan);

    //Normalization
    for (int x=0; x<N; x++) {
        for (int y=0; y<N; y++) {
            for (int z=0; z<=N/2; z++) {
                k_box[half_box_idx(N, x, y, z)][0] *= box_volume/(N*N*N);
                k_box[half_box_idx(N, x, y, z)][1] *= box_volume/(N*N*N);
            }
        }
    }


    //Next, we either load particle positions (e.g. from a glass) or
    // generate them from a grid
    bool gridgen = true;
    const long int particle_num = PARTICLE_NUM;
    std::vector<corpuscle> bodies(particle_num);
    if (gridgen) {
        std::cout << "PHASE 2A - Placing cold particles on a grid" << std::endl;
        for (int x=0; x<NP; x++) {
            for (int y=0; y<NP; y++) {
                for (int z=0; z<NP; z++) {
                    corpuscle body;

                    body.id = (long int) box_idx(NP, x, y, z);
                    body.X = x*(box_len/NP);
                    body.Y = y*(box_len/NP);
                    body.Z = z*(box_len/NP);

                    bodies[body.id] = body;
                }
            }
        }
        std::cout << "1) Done with placing " << NP << "^3 cold particles." << std::endl << std::endl;
    }


    //Now determine the displacement field psi
    double *psi_x_box = (double*) fftw_malloc(sizeof(double)*N*N*N);
    double *psi_y_box = (double*) fftw_malloc(sizeof(double)*N*N*N);
    double *psi_z_box = (double*) fftw_malloc(sizeof(double)*N*N*N);
    fftw_complex *psi_k_x_box = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*N*(N/2+1));
    fftw_complex *psi_k_y_box = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*N*(N/2+1));
    fftw_complex *psi_k_z_box = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*N*(N/2+1));

    //If necessary, solve the Monge-Ampère equation for infinite LPT
    if (USE_INFINI_LPT) {
        double *phi_box = (double*) fftw_malloc(sizeof(double)*N*N*N);
        fftw_complex *phi_k_box = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*N*(N/2+1));

        std::cout << "PHASE 2B.i - Solve the Monge-Ampère equation for infinite LPT" << std::endl;
        double tol = 1e-15;
        do_infini_lpt(phi_box, primordial_box, N, box_len, tol);
        std::cout << std::endl;

        //Write GRF to binary file
        write_array_to_disk(std::string(OUTPUT_DIR) + "phi.box", phi_box, N);
        std::cout << "   Output written to " << std::string(OUTPUT_DIR) + "phi.box" << std::endl;
        std::cout << std::endl;

        std::cout << "PHASE 2B.ii - Compute the displacement vector field" << std::endl;
        std::cout << "1) Apply kernel to the Fourier transform of the potential field." << std::endl;

        //FFT the potential phi
        fftw_plan phi_plan = fftw_plan_dft_r2c_3d(N, N, N, phi_box, phi_k_box, FFTW_ESTIMATE);
    	fftw_execute(phi_plan);

        //Normalization
        for (int x=0; x<N; x++) {
            for (int y=0; y<N; y++) {
                for (int z=0; z<=N/2; z++) {
                    phi_k_box[half_box_idx(N, x, y, z)][0] *= box_volume/(N*N*N);
                    phi_k_box[half_box_idx(N, x, y, z)][1] *= box_volume/(N*N*N);
                }
            }
        }

        //Now differentiate to find psi
        for (int x=0; x<N; x++) {
            for (int y=0; y<N; y++) {
                for (int z=0; z<=N/2; z++) { //note that we stop at the (N/2+1)th entry
                    double k_x = (x > N/2) ? (x - N)*delta_k : x*delta_k; //Mpc^-1
                    double k_y = (y > N/2) ? (y - N)*delta_k : y*delta_k; //Mpc^-1
                    double k_z = (z > N/2) ? (z - N)*delta_k : z*delta_k; //Mpc^-1

                    double k = sqrt(k_x*k_x + k_y*k_y + k_z*k_z);

                    if (k>0) {
                        psi_k_x_box[half_box_idx(N, x, y, z)][0] = phi_k_box[half_box_idx(N, x, y, z)][1] * k_x;
                        psi_k_x_box[half_box_idx(N, x, y, z)][1] = -phi_k_box[half_box_idx(N, x, y, z)][0] * k_x;
                        psi_k_y_box[half_box_idx(N, x, y, z)][0] = phi_k_box[half_box_idx(N, x, y, z)][1] * k_y;
                        psi_k_y_box[half_box_idx(N, x, y, z)][1] = -phi_k_box[half_box_idx(N, x, y, z)][0] * k_y;
                        psi_k_z_box[half_box_idx(N, x, y, z)][0] = phi_k_box[half_box_idx(N, x, y, z)][1] * k_z;
                        psi_k_z_box[half_box_idx(N, x, y, z)][1] = -phi_k_box[half_box_idx(N, x, y, z)][0] * k_z;
                    } else {
                        psi_k_x_box[half_box_idx(N, x, y, z)][0] = 0;
                        psi_k_x_box[half_box_idx(N, x, y, z)][1] = 0;
                        psi_k_y_box[half_box_idx(N, x, y, z)][0] = 0;
                        psi_k_y_box[half_box_idx(N, x, y, z)][1] = 0;
                        psi_k_z_box[half_box_idx(N, x, y, z)][0] = 0;
                        psi_k_z_box[half_box_idx(N, x, y, z)][1] = 0;
                    }
                }
            }
        }
    } else {
        //Compute the displacement field from the random field
        std::cout << "PHASE 2B - Compute the displacement vector field" << std::endl;
        std::cout << "1) Apply kernel to the Fourier transform of the primordial field." << std::endl;

        //Multiply with the inverse Poisson kernel and differentiate
        for (int x=0; x<N; x++) {
            for (int y=0; y<N; y++) {
                for (int z=0; z<=N/2; z++) { //note that we stop at the (N/2+1)th entry
                    double k_x = (x > N/2) ? (x - N)*delta_k : x*delta_k; //Mpc^-1
                    double k_y = (y > N/2) ? (y - N)*delta_k : y*delta_k; //Mpc^-1
                    double k_z = (z > N/2) ? (z - N)*delta_k : z*delta_k; //Mpc^-1

                    double k = sqrt(k_x*k_x + k_y*k_y + k_z*k_z);

                    if (k>0) {
                        psi_k_x_box[half_box_idx(N, x, y, z)][0] = k_box[half_box_idx(N, x, y, z)][1] * k_x / (k*k);
                        psi_k_x_box[half_box_idx(N, x, y, z)][1] = -k_box[half_box_idx(N, x, y, z)][0] * k_x / (k*k);
                        psi_k_y_box[half_box_idx(N, x, y, z)][0] = k_box[half_box_idx(N, x, y, z)][1] * k_y / (k*k);
                        psi_k_y_box[half_box_idx(N, x, y, z)][1] = -k_box[half_box_idx(N, x, y, z)][0] * k_y / (k*k);
                        psi_k_z_box[half_box_idx(N, x, y, z)][0] = k_box[half_box_idx(N, x, y, z)][1] * k_z / (k*k);
                        psi_k_z_box[half_box_idx(N, x, y, z)][1] = -k_box[half_box_idx(N, x, y, z)][0] * k_z / (k*k);
                    } else {
                        psi_k_x_box[half_box_idx(N, x, y, z)][0] = 0;
                        psi_k_x_box[half_box_idx(N, x, y, z)][1] = 0;
                        psi_k_y_box[half_box_idx(N, x, y, z)][0] = 0;
                        psi_k_y_box[half_box_idx(N, x, y, z)][1] = 0;
                        psi_k_z_box[half_box_idx(N, x, y, z)][0] = 0;
                        psi_k_z_box[half_box_idx(N, x, y, z)][1] = 0;
                    }
                }
            }
        }
    }

    std::cout << "2) Fourier transform back to real coordinates." << std::endl;

	//Do the IFFTs
	fftw_plan px  = fftw_plan_dft_c2r_3d(N, N, N, psi_k_x_box, psi_x_box, FFTW_ESTIMATE);
	fftw_execute(px);

	fftw_plan py  = fftw_plan_dft_c2r_3d(N, N, N, psi_k_y_box, psi_y_box, FFTW_ESTIMATE);
	fftw_execute(py);

	fftw_plan pz  = fftw_plan_dft_c2r_3d(N, N, N, psi_k_z_box, psi_z_box, FFTW_ESTIMATE);
	fftw_execute(pz);

    //Normalization
	for (int x=0; x<N; x++) {
		for (int y=0; y<N; y++) {
			for (int z=0; z<N; z++) {
				psi_x_box[box_idx(N, x, y, z)] /= box_volume;
				psi_y_box[box_idx(N, x, y, z)] /= box_volume;
				psi_z_box[box_idx(N, x, y, z)] /= box_volume;
			}
		}
	}

    std::cout << "3) Done. The result consists of 3x" << N << "^3 real numbers." << std::endl;

    write_array_to_disk(std::string(OUTPUT_DIR) + "psi_x.box", psi_x_box, N);
	write_array_to_disk(std::string(OUTPUT_DIR) + "psi_y.box", psi_y_box, N);
	write_array_to_disk(std::string(OUTPUT_DIR) + "psi_z.box", psi_z_box, N);

    std::cout << "   Output written to " << std::string(OUTPUT_DIR) + "psi_x.box" << std::endl;
    std::cout << "   Output written to " << std::string(OUTPUT_DIR) + "psi_y.box" << std::endl;
    std::cout << "   Output written to " << std::string(OUTPUT_DIR) + "psi_z.box" << std::endl;
    std::cout << "" << std::endl;


    std::cout << "PHASE 2C - Displace the cold particles" << std::endl;

    for (auto body : bodies) {
        double X = body.X*N/box_len;
        double Y = body.Y*N/box_len;
        double Z = body.Z*N/box_len;

        //Grid positions
        int iX = (int) floor(X);
        int iY = (int) floor(Y);
        int iZ = (int) floor(Z);

        //Intepolate the necessary fields with TSC
        float lookLength = 1.0;
        int lookLftX = (int) floor((X-iX) - lookLength);
        int lookRgtX = (int) floor((X-iX) + lookLength);
        int lookLftY = (int) floor((Y-iY) - lookLength);
        int lookRgtY = (int) floor((Y-iY) + lookLength);
        int lookLftZ = (int) floor((Z-iZ) - lookLength);
        int lookRgtZ = (int) floor((Z-iZ) + lookLength);

        //Accumulate interpolated values in psi_{xyz}
        double psi_x = 0, psi_y = 0, psi_z = 0;

        for (int i=lookLftX; i<=lookRgtX; i++) {
            for (int j=lookLftY; j<=lookRgtY; j++) {
                for (int k=lookLftZ; k<=lookRgtZ; k++) {
                    //Pull the interpolated long-range force from the mesh
                    double xx = abs(X - (iX+i));
                    double yy = abs(Y - (iY+j));
                    double zz = abs(Z - (iZ+k));

                    double part_x = xx <= 1 ? 1-xx : 0;
                    double part_y = yy <= 1 ? 1-yy : 0;
                    double part_z = zz <= 1 ? 1-zz : 0;

                    psi_x += psi_x_box[box_wrap_idx(N, iX+i, iY+j, iZ+k)] * (part_x*part_y*part_z);
                    psi_y += psi_y_box[box_wrap_idx(N, iX+i, iY+j, iZ+k)] * (part_x*part_y*part_z);
                    psi_z += psi_z_box[box_wrap_idx(N, iX+i, iY+j, iZ+k)] * (part_x*part_y*part_z);
                }
            }
        }

        body.delta_X = -psi_x;
        body.delta_Y = -psi_y;
        body.delta_Z = -psi_z;

        //Displace the particles before assigning the velocities
        if (VELOCITY_AT_DISPLACED_POS) {
            body.X += body.delta_X;
            body.Y += body.delta_Y;
            body.Z += body.delta_Z;
        }

        body.mass = Omega_m*rho_crit*pow(Mpc,3)*box_volume/particle_num;

        bodies[body.id] = body;
    }

    if (VELOCITY_AT_DISPLACED_POS) {
        std::cout << "1) Displaced " << particle_num << " cold particles." << std::endl;
    } else {
        std::cout << "1) Determined " << particle_num << " displacement vectors, particles will be moved later." << std::endl;
    }
    std::cout << "2) Particle mass " << bodies[0].mass << " kg." << std::endl;
    std::cout << "  " << std::endl;


    std::cout << "PHASE 2D - Undo cdm transfer function & export pure Gaussian field" << std::endl;

    //Undo the CDM density transfer function and obtain a pure Gaussian field (no transfer function)
    for (int x=0; x<N; x++) {
        for (int y=0; y<N; y++) {
            for (int z=0; z<=N/2; z++) { //note that we stop at the (N/2+1)th entry
                double k_x = (x > N/2) ? (x - N)*delta_k : x*delta_k; //Mpc^-1
                double k_y = (y > N/2) ? (y - N)*delta_k : y*delta_k; //Mpc^-1
                double k_z = (z > N/2) ? (z - N)*delta_k : z*delta_k; //Mpc^-1

                double k = sqrt(k_x*k_x + k_y*k_y + k_z*k_z);

                if (k > 0) {
                    k_box[half_box_idx(N, x, y, z)][0] *= sigma_func_pure(k)/sigma_func_cdm(k);
                    k_box[half_box_idx(N, x, y, z)][1] *= sigma_func_pure(k)/sigma_func_cdm(k);
                }
            }
        }
    }

    //FTT back
    fftw_execute(c2r_plan);

    //Normalization (from Fourier conventions alone)
    for (int x=0; x<N; x++) {
        for (int y=0; y<N; y++) {
            for (int z=0; z<N; z++) {
                primordial_box[box_idx(N, x, y, z)] /= box_volume;
            }
        }
    }

    //Write GRF to binary file
    write_array_to_disk(std::string(OUTPUT_DIR) + "gaussian_pure.box", primordial_box, N);
    //Write GRF to H5 format
    writeGRF_H5(primordial_box, N, box_len, std::string(OUTPUT_DIR) + "gaussian_pure.hdf5");

    std::cout << "1) The result has been written to " << std::string(OUTPUT_DIR) << "gaussian_pure.box" << std::endl;
    std::cout << "2) The result has been written to " << std::string(OUTPUT_DIR) << "gaussian_pure.hdf5" << std::endl;
    std::cout << "  " << std::endl;

    //FFT the primordial box for later
    fftw_execute(r2c_plan);

    //Normalization
    for (int x=0; x<N; x++) {
        for (int y=0; y<N; y++) {
            for (int z=0; z<=N/2; z++) {
                k_box[half_box_idx(N, x, y, z)][0] *= box_volume/(N*N*N);
                k_box[half_box_idx(N, x, y, z)][1] *= box_volume/(N*N*N);
            }
        }
    }

    double dVdX;
    if (VELOCITY_METHOD == VEL_ZELDOVICH) {
        std::cout << "PHASE 2E - Compute the velocity proportionality constant for the Zel'dovich method" << std::endl;

        //Compute the constant of proportionality
        double a = a_scale_factor_of_z(z_start);
        double H = H_hubble_of_z(z_start);
        double f = logarithmic_derivative_f_1(z_start);

        /* N.B. The reason we have a^2Hf, is because HeWon and Swift (although
         * Swift wants the peculiar velocity v=a(dx/dt) in the ic file) both use
         * generalized velocity coordinates p = a^2 (dx/dt) where x is the
         * comoving position. This is also the convention of Quinn et al. (1997)
         * astro-ph/9710043 where kick and drift operators are derived.
         */
        dVdX = pow(a,2)*H*f;

        std::cout << "1) Velocity convention v = a^2 (dx/dt), where x=r/a is comoving." << std::endl;
        std::cout << "2) Expansion factor a(z) = " << a << "." << std::endl;
        std::cout << "3) Hubble rate H(z) = " << H << " Gyr^-1." << std::endl;
        std::cout << "4) Logarithmic derivative of growth factor f(z) = " << f << "." << std::endl;
        std::cout << "5) Proportionality constant dVdX = a^2Hf = " << dVdX << " Gyr^-1." << std::endl;
        std::cout << "  " << std::endl;

        std::cout << "PHASE 2F - Assign initial velocities to the cold particles using the Zel'dovich method" << std::endl;

        for (auto body : bodies) {
            double X = body.X*N/box_len;
            double Y = body.Y*N/box_len;
            double Z = body.Z*N/box_len;

            //Grid positions
            int iX = (int) floor(X);
            int iY = (int) floor(Y);
            int iZ = (int) floor(Z);

            //Intepolate the necessary fields with TSC
            float lookLength = 1.0;
            int lookLftX = (int) floor((X-iX) - lookLength);
            int lookRgtX = (int) floor((X-iX) + lookLength);
            int lookLftY = (int) floor((Y-iY) - lookLength);
            int lookRgtY = (int) floor((Y-iY) + lookLength);
            int lookLftZ = (int) floor((Z-iZ) - lookLength);
            int lookRgtZ = (int) floor((Z-iZ) + lookLength);

            //Accumulate interpolated values in psi_{xyz}
            double psi_x = 0, psi_y = 0, psi_z = 0;

            for (int i=lookLftX; i<=lookRgtX; i++) {
                for (int j=lookLftY; j<=lookRgtY; j++) {
                    for (int k=lookLftZ; k<=lookRgtZ; k++) {
                        //Pull the interpolated long-range force from the mesh
                        double xx = abs(X - (iX+i));
                        double yy = abs(Y - (iY+j));
                        double zz = abs(Z - (iZ+k));

                        double part_x = xx <= 1 ? 1-xx : 0;
                        double part_y = yy <= 1 ? 1-yy : 0;
                        double part_z = zz <= 1 ? 1-zz : 0;

                        psi_x += psi_x_box[box_wrap_idx(N, iX+i, iY+j, iZ+k)] * (part_x*part_y*part_z);
                        psi_y += psi_y_box[box_wrap_idx(N, iX+i, iY+j, iZ+k)] * (part_x*part_y*part_z);
                        psi_z += psi_z_box[box_wrap_idx(N, iX+i, iY+j, iZ+k)] * (part_x*part_y*part_z);
                    }
                }
            }

            body.v_X = -dVdX * psi_x;
            body.v_Y = -dVdX * psi_y;
            body.v_Z = -dVdX * psi_z;

            //Displace the particles if this has not yet been done
            if (!VELOCITY_AT_DISPLACED_POS) {
                body.X += body.delta_X;
                body.Y += body.delta_Y;
                body.Z += body.delta_Z;
            }

            bodies[body.id] = body;
        }

        std::cout << "1) Assigned velocities to " << particle_num << " cold particles." << std::endl;
        if (!VELOCITY_AT_DISPLACED_POS) {
            std::cout << "2) Displaced " << particle_num << " cold particles." << std::endl;
        }

    } else if (VELOCITY_METHOD == VEL_TRANSFER) {
        std::cout << "PHASE 2E - Applying the cold velocity transfer function" << std::endl;

        //Undo the CDM density transfer function and apply the CDM velocity transfer function
        for (int x=0; x<N; x++) {
            for (int y=0; y<N; y++) {
                for (int z=0; z<=N/2; z++) { //note that we stop at the (N/2+1)th entry
                    double k_x = (x > N/2) ? (x - N)*delta_k : x*delta_k; //Mpc^-1
                    double k_y = (y > N/2) ? (y - N)*delta_k : y*delta_k; //Mpc^-1
                    double k_z = (z > N/2) ? (z - N)*delta_k : z*delta_k; //Mpc^-1

                    double k = sqrt(k_x*k_x + k_y*k_y + k_z*k_z);

                    if (k > 0) {
                        k_box[half_box_idx(N, x, y, z)][0] *= sigma_func_vel_cdm(k)/sigma_func_pure(k);
                        k_box[half_box_idx(N, x, y, z)][1] *= sigma_func_vel_cdm(k)/sigma_func_pure(k);
                    }
                }
            }
        }

        //FTT back
        fftw_execute(c2r_plan);

        //Normalization (from Fourier conventions alone)
        for (int x=0; x<N; x++) {
            for (int y=0; y<N; y++) {
                for (int z=0; z<N; z++) {
                    primordial_box[box_idx(N, x, y, z)] /= box_volume;
                }
            }
        }

        //Write GRF to binary file
        write_array_to_disk(std::string(OUTPUT_DIR) + "gaussian_theta.box", primordial_box, N);
        //Write GRF to H5 format
        writeGRF_H5(primordial_box, N, box_len, std::string(OUTPUT_DIR) + "gaussian_theta.hdf5");

        std::cout << "1) The result has been written to " << std::string(OUTPUT_DIR) << "gaussian_theta.box" << std::endl;
        std::cout << "2) The result has been written to " << std::string(OUTPUT_DIR) << "gaussian_theta.hdf5" << std::endl;
        std::cout << "3) Theta is the divergence of the peculiar velocity in Mpc^-1," << std::endl;
        std::cout << "   which is (dx/dtau) = a*(dx/dt) with x=r/a comoving." << std::endl;
        std::cout << std::endl;

        //Compute the velocity field from the random field
        std::cout << "PHASE 2F - Compute the cold velocity field" << std::endl;
        std::cout << "1) Velocity convention v = a^2 (dx/dt), where x=r/a is comoving." << std::endl;
        std::cout << "2) Apply kernel to the Fourier transform of the divergence field." << std::endl;

        //FTT the primordial box
        fftw_execute(r2c_plan);

        //Normalization
        for (int x=0; x<N; x++) {
            for (int y=0; y<N; y++) {
                for (int z=0; z<=N/2; z++) {
                    k_box[half_box_idx(N, x, y, z)][0] *= box_volume/(N*N*N);
                    k_box[half_box_idx(N, x, y, z)][1] *= box_volume/(N*N*N);
                }
            }
        }

        for (int x=0; x<N; x++) {
            for (int y=0; y<N; y++) {
                for (int z=0; z<=N/2; z++) { //note that we stop at the (N/2+1)th entry
                    double k_x = (x > N/2) ? (x - N)*delta_k : x*delta_k; //Mpc^-1
                    double k_y = (y > N/2) ? (y - N)*delta_k : y*delta_k; //Mpc^-1
                    double k_z = (z > N/2) ? (z - N)*delta_k : z*delta_k; //Mpc^-1

                    double k = sqrt(k_x*k_x + k_y*k_y + k_z*k_z);

                    if (k>0) {
                        psi_k_x_box[half_box_idx(N, x, y, z)][0] = k_box[half_box_idx(N, x, y, z)][1] * k_x / (k*k);
                        psi_k_x_box[half_box_idx(N, x, y, z)][1] = -k_box[half_box_idx(N, x, y, z)][0] * k_x / (k*k);
                        psi_k_y_box[half_box_idx(N, x, y, z)][0] = k_box[half_box_idx(N, x, y, z)][1] * k_y / (k*k);
                        psi_k_y_box[half_box_idx(N, x, y, z)][1] = -k_box[half_box_idx(N, x, y, z)][0] * k_y / (k*k);
                        psi_k_z_box[half_box_idx(N, x, y, z)][0] = k_box[half_box_idx(N, x, y, z)][1] * k_z / (k*k);
                        psi_k_z_box[half_box_idx(N, x, y, z)][1] = -k_box[half_box_idx(N, x, y, z)][0] * k_z / (k*k);
                    } else {
                        psi_k_x_box[half_box_idx(N, x, y, z)][0] = 0;
                        psi_k_x_box[half_box_idx(N, x, y, z)][1] = 0;
                        psi_k_y_box[half_box_idx(N, x, y, z)][0] = 0;
                        psi_k_y_box[half_box_idx(N, x, y, z)][1] = 0;
                        psi_k_z_box[half_box_idx(N, x, y, z)][0] = 0;
                        psi_k_z_box[half_box_idx(N, x, y, z)][1] = 0;
                    }
                }
            }
        }

        std::cout << "3) Fourier transform back to real coordinates." << std::endl;

        //Do the IFFTs
        fftw_execute(px);
        fftw_execute(py);
        fftw_execute(pz);

        //Normalization
        for (int x=0; x<N; x++) {
            for (int y=0; y<N; y++) {
                for (int z=0; z<N; z++) {
                    psi_x_box[box_idx(N, x, y, z)] /= box_volume;
                    psi_y_box[box_idx(N, x, y, z)] /= box_volume;
                    psi_z_box[box_idx(N, x, y, z)] /= box_volume;
                }
            }
        }


        std::cout << "4) Done. The result consists of 3x" << N << "^3 real numbers." << std::endl;
        std::cout << "5) Make dimensionful by inserting c = " << c_vel << " Mpc/Gyr." << std::endl;

        //Insert units (convert from c=1 to dimensionful quantity)
        for (int x=0; x<N; x++) {
            for (int y=0; y<N; y++) {
                for (int z=0; z<N; z++) {
                    psi_x_box[box_idx(N, x, y, z)] *= c_vel;
                    psi_y_box[box_idx(N, x, y, z)] *= c_vel;
                    psi_z_box[box_idx(N, x, y, z)] *= c_vel;
                }
            }
        }

        double a = a_scale_factor_of_z(z_start);
        std::cout << "6) Multiply by the scale factor a = " << a << " to get v = a*(dx/tau) = a^2*(dx/dt)." << std::endl;

        //Insert scale factor
        for (int x=0; x<N; x++) {
            for (int y=0; y<N; y++) {
                for (int z=0; z<N; z++) {
                    psi_x_box[box_idx(N, x, y, z)] *= a;
                    psi_y_box[box_idx(N, x, y, z)] *= a;
                    psi_z_box[box_idx(N, x, y, z)] *= a;
                }
            }
        }

        write_array_to_disk(std::string(OUTPUT_DIR) + "v_cdm_x.box", psi_x_box, N);
        write_array_to_disk(std::string(OUTPUT_DIR) + "v_cdm_y.box", psi_y_box, N);
        write_array_to_disk(std::string(OUTPUT_DIR) + "v_cdm_z.box", psi_z_box, N);

        std::cout << "   Output written to " << std::string(OUTPUT_DIR) + "v_cdm_x.box" << std::endl;
        std::cout << "   Output written to " << std::string(OUTPUT_DIR) + "v_cdm_y.box" << std::endl;
        std::cout << "   Output written to " << std::string(OUTPUT_DIR) + "v_cdm_z.box" << std::endl;
        std::cout << "" << std::endl;


        std::cout << "PHASE 2G - Assign initial velocities to the cold particles" << std::endl;

        for (auto body : bodies) {
            double X = body.X*N/box_len;
            double Y = body.Y*N/box_len;
            double Z = body.Z*N/box_len;

            //Grid positions
            int iX = (int) floor(X);
            int iY = (int) floor(Y);
            int iZ = (int) floor(Z);

            //Intepolate the necessary fields with TSC
            float lookLength = 1.0;
            int lookLftX = (int) floor((X-iX) - lookLength);
            int lookRgtX = (int) floor((X-iX) + lookLength);
            int lookLftY = (int) floor((Y-iY) - lookLength);
            int lookRgtY = (int) floor((Y-iY) + lookLength);
            int lookLftZ = (int) floor((Z-iZ) - lookLength);
            int lookRgtZ = (int) floor((Z-iZ) + lookLength);

            //Accumulate interpolated values in psi_{xyz}
            double psi_x = 0, psi_y = 0, psi_z = 0;

            for (int i=lookLftX; i<=lookRgtX; i++) {
                for (int j=lookLftY; j<=lookRgtY; j++) {
                    for (int k=lookLftZ; k<=lookRgtZ; k++) {
                        //Pull the interpolated long-range force from the mesh
                        double xx = abs(X - (iX+i));
                        double yy = abs(Y - (iY+j));
                        double zz = abs(Z - (iZ+k));

                        double part_x = xx <= 1 ? 1-xx : 0;
                        double part_y = yy <= 1 ? 1-yy : 0;
                        double part_z = zz <= 1 ? 1-zz : 0;

                        psi_x += psi_x_box[box_wrap_idx(N, iX+i, iY+j, iZ+k)] * (part_x*part_y*part_z);
                        psi_y += psi_y_box[box_wrap_idx(N, iX+i, iY+j, iZ+k)] * (part_x*part_y*part_z);
                        psi_z += psi_z_box[box_wrap_idx(N, iX+i, iY+j, iZ+k)] * (part_x*part_y*part_z);
                    }
                }
            }

            body.v_X = psi_x;
            body.v_Y = psi_y;
            body.v_Z = psi_z;

            //Displace the particles if this has not yet been done
            if (!VELOCITY_AT_DISPLACED_POS) {
                body.X += body.delta_X;
                body.Y += body.delta_Y;
                body.Z += body.delta_Z;
            }

            bodies[body.id] = body;
        }

        std::cout << "1) Assigned velocities to " << particle_num << " cold particles." << std::endl;
        if (!VELOCITY_AT_DISPLACED_POS) {
            std::cout << "2) Displaced " << particle_num << " cold particles." << std::endl;
        }
    } else {
        std::cout << "No valid method to determine the initial velocities. Exiting" << std::endl;
        return 0;
    }


    //Determine the average & maximum cold velocity for diagnostic purposes
    double avg_cb_speed = 0;
    double max_cb_speed = 0;
    double V = 0;

    for (auto body : bodies) {
        double V = sqrt(body.v_X*body.v_X + body.v_Y*body.v_Y + body.v_Z*body.v_Z);
        avg_cb_speed += V/particle_num;
        if (V > max_cb_speed) {
            max_cb_speed = V;
        }
    }

    std::cout << "  " << std::endl;
    std::cout << "PHASE 2G - Summary" << std::endl;
    std::cout << "1) Average speed " << avg_cb_speed << " comoving Mpc/Gyr." << std::endl;
    std::cout << "2) Maximum speed " << max_cb_speed << " comoving Mpc/Gyr." << std::endl;
    std::cout << "  " << std::endl;

    std::cout << "  ...  " << std::endl;
    std::cout << "  " << std::endl;

    /* NEXT, the neutrinos */

    //Next, we either load particle positions (e.g. from a glass) or
    // generate them from a grid
    bool gridgen_nu = true;
    const long int neutrino_num = NEUTRINO_NUM;
    std::vector<corpuscle> bodies_nu(neutrino_num);
    if (gridgen_nu) {
        std::cout << "PHASE 3A - Placing neutrino particles on a grid" << std::endl;
        for (int x=0; x<NNUP; x++) {
            for (int y=0; y<NNUP; y++) {
                for (int z=0; z<NNUP; z++) {
                    corpuscle body;

                    body.id = (long int) box_idx(NNUP, x, y, z);
                    body.X = (x+0.5)*(box_len/NNUP);
                    body.Y = (y+0.5)*(box_len/NNUP);
                    body.Z = (z+0.5)*(box_len/NNUP);

                    bodies_nu[body.id] = body;
                }
            }
        }
        std::cout << "1) Done with placing " << NNUP << "^3 neutrino particles." << std::endl << std::endl;
    }

    std::cout << "PHASE 3B - Applying the neutrino transfer function" << std::endl;

    //Undo the CDM transfer function and apply the neutrino transfer function
    for (int x=0; x<N; x++) {
        for (int y=0; y<N; y++) {
            for (int z=0; z<=N/2; z++) { //note that we stop at the (N/2+1)th entry
                double k_x = (x > N/2) ? (x - N)*delta_k : x*delta_k; //Mpc^-1
                double k_y = (y > N/2) ? (y - N)*delta_k : y*delta_k; //Mpc^-1
                double k_z = (z > N/2) ? (z - N)*delta_k : z*delta_k; //Mpc^-1

                double k = sqrt(k_x*k_x + k_y*k_y + k_z*k_z);

                if (k > 0) {
                    if (VELOCITY_METHOD == VEL_ZELDOVICH) {
                        k_box[half_box_idx(N, x, y, z)][0] *= sigma_func_neutrino(k)/sigma_func_pure(k);
                        k_box[half_box_idx(N, x, y, z)][1] *= sigma_func_neutrino(k)/sigma_func_pure(k);
                    } else if (VELOCITY_METHOD == VEL_TRANSFER) {
                        k_box[half_box_idx(N, x, y, z)][0] *= sigma_func_neutrino(k)/sigma_func_vel_cdm(k);
                        k_box[half_box_idx(N, x, y, z)][1] *= sigma_func_neutrino(k)/sigma_func_vel_cdm(k);
                    } else {
                        std::cout << "No valid method to determine the initial velocities. Exiting" << std::endl;
                        return 0;
                    }
                }
            }
        }
    }

    //FTT back
    fftw_execute(c2r_plan);

    //Normalization (from Fourier conventions alone)
    for (int x=0; x<N; x++) {
        for (int y=0; y<N; y++) {
            for (int z=0; z<N; z++) {
                primordial_box[box_idx(N, x, y, z)] /= box_volume;
            }
        }
    }

    //Write GRF to binary file
    write_array_to_disk(std::string(OUTPUT_DIR) + "gaussian_nu.box", primordial_box, N);
    //Write GRF to H5 format
    writeGRF_H5(primordial_box, N, box_len, std::string(OUTPUT_DIR) + "gaussian_nu.hdf5");

    std::cout << "1) The result has been written to " << std::string(OUTPUT_DIR) << "gaussian_nu.box" << std::endl;
    std::cout << "2) The result has been written to " << std::string(OUTPUT_DIR) << "gaussian_nu.hdf5" << std::endl;
    std::cout << std::endl;

    //FTT the primordial box with (neutrino transfer function) for later
    fftw_execute(r2c_plan);

    //Normalization
    for (int x=0; x<N; x++) {
        for (int y=0; y<N; y++) {
            for (int z=0; z<=N/2; z++) {
                k_box[half_box_idx(N, x, y, z)][0] *= box_volume/(N*N*N);
                k_box[half_box_idx(N, x, y, z)][1] *= box_volume/(N*N*N);
            }
        }
    }

    //Compute the neutrino displacement field from the random field

    //If necessary, solve the Monge-Ampère equation for infinite LPT
    if (USE_INFINI_LPT) {
        double *phi_box = (double*) fftw_malloc(sizeof(double)*N*N*N);
        fftw_complex *phi_k_box = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*N*(N/2+1));

        std::cout << "PHASE 3C.i - Solve the Monge-Ampère equation for infinite LPT" << std::endl;
        double tol = 1e-12;
        do_infini_lpt(phi_box, primordial_box, N, box_len, tol);
        std::cout << std::endl;

        //Write GRF to binary file
        write_array_to_disk(std::string(OUTPUT_DIR) + "phi_nu.box", phi_box, N);
        std::cout << "   Output written to " << std::string(OUTPUT_DIR) + "phi_nu.box" << std::endl;
        std::cout << std::endl;

        std::cout << "PHASE 3C.ii - Compute the neutrino displacement vector field" << std::endl;
        std::cout << "1) Apply kernel to the Fourier transform of the potential field." << std::endl;

        //FFT the potential phi
        fftw_plan phi_plan = fftw_plan_dft_r2c_3d(N, N, N, phi_box, phi_k_box, FFTW_ESTIMATE);
        fftw_execute(phi_plan);

        //Normalization
        for (int x=0; x<N; x++) {
            for (int y=0; y<N; y++) {
                for (int z=0; z<=N/2; z++) {
                    phi_k_box[half_box_idx(N, x, y, z)][0] *= box_volume/(N*N*N);
                    phi_k_box[half_box_idx(N, x, y, z)][1] *= box_volume/(N*N*N);
                }
            }
        }

        for (int x=0; x<N; x++) {
            for (int y=0; y<N; y++) {
                for (int z=0; z<=N/2; z++) { //note that we stop at the (N/2+1)th entry
                    double k_x = (x > N/2) ? (x - N)*delta_k : x*delta_k; //Mpc^-1
                    double k_y = (y > N/2) ? (y - N)*delta_k : y*delta_k; //Mpc^-1
                    double k_z = (z > N/2) ? (z - N)*delta_k : z*delta_k; //Mpc^-1

                    double k = sqrt(k_x*k_x + k_y*k_y + k_z*k_z);

                    if (k>0) {
                        psi_k_x_box[half_box_idx(N, x, y, z)][0] = phi_k_box[half_box_idx(N, x, y, z)][1] * k_x;
                        psi_k_x_box[half_box_idx(N, x, y, z)][1] = -phi_k_box[half_box_idx(N, x, y, z)][0] * k_x;
                        psi_k_y_box[half_box_idx(N, x, y, z)][0] = phi_k_box[half_box_idx(N, x, y, z)][1] * k_y;
                        psi_k_y_box[half_box_idx(N, x, y, z)][1] = -phi_k_box[half_box_idx(N, x, y, z)][0] * k_y;
                        psi_k_z_box[half_box_idx(N, x, y, z)][0] = phi_k_box[half_box_idx(N, x, y, z)][1] * k_z;
                        psi_k_z_box[half_box_idx(N, x, y, z)][1] = -phi_k_box[half_box_idx(N, x, y, z)][0] * k_z;
                    } else {
                        psi_k_x_box[half_box_idx(N, x, y, z)][0] = 0;
                        psi_k_x_box[half_box_idx(N, x, y, z)][1] = 0;
                        psi_k_y_box[half_box_idx(N, x, y, z)][0] = 0;
                        psi_k_y_box[half_box_idx(N, x, y, z)][1] = 0;
                        psi_k_z_box[half_box_idx(N, x, y, z)][0] = 0;
                        psi_k_z_box[half_box_idx(N, x, y, z)][1] = 0;
                    }
                }
            }
        }
    } else {
        std::cout << "PHASE 3C - Compute the neutrino displacement vector field" << std::endl;
        std::cout << "1) Apply kernel to the Fourier transform of the primordial field." << std::endl;

        for (int x=0; x<N; x++) {
            for (int y=0; y<N; y++) {
                for (int z=0; z<=N/2; z++) { //note that we stop at the (N/2+1)th entry
                    double k_x = (x > N/2) ? (x - N)*delta_k : x*delta_k; //Mpc^-1
                    double k_y = (y > N/2) ? (y - N)*delta_k : y*delta_k; //Mpc^-1
                    double k_z = (z > N/2) ? (z - N)*delta_k : z*delta_k; //Mpc^-1

                    double k = sqrt(k_x*k_x + k_y*k_y + k_z*k_z);

                    if (k>0) {
                        psi_k_x_box[half_box_idx(N, x, y, z)][0] = k_box[half_box_idx(N, x, y, z)][1] * k_x / (k*k);
                        psi_k_x_box[half_box_idx(N, x, y, z)][1] = -k_box[half_box_idx(N, x, y, z)][0] * k_x / (k*k);
                        psi_k_y_box[half_box_idx(N, x, y, z)][0] = k_box[half_box_idx(N, x, y, z)][1] * k_y / (k*k);
                        psi_k_y_box[half_box_idx(N, x, y, z)][1] = -k_box[half_box_idx(N, x, y, z)][0] * k_y / (k*k);
                        psi_k_z_box[half_box_idx(N, x, y, z)][0] = k_box[half_box_idx(N, x, y, z)][1] * k_z / (k*k);
                        psi_k_z_box[half_box_idx(N, x, y, z)][1] = -k_box[half_box_idx(N, x, y, z)][0] * k_z / (k*k);
                    } else {
                        psi_k_x_box[half_box_idx(N, x, y, z)][0] = 0;
                        psi_k_x_box[half_box_idx(N, x, y, z)][1] = 0;
                        psi_k_y_box[half_box_idx(N, x, y, z)][0] = 0;
                        psi_k_y_box[half_box_idx(N, x, y, z)][1] = 0;
                        psi_k_z_box[half_box_idx(N, x, y, z)][0] = 0;
                        psi_k_z_box[half_box_idx(N, x, y, z)][1] = 0;
                    }
                }
            }
        }
    }

    std::cout << "2) Fourier transform back to real coordinates." << std::endl;

	//Do the IFFTs
	fftw_execute(px);
	// fftw_destroy_plan(px);
    // fftw_free(psi_k_x_box);

	fftw_execute(py);
	// fftw_destroy_plan(py);
    // fftw_free(psi_k_y_box);

	fftw_execute(pz);
	// fftw_destroy_plan(pz);
    // fftw_free(psi_k_z_box);

    //Normalization
	for (int x=0; x<N; x++) {
		for (int y=0; y<N; y++) {
			for (int z=0; z<N; z++) {
				psi_x_box[box_idx(N, x, y, z)] /= box_volume;
				psi_y_box[box_idx(N, x, y, z)] /= box_volume;
				psi_z_box[box_idx(N, x, y, z)] /= box_volume;
			}
		}
	}

    std::cout << "3) Done. The result consists of 3x" << N << "^3 real numbers." << std::endl;

    write_array_to_disk(std::string(OUTPUT_DIR) + "psi_nu_x.box", psi_x_box, N);
	write_array_to_disk(std::string(OUTPUT_DIR) + "psi_nu_y.box", psi_y_box, N);
	write_array_to_disk(std::string(OUTPUT_DIR) + "psi_nu_z.box", psi_z_box, N);

    std::cout << "   Output written to " << std::string(OUTPUT_DIR) + "psi_nu_x.box" << std::endl;
    std::cout << "   Output written to " << std::string(OUTPUT_DIR) + "psi_nu_y.box" << std::endl;
    std::cout << "   Output written to " << std::string(OUTPUT_DIR) + "psi_nu_z.box" << std::endl;
    std::cout << "" << std::endl;

    std::cout << "PHASE 3D - Displace the neutrino particles" << std::endl;

    for (auto body : bodies_nu) {
		double X = body.X*N/box_len;
		double Y = body.Y*N/box_len;
		double Z = body.Z*N/box_len;

        //Grid positions
		int iX = (int) floor(X);
		int iY = (int) floor(Y);
		int iZ = (int) floor(Z);

		//Intepolate the necessary fields with TSC
		float lookLength = 1.0;
		int lookLftX = (int) floor((X-iX) - lookLength);
		int lookRgtX = (int) floor((X-iX) + lookLength);
		int lookLftY = (int) floor((Y-iY) - lookLength);
		int lookRgtY = (int) floor((Y-iY) + lookLength);
		int lookLftZ = (int) floor((Z-iZ) - lookLength);
		int lookRgtZ = (int) floor((Z-iZ) + lookLength);

        //Accumulate interpolated values in psi_{xyz}
		double psi_x = 0, psi_y = 0, psi_z = 0;

		for (int i=lookLftX; i<=lookRgtX; i++) {
			for (int j=lookLftY; j<=lookRgtY; j++) {
				for (int k=lookLftZ; k<=lookRgtZ; k++) {
					//Pull the interpolated long-range force from the mesh
					double xx = abs(X - (iX+i));
					double yy = abs(Y - (iY+j));
					double zz = abs(Z - (iZ+k));

					double part_x = xx <= 1 ? 1-xx : 0;
					double part_y = yy <= 1 ? 1-yy : 0;
					double part_z = zz <= 1 ? 1-zz : 0;

					psi_x += psi_x_box[box_wrap_idx(N, iX+i, iY+j, iZ+k)] * (part_x*part_y*part_z);
					psi_y += psi_y_box[box_wrap_idx(N, iX+i, iY+j, iZ+k)] * (part_x*part_y*part_z);
					psi_z += psi_z_box[box_wrap_idx(N, iX+i, iY+j, iZ+k)] * (part_x*part_y*part_z);
				}
			}
		}

        body.delta_X = -psi_x;
        body.delta_Y = -psi_y;
        body.delta_Z = -psi_z;

        //Displace the particles before assigning the velocities
        if (VELOCITY_AT_DISPLACED_POS) {
            body.X += body.delta_X;
            body.Y += body.delta_Y;
            body.Z += body.delta_Z;
        }

        //We are treating the three species as degenerate for now
        body.mass = Omega_nu*rho_crit*pow(Mpc,3)*box_volume/neutrino_num;

		bodies_nu[body.id] = body;
	}

    if (VELOCITY_AT_DISPLACED_POS) {
        std::cout << "1) Displaced " << neutrino_num << " neutrino particles." << std::endl;
    } else {
        std::cout << "1) Determined " << neutrino_num << " displacement vectors, particles will be moved later." << std::endl;
    }

    std::cout << "2) Particle mass " << bodies_nu[0].mass << " kg." << std::endl;
    std::cout << "  " << std::endl;

    if (VELOCITY_METHOD == VEL_ZELDOVICH) {
        std::cout << "PHASE 3E - Assign initial velocities to the neutrinos using the Zel'dovich method" << std::endl;
        std::cout << "1) Velocity convention v = a^2 (dx/dt), where x=r/a is comoving." << std::endl;

        for (auto body : bodies_nu) {
            double X = body.X*N/box_len;
            double Y = body.Y*N/box_len;
            double Z = body.Z*N/box_len;

            //Grid positions
            int iX = (int) floor(X);
            int iY = (int) floor(Y);
            int iZ = (int) floor(Z);

            //Intepolate the necessary fields with TSC
            float lookLength = 1.0;
            int lookLftX = (int) floor((X-iX) - lookLength);
            int lookRgtX = (int) floor((X-iX) + lookLength);
            int lookLftY = (int) floor((Y-iY) - lookLength);
            int lookRgtY = (int) floor((Y-iY) + lookLength);
            int lookLftZ = (int) floor((Z-iZ) - lookLength);
            int lookRgtZ = (int) floor((Z-iZ) + lookLength);

            //Accumulate interpolated values in psi_{xyz}
            double psi_x = 0, psi_y = 0, psi_z = 0;

            for (int i=lookLftX; i<=lookRgtX; i++) {
                for (int j=lookLftY; j<=lookRgtY; j++) {
                    for (int k=lookLftZ; k<=lookRgtZ; k++) {
                        //Pull the interpolated long-range force from the mesh
                        double xx = abs(X - (iX+i));
                        double yy = abs(Y - (iY+j));
                        double zz = abs(Z - (iZ+k));

                        double part_x = xx <= 1 ? 1-xx : 0;
                        double part_y = yy <= 1 ? 1-yy : 0;
                        double part_z = zz <= 1 ? 1-zz : 0;

                        psi_x += psi_x_box[box_wrap_idx(N, iX+i, iY+j, iZ+k)] * (part_x*part_y*part_z);
                        psi_y += psi_y_box[box_wrap_idx(N, iX+i, iY+j, iZ+k)] * (part_x*part_y*part_z);
                        psi_z += psi_z_box[box_wrap_idx(N, iX+i, iY+j, iZ+k)] * (part_x*part_y*part_z);
                    }
                }
            }

            body.v_X = -dVdX * psi_x;
            body.v_Y = -dVdX * psi_y;
            body.v_Z = -dVdX * psi_z;

            //Displace the particles after assigning velocities
            if (!VELOCITY_AT_DISPLACED_POS) {
                body.X += body.delta_X;
                body.Y += body.delta_Y;
                body.Z += body.delta_Z;
            }

            bodies_nu[body.id] = body;
        }

        std::cout << "2) Assigned velocities to " << neutrino_num << " neutrino particles." << std::endl;
        if (!VELOCITY_AT_DISPLACED_POS) {
            std::cout << "3) Displaced " << particle_num << " cold particles." << std::endl;
        }
        std::cout << "  " << std::endl;
    } else if (VELOCITY_METHOD == VEL_TRANSFER) {
        std::cout << "PHASE 3E - Applying the neutrino velocity transfer function" << std::endl;

        //Undo the neutrino density transfer function and apply the neutrino velocity transfer function
        for (int x=0; x<N; x++) {
            for (int y=0; y<N; y++) {
                for (int z=0; z<=N/2; z++) { //note that we stop at the (N/2+1)th entry
                    double k_x = (x > N/2) ? (x - N)*delta_k : x*delta_k; //Mpc^-1
                    double k_y = (y > N/2) ? (y - N)*delta_k : y*delta_k; //Mpc^-1
                    double k_z = (z > N/2) ? (z - N)*delta_k : z*delta_k; //Mpc^-1

                    double k = sqrt(k_x*k_x + k_y*k_y + k_z*k_z);

                    if (k > 0) {
                        k_box[half_box_idx(N, x, y, z)][0] *= sigma_func_vel_neutrino(k)/sigma_func_neutrino(k);
                        k_box[half_box_idx(N, x, y, z)][1] *= sigma_func_vel_neutrino(k)/sigma_func_neutrino(k);
                    }
                }
            }
        }

        //FTT back
        fftw_execute(c2r_plan);

        //Normalization (from Fourier conventions alone)
        for (int x=0; x<N; x++) {
            for (int y=0; y<N; y++) {
                for (int z=0; z<N; z++) {
                    primordial_box[box_idx(N, x, y, z)] /= box_volume;
                }
            }
        }

        //Write GRF to binary file
        write_array_to_disk(std::string(OUTPUT_DIR) + "gaussian_theta_nu.box", primordial_box, N);
        //Write GRF to H5 format
        writeGRF_H5(primordial_box, N, box_len, std::string(OUTPUT_DIR) + "gaussian_theta_nu.hdf5");

        std::cout << "1) The result has been written to " << std::string(OUTPUT_DIR) << "gaussian_theta_nu.box" << std::endl;
        std::cout << "2) The result has been written to " << std::string(OUTPUT_DIR) << "gaussian_theta_nu.hdf5" << std::endl;
        std::cout << "3) Theta is the divergence of the peculiar velocity in Mpc^-1," << std::endl;
        std::cout << "   which is (dx/dtau) = a*(dx/dt) with x=r/a comoving." << std::endl;
        std::cout << std::endl;

        //Compute the velocity field from the random field
        std::cout << "PHASE 3F - Compute the neutrino velocity field" << std::endl;
        std::cout << "1) Velocity convention v = a^2 (dx/dt), where x=r/a is comoving." << std::endl;
        std::cout << "2) Apply kernel to the Fourier transform of the divergence field." << std::endl;

        //FTT the primordial box
        fftw_execute(r2c_plan);

        //Normalization
        for (int x=0; x<N; x++) {
            for (int y=0; y<N; y++) {
                for (int z=0; z<=N/2; z++) {
                    k_box[half_box_idx(N, x, y, z)][0] *= box_volume/(N*N*N);
                    k_box[half_box_idx(N, x, y, z)][1] *= box_volume/(N*N*N);
                }
            }
        }

        for (int x=0; x<N; x++) {
            for (int y=0; y<N; y++) {
                for (int z=0; z<=N/2; z++) { //note that we stop at the (N/2+1)th entry
                    double k_x = (x > N/2) ? (x - N)*delta_k : x*delta_k; //Mpc^-1
                    double k_y = (y > N/2) ? (y - N)*delta_k : y*delta_k; //Mpc^-1
                    double k_z = (z > N/2) ? (z - N)*delta_k : z*delta_k; //Mpc^-1

                    double k = sqrt(k_x*k_x + k_y*k_y + k_z*k_z);

                    if (k>0) {
                        psi_k_x_box[half_box_idx(N, x, y, z)][0] = k_box[half_box_idx(N, x, y, z)][1] * k_x / (k*k);
                        psi_k_x_box[half_box_idx(N, x, y, z)][1] = -k_box[half_box_idx(N, x, y, z)][0] * k_x / (k*k);
                        psi_k_y_box[half_box_idx(N, x, y, z)][0] = k_box[half_box_idx(N, x, y, z)][1] * k_y / (k*k);
                        psi_k_y_box[half_box_idx(N, x, y, z)][1] = -k_box[half_box_idx(N, x, y, z)][0] * k_y / (k*k);
                        psi_k_z_box[half_box_idx(N, x, y, z)][0] = k_box[half_box_idx(N, x, y, z)][1] * k_z / (k*k);
                        psi_k_z_box[half_box_idx(N, x, y, z)][1] = -k_box[half_box_idx(N, x, y, z)][0] * k_z / (k*k);
                    } else {
                        psi_k_x_box[half_box_idx(N, x, y, z)][0] = 0;
                        psi_k_x_box[half_box_idx(N, x, y, z)][1] = 0;
                        psi_k_y_box[half_box_idx(N, x, y, z)][0] = 0;
                        psi_k_y_box[half_box_idx(N, x, y, z)][1] = 0;
                        psi_k_z_box[half_box_idx(N, x, y, z)][0] = 0;
                        psi_k_z_box[half_box_idx(N, x, y, z)][1] = 0;
                    }
                }
            }
        }

        std::cout << "3) Fourier transform back to real coordinates." << std::endl;

        //Do the IFFTs
        fftw_execute(px);
        fftw_execute(py);
        fftw_execute(pz);

        //Normalization
        for (int x=0; x<N; x++) {
            for (int y=0; y<N; y++) {
                for (int z=0; z<N; z++) {
                    psi_x_box[box_idx(N, x, y, z)] /= box_volume;
                    psi_y_box[box_idx(N, x, y, z)] /= box_volume;
                    psi_z_box[box_idx(N, x, y, z)] /= box_volume;
                }
            }
        }

        std::cout << "4) Done. The result consists of 3x" << N << "^3 real numbers." << std::endl;
        std::cout << "5) Make dimensionful by inserting c = " << c_vel << " Mpc/Gyr." << std::endl;

        //Insert units (convert from c=1 to dimensionful quantity)
        for (int x=0; x<N; x++) {
            for (int y=0; y<N; y++) {
                for (int z=0; z<N; z++) {
                    psi_x_box[box_idx(N, x, y, z)] *= c_vel;
                    psi_y_box[box_idx(N, x, y, z)] *= c_vel;
                    psi_z_box[box_idx(N, x, y, z)] *= c_vel;
                }
            }
        }

        double a = a_scale_factor_of_z(z_start);
        std::cout << "6) Multiply by the scale factor a = " << a << " to get v = a*(dx/tau) = a^2*(dx/dt)." << std::endl;

        //Insert scale factor
        for (int x=0; x<N; x++) {
            for (int y=0; y<N; y++) {
                for (int z=0; z<N; z++) {
                    psi_x_box[box_idx(N, x, y, z)] *= a;
                    psi_y_box[box_idx(N, x, y, z)] *= a;
                    psi_z_box[box_idx(N, x, y, z)] *= a;
                }
            }
        }

        write_array_to_disk(std::string(OUTPUT_DIR) + "v_nu_x.box", psi_x_box, N);
        write_array_to_disk(std::string(OUTPUT_DIR) + "v_nu_y.box", psi_y_box, N);
        write_array_to_disk(std::string(OUTPUT_DIR) + "v_nu_z.box", psi_z_box, N);

        std::cout << "   Output written to " << std::string(OUTPUT_DIR) + "v_nu_x.box" << std::endl;
        std::cout << "   Output written to " << std::string(OUTPUT_DIR) + "v_nu_y.box" << std::endl;
        std::cout << "   Output written to " << std::string(OUTPUT_DIR) + "v_nu_z.box" << std::endl;
        std::cout << "" << std::endl;


        std::cout << "PHASE 3G - Assign initial gravitational flow velocities to the neutrinos" << std::endl;

        for (auto body : bodies_nu) {
            double X = body.X*N/box_len;
            double Y = body.Y*N/box_len;
            double Z = body.Z*N/box_len;

            //Grid positions
            int iX = (int) floor(X);
            int iY = (int) floor(Y);
            int iZ = (int) floor(Z);

            //Intepolate the necessary fields with TSC
            float lookLength = 1.0;
            int lookLftX = (int) floor((X-iX) - lookLength);
            int lookRgtX = (int) floor((X-iX) + lookLength);
            int lookLftY = (int) floor((Y-iY) - lookLength);
            int lookRgtY = (int) floor((Y-iY) + lookLength);
            int lookLftZ = (int) floor((Z-iZ) - lookLength);
            int lookRgtZ = (int) floor((Z-iZ) + lookLength);

            //Accumulate interpolated values in psi_{xyz}
            double psi_x = 0, psi_y = 0, psi_z = 0;

            for (int i=lookLftX; i<=lookRgtX; i++) {
                for (int j=lookLftY; j<=lookRgtY; j++) {
                    for (int k=lookLftZ; k<=lookRgtZ; k++) {
                        //Pull the interpolated long-range force from the mesh
                        double xx = abs(X - (iX+i));
                        double yy = abs(Y - (iY+j));
                        double zz = abs(Z - (iZ+k));

                        double part_x = xx <= 1 ? 1-xx : 0;
                        double part_y = yy <= 1 ? 1-yy : 0;
                        double part_z = zz <= 1 ? 1-zz : 0;

                        psi_x += psi_x_box[box_wrap_idx(N, iX+i, iY+j, iZ+k)] * (part_x*part_y*part_z);
                        psi_y += psi_y_box[box_wrap_idx(N, iX+i, iY+j, iZ+k)] * (part_x*part_y*part_z);
                        psi_z += psi_z_box[box_wrap_idx(N, iX+i, iY+j, iZ+k)] * (part_x*part_y*part_z);
                    }
                }
            }

            body.v_X = psi_x;
            body.v_Y = psi_y;
            body.v_Z = psi_z;

            //Displace the particles if this has not yet been done
            if (!VELOCITY_AT_DISPLACED_POS) {
                body.X += body.delta_X;
                body.Y += body.delta_Y;
                body.Z += body.delta_Z;
            }

            bodies_nu[body.id] = body;
        }

        std::cout << "2) Assigned velocities to " << neutrino_num << " neutrino particles." << std::endl;
        if (!VELOCITY_AT_DISPLACED_POS) {
            std::cout << "3) Displaced " << particle_num << " cold particles." << std::endl;
        }
        std::cout << "  " << std::endl;
    } else {
        std::cout << "No valid method to determine the initial velocities. Exiting" << std::endl;
        return 0;
    }


    std::cout << "PHASE 4A - Add thermal motion to the neutrinos" << std::endl;

    double T = T_nu;
    double mu = mu_nu;
    sampler s;

    //Provide a reproducible random seed to the Fermi-Dirac sampler
    seed_rng(&s, seed+1);

    //Prepare the interpolation intevals of the quantile function
    prepare_intervals(&s, k_b*T, mu);

    double testdraw = sampler_draw(&s);

    std::cout << "1) Loading Fermi-Dirac sampler with T = " << T << " K, mu = " << mu << ", kT = " << k_b*T << " eV." << std::endl;
    std::cout << "2) The normalization is " << s.norm << "." << std::endl;
    std::cout << "3) Test draw: " << testdraw << " eV." << std::endl;

    double a_start = 1.0 / (1+z_start);
    double avg_speed = 0;
    double max_speed = 0;
    double avg_thermal_speed = 0;
    double max_thermal_speed = 0;
    double avg_bulk_speed = 0;

    for (auto body : bodies_nu) {
        //Generate a random speed V
        double draw = sampler_draw(&s); // E=pc in eV
        double p0 = draw * eV / c_vel * pow(Gyr/Mpc,2); // momentum in kg*Mpc/Gyr
        double p = p0/a_start; // redshifted momentum in kg*Mpc/Gyr
        double gamma = sqrt(1 + pow(p / (M_nu_kg*c_vel), 2)); // Lorentz factor
        double V = p/(gamma*M_nu_kg); // physical speed in Mpc/Gyr

        //Recall that our internal velocity variable is V = a^2(dx/dt),
        //where x=r/a is comoving. We therefore multiply by a.
        V *= a_start;

        //Just for diagnostics, look at the bulk speed before adding the thermal component
        double V_bulk = sqrt(body.v_X*body.v_X + body.v_Y*body.v_Y + body.v_Z*body.v_Z);

        //Generate a random point on the sphere
        double x = Gaussian(oracle);
        double y = Gaussian(oracle);
        double z = Gaussian(oracle);

        //And normalize
        double length = sqrt(x*x + y*y + z*z);
        if (length > 0) {
            x /= length;
            y /= length;
            z /= length;
        } else {
            //done since x=y=z=0
        }

        body.v_X += V*x;
        body.v_Y += V*y;
        body.v_Z += V*z;

        //For diagnostics, look at the new total speed
        double V_tot = sqrt(body.v_X*body.v_X + body.v_Y*body.v_Y + body.v_Z*body.v_Z);

        //Record statistics for diagnostics
        avg_speed += V_tot/neutrino_num;
        avg_thermal_speed += V/neutrino_num;
        avg_bulk_speed += V_bulk/neutrino_num;
        if (V > max_thermal_speed) {
            max_thermal_speed = V;
        }
        if (V_tot > max_speed) {
            max_speed = V_tot;
        }

        bodies_nu[body.id] = body;
    }

    std::cout << "4) Added thermal motion to " << neutrino_num << " particles." << std::endl;
    std::cout << "  " << std::endl;
    std::cout << "PHASE 4B - Summary" << std::endl;
    std::cout << "1) Average total speed " << avg_speed << " (internal) Mpc/Gyr or " << avg_speed/a_start/c_vel << " c physical." << std::endl;
    std::cout << "2) Maximum total speed " << max_speed << " (internal) Mpc/Gyr or " << max_speed/a_start/c_vel << " c physical." << std::endl;
    std::cout << "3) Average thermal speed " << avg_thermal_speed << " (internal) Mpc/Gyr or " << avg_thermal_speed/a_start/c_vel << " c physical." << std::endl;
    std::cout << "4) Maximum thermal speed " << max_thermal_speed << " (internal) Mpc/Gyr or " << max_thermal_speed/a_start/c_vel << " c physical." << std::endl;
    std::cout << "5) Average bulk flow speed " << avg_bulk_speed << " (internal) Mpc/Gyr." << std::endl;
    std::cout << "6) Internal speed means v = a^2 |dx/dt|, where x=r/a is comoving." << std::endl;
    std::cout << "  " << std::endl;

    /* NEXT, exporting the data in HDF5 format */

    std::cout << "PHASE 5 - Exporting initial conditions to HDF5 files." << std::endl;
    std::cout << "  " << std::endl;

    //Convert to Swift units. The below values are conversions to cgs
    double swift_unitcurrent = 1;
    double swift_unitlength = Mpc/cm;
    double swift_unitmass = M_swift;
    double swift_unittemp = 1;
    double swift_unittime = Gyr;

    std::cout << "Output length unit: Mpc (" << Mpc/cm << " cm)." << std::endl;
    std::cout << "Output time unit: Gyr (" << Gyr << " s)." << std::endl;
    std::cout << "Output mass unit: 10^10 M_sol (" << swift_unitmass << " g)." << std::endl;
    std::cout << "  " << std::endl;

    //Convert the masses (lengths and times are unchanged)
    for (int i=0; i<particle_num; i++) {
        bodies[i].mass /= (swift_unitmass/kg); //from kg to swift unit
    }

    //Convert the neutrino masses (lengths and times are unchanged)
    for (int i=0; i<neutrino_num; i++) {
        bodies_nu[i].mass /= (swift_unitmass/kg); //from kg to swift unit
    }

    std::cout << "Converted masses from kg to U_M = " << swift_unitmass << " g." << std::endl;

    //Convert generalized velocities a^2(dx/dt) to peculiar velocities a(dx/dt)
    if (OUTPUT_VELOCITY_FORMAT == PECULIAR_VELOCITY) {
        double a = 1.0/(1.0 + z_start);
        for (int i=0; i<particle_num; i++) {
            bodies[i].v_X /= a;
            bodies[i].v_Y /= a;
            bodies[i].v_Z /= a;
        }

        for (int i=0; i<neutrino_num; i++) {
            bodies_nu[i].v_X /= a;
            bodies_nu[i].v_Y /= a;
            bodies_nu[i].v_Z /= a;
        }

        std::cout << "Converted generalized velocities a^2(dx/dt) to peculiar velocities a*(dx/dt)." << std::endl;
        std::cout << "  " << std::endl;
    }

    double swift_rho_crit = rho_crit * pow(Mpc,3) / (swift_unitmass/kg);
    std::cout << "SUMMARY" << std::endl;
    std::cout << "Inferred Omega_m = " << particle_num * bodies[0].mass / (swift_rho_crit*box_volume) << std::endl;
    std::cout << "Inferred Omega_nu = " << neutrino_num * bodies_nu[0].mass / (swift_rho_crit*box_volume) << std::endl;
    std::cout << "Total cold mass (cdm+b) " << particle_num * bodies[0].mass << std::endl;
    std::cout << "Total neutrino mass " << neutrino_num * bodies_nu[0].mass << std::endl;
    std::cout << "Total mass " << (particle_num * bodies[0].mass + neutrino_num * bodies_nu[0].mass) << std::endl;
    std::cout << "Volume " << box_volume << std::endl;
    std::cout << "Rho crit " << swift_rho_crit << std::endl;

    //The particle file
    std::string fname = std::string(OUTPUT_DIR) + "particles.hdf5";

    //Convert to char array to please H5
    char fname_chars[fname.size() + 1];
    fname.copy(fname_chars, fname.size() + 1);
    fname_chars[fname.size()] = '\0';

    //Export the particles to an HDF5 file
    H5::H5File file(fname_chars, H5F_ACC_TRUNC);

    //Write Header group
    H5::Group headerGroup(file.createGroup("/Header"));

    //Dataspace for a single number
    H5::DataSpace scalarSpace(H5S_SCALAR);
    //Data space for a row of 7 numbers (one per particle type)
    const std::size_t ptNDIMS = 1;
    hsize_t ptdims[ptNDIMS] = {7};
    H5::DataSpace row6Space(ptNDIMS, ptdims);
    //Data space for a row of 3 numbers (one per coordinate)
    const std::size_t cNDIMS = 1;
    hsize_t cdims[cNDIMS] = {3};
    H5::DataSpace row3Space(cNDIMS, cdims);

    //Write all the header attributes to the file
    H5::Attribute att0 = headerGroup.createAttribute("BoxSize", H5::PredType::NATIVE_FLOAT, row3Space);
    float boxsizes[3] = {box_len, box_len, box_len};
    att0.write(H5::PredType::NATIVE_FLOAT, boxsizes);

    H5::Attribute att1 = headerGroup.createAttribute("Dimension", H5::PredType::NATIVE_INT, scalarSpace);
    int coordinateDim = 3;
    att1.write(H5::PredType::NATIVE_INT, &coordinateDim);

    H5::Attribute att2 = headerGroup.createAttribute("Flag_Entropy_ICs", H5::PredType::NATIVE_INT, scalarSpace);
    int entropyFlag = 0;
    att2.write(H5::PredType::NATIVE_INT, &entropyFlag);

    H5::Attribute att3 = headerGroup.createAttribute("MassTable", H5::PredType::NATIVE_FLOAT, row6Space);
    float massTable[7] = {0,0,0,0,0,0,0};
    att3.write(H5::PredType::NATIVE_FLOAT, massTable);

    H5::Attribute att4 = headerGroup.createAttribute("NumFilesPerSnapshot", H5::PredType::NATIVE_INT, scalarSpace);
    int NumFilesPerSnapshot = 1;
    att4.write(H5::PredType::NATIVE_INT, &NumFilesPerSnapshot);

    H5::Attribute att5 = headerGroup.createAttribute("NumPart_Total", H5::PredType::NATIVE_INT, row6Space);
    int numparttotal[7] = {0,PARTICLE_NUM,0,0,0,0,NEUTRINO_NUM}; //number particles of type 0,1,2,3,4,5,6 (3 is dummy, doesn't exist)
    att5.write(H5::PredType::NATIVE_INT, &numparttotal);

    H5::Attribute att6 = headerGroup.createAttribute("NumPart_Total_HighWord", H5::PredType::NATIVE_INT, row6Space);
    int numparttotalhighword[7] = {0,0,0,0,0,0,0};
    att6.write(H5::PredType::NATIVE_INT, &numparttotalhighword);

    H5::Attribute att7 = headerGroup.createAttribute("Time", H5::PredType::NATIVE_FLOAT, scalarSpace);
    float Time = 0;
    att7.write(H5::PredType::NATIVE_FLOAT, &Time);


    //Choose particle group (dark matter is particle type 1)
    H5::Group group(file.createGroup("/PartType1"));

    float positions[3*PARTICLE_NUM];
    float internal_energy[PARTICLE_NUM];
    float masses[PARTICLE_NUM];
    long int particleIDs[PARTICLE_NUM];
    float smoothingLength[PARTICLE_NUM];
    float velocities[3*PARTICLE_NUM];

    for (int i=0; i<PARTICLE_NUM; i++) {
        positions[0+i*3] = bodies[i].X;
        positions[1+i*3] = bodies[i].Y;
        positions[2+i*3] = bodies[i].Z;

        internal_energy[i] = 1.23734;
        masses[i] = bodies[i].mass;
        particleIDs[i] = bodies[i].id;
        smoothingLength[i] = 2.4696;

        velocities[0+i*3] = bodies[i].v_X;
        velocities[1+i*3] = bodies[i].v_Y;
        velocities[2+i*3] = bodies[i].v_Z;
    }

    //std::cout << bodies[PARTICLE_NUM-1].X << " " << bodies[PARTICLE_NUM-1].Y << " " << bodies[PARTICLE_NUM-1].Z << std::endl;

    //Create two dataspaces (one 3xPARTICLE_NUM and one row 1xPARTICLE_NUM)
    const std::size_t NDIMS = 2;
    const std::size_t pNDIMS = 1;
    hsize_t dims[NDIMS] = {PARTICLE_NUM,3};
    hsize_t pdims[pNDIMS] = {PARTICLE_NUM};
    H5::DataSpace dataspace(NDIMS, dims);
    H5::DataSpace particle_dataspace(pNDIMS, pdims);

    //Write the position array to disk
    H5::DataSet dataset = group.createDataSet("Coordinates", H5::PredType::NATIVE_FLOAT, dataspace);
    dataset.write(positions, H5::PredType::NATIVE_FLOAT);

    // //Write the internal_energy array to disk
    // H5::DataSet dataset_energy = group.createDataSet("InternalEnergy", H5::PredType::NATIVE_FLOAT, particle_dataspace);
    // dataset_energy.write(internal_energy, H5::PredType::NATIVE_FLOAT);

    //Write the masses array to disk
    H5::DataSet dataset_mass = group.createDataSet("Masses", H5::PredType::NATIVE_FLOAT, particle_dataspace);
    dataset_mass.write(masses, H5::PredType::NATIVE_FLOAT);

    //Write the particleIDs array to disk
    H5::DataSet dataset_pid = group.createDataSet("ParticleIDs", H5::PredType::NATIVE_LONG, particle_dataspace);
    dataset_pid.write(particleIDs, H5::PredType::NATIVE_LONG);

    // //Write the smoothingLength array to disk
    // H5::DataSet dataset_sl = group.createDataSet("SmoothingLength", H5::PredType::NATIVE_FLOAT, particle_dataspace);
    // dataset_sl.write(smoothingLength, H5::PredType::NATIVE_FLOAT);

    //Write the velocity array to disk
    H5::DataSet dataset_v = group.createDataSet("Velocities", H5::PredType::NATIVE_FLOAT, dataspace);
    dataset_v.write(velocities, H5::PredType::NATIVE_FLOAT);


    /* Next, export the neutrinos */

    //Choose particle group (neutrinos is particle type 6)
    H5::Group group_nu(file.createGroup("/PartType6"));

    float positions_nu[3*NEUTRINO_NUM];
    float internal_energy_nu[NEUTRINO_NUM];
    float masses_nu[NEUTRINO_NUM];
    long int particleIDs_nu[NEUTRINO_NUM];
    float smoothingLength_nu[NEUTRINO_NUM];
    float velocities_nu[3*NEUTRINO_NUM];

    for (int i=0; i<NEUTRINO_NUM; i++) {
        positions_nu[0+i*3] = bodies_nu[i].X;
        positions_nu[1+i*3] = bodies_nu[i].Y;
        positions_nu[2+i*3] = bodies_nu[i].Z;

        internal_energy_nu[i] = 1.23734;
        masses_nu[i] = bodies_nu[i].mass;
        particleIDs_nu[i] = bodies_nu[i].id + PARTICLE_NUM;
        smoothingLength_nu[i] = 2.4696;

        velocities_nu[0+i*3] = bodies_nu[i].v_X;
        velocities_nu[1+i*3] = bodies_nu[i].v_Y;
        velocities_nu[2+i*3] = bodies_nu[i].v_Z;
    }

    //std::cout << bodies[PARTICLE_NUM-1].X << " " << bodies[PARTICLE_NUM-1].Y << " " << bodies[PARTICLE_NUM-1].Z << std::endl;

    //Create two dataspaces (one 3xPARTICLE_NUM and one row 1xPARTICLE_NUM)
    hsize_t dims_nu[NDIMS] = {NEUTRINO_NUM,3};
    hsize_t pdims_nu[pNDIMS] = {NEUTRINO_NUM};
    H5::DataSpace dataspace_nu(NDIMS, dims_nu);
    H5::DataSpace particle_dataspace_nu(pNDIMS, pdims_nu);

    //Write the position array to disk
    H5::DataSet dataset_nu = group_nu.createDataSet("Coordinates", H5::PredType::NATIVE_FLOAT, dataspace_nu);
    dataset_nu.write(positions_nu, H5::PredType::NATIVE_FLOAT);

    //Write the masses array to disk
    H5::DataSet dataset_mass_nu = group_nu.createDataSet("Masses", H5::PredType::NATIVE_FLOAT, particle_dataspace_nu);
    dataset_mass_nu.write(masses_nu, H5::PredType::NATIVE_FLOAT);

    //Write the particleIDs array to disk
    H5::DataSet dataset_pid_nu = group_nu.createDataSet("ParticleIDs", H5::PredType::NATIVE_LONG, particle_dataspace_nu);
    dataset_pid_nu.write(particleIDs_nu, H5::PredType::NATIVE_LONG);

    //Write the velocity array to disk
    H5::DataSet dataset_v_nu = group_nu.createDataSet("Velocities", H5::PredType::NATIVE_FLOAT, dataspace_nu);
    dataset_v_nu.write(velocities_nu, H5::PredType::NATIVE_FLOAT);

    //Write more attributes

    //Write RuntimePars group
    H5::Group runtimeGroup(file.createGroup("/RuntimePars"));

    H5::Attribute rt_att0 = runtimeGroup.createAttribute("PeriodicBoundariesOn", H5::PredType::NATIVE_INT, scalarSpace);
    int PeriodicBoundariesOn = 1;
    rt_att0.write(H5::PredType::NATIVE_INT, &PeriodicBoundariesOn);

    //Write Units group
    H5::Group unitsGroup(file.createGroup("/Units"));

    H5::Attribute units_att0 = unitsGroup.createAttribute("Unit current in cgs (U_I)", H5::PredType::NATIVE_DOUBLE, scalarSpace);
    units_att0.write(H5::PredType::NATIVE_DOUBLE, &swift_unitcurrent);

    H5::Attribute units_att1 = unitsGroup.createAttribute("Unit length in cgs (U_L)", H5::PredType::NATIVE_DOUBLE, scalarSpace);
    units_att1.write(H5::PredType::NATIVE_DOUBLE, &swift_unitlength);

    H5::Attribute units_att2 = unitsGroup.createAttribute("Unit mass in cgs (U_M)", H5::PredType::NATIVE_DOUBLE, scalarSpace);
    units_att2.write(H5::PredType::NATIVE_DOUBLE, &swift_unitmass);

    H5::Attribute units_att3 = unitsGroup.createAttribute("Unit temperature in cgs (U_T)", H5::PredType::NATIVE_DOUBLE, scalarSpace);
    units_att3.write(H5::PredType::NATIVE_DOUBLE, &swift_unittemp);

    H5::Attribute units_att4 = unitsGroup.createAttribute("Unit time in cgs (U_t)", H5::PredType::NATIVE_DOUBLE, scalarSpace);
    units_att4.write(H5::PredType::NATIVE_DOUBLE, &swift_unittime);

    file.close();
}


void writeGRF_H5(double *box, size_t N, float box_len, std::string fname) {
    //Dataspace for a single number
    H5::DataSpace scalarSpace(H5S_SCALAR);
    //Data space for a row of 3 numbers (one per coordinate)
    const std::size_t cNDIMS = 1;
    hsize_t cdims[cNDIMS] = {3};
    H5::DataSpace row3Space(cNDIMS, cdims);

    //Convert to char array to please H5
    char fname_chars[fname.size() + 1];
    fname.copy(fname_chars, fname.size() + 1);
    fname_chars[fname.size()] = '\0';

    //Export the primordial Gaussian random field to another HDF5 file
    H5::H5File grf_file(fname_chars, H5F_ACC_TRUNC);

    //Write Header group
    H5::Group grf_headerGroup(grf_file.createGroup("/Header"));

    //Write all the header attributes to the file
    H5::Attribute grf_att0 = grf_headerGroup.createAttribute("BoxSize", H5::PredType::NATIVE_FLOAT, row3Space);
    float grf_boxsizes[3] = {box_len, box_len, box_len};
    grf_att0.write(H5::PredType::NATIVE_FLOAT, grf_boxsizes);

    H5::Attribute grf_att1 = grf_headerGroup.createAttribute("Dimension", H5::PredType::NATIVE_INT, scalarSpace);
    int grf_coordinateDim = 3;
    grf_att1.write(H5::PredType::NATIVE_INT, &grf_coordinateDim);

    //Write Field group
    H5::Group grf_fieldGroup(grf_file.createGroup("/Field"));

    //Create an N^3 dataspace
    const std::size_t grf_NDIMS = 3;
    hsize_t grf_dims[grf_NDIMS] = {N,N,N};
    H5::DataSpace grf_dataspace(grf_NDIMS, grf_dims);


    //Copy the Gaussian random field from FFTW format to a suitable formatted array
    float gaussian_arr[N*N*N];
    for (int x=0; x<N; x++) {
        for (int y=0; y<N; y++) {
            for (int z=0; z<N; z++) {
                gaussian_arr[box_idx(N, x, y, z)] = box[box_idx(N, x, y, z)];
            }
        }
    }

    //Write the gaussian field array to disk
    H5::DataSet grf_dataset = grf_fieldGroup.createDataSet("GaussianRandomField", H5::PredType::NATIVE_FLOAT, grf_dataspace);
    grf_dataset.write(gaussian_arr, H5::PredType::NATIVE_FLOAT);

    grf_file.close();
}
