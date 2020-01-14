/**
 * @infini_lpt.cpp
 * @author  Willem Elbers <whe@willemelbers.com>
 * @version 1.0
 * @license MIT (https://opensource.org/licenses/MIT)
 *
 * @section DESCRIPTION
 *
 * Infinite Lagrangian perturbation theory.
 */

 #include "infini_lpt.h"




/**
 * Function do_infini_lpt()
 * Iteratively solve the Monge-Ampère equation for use in infinite LPT.
 *
 * @param phi_box Array of the potential which is to be solved for
 * @param src_box Array of the source function of the Monche-Ampère equation,
 *        i.e. the overdensity
 * @param width The dimension of the grid (=N)
 * @param box_len Physical dimension of the box in Mpc
 * @param tol Tolerance parameter
 */
void do_infini_lpt(double* phi_box, double* src_box, int width, double box_len, double tol) {
	double box_volume = pow(box_len, 3);
	double factor = sqrt(box_volume/2.0); //see fourier.pdf
	double delta_k = 2*M_PI/box_len; //Mpc^-1
	int half_width = width/2;
    int N = width;

	if (width%2 == 1) {
		throw std::invalid_argument("We expected an even dimension.");
	}

    //Boxes and FFT plans for the potential phi
    fftw_complex *phi_k_box = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*N*(N/2+1));
    fftw_plan phi_c2r  = fftw_plan_dft_c2r_3d(N, N, N, phi_k_box, phi_box, FFTW_ESTIMATE);

    //Boxes and FFT plans for the difference of (1+source) - det[I + Hessian(phi)]
    double *dif_box = (double*) fftw_malloc(sizeof(double)*N*N*N);
    fftw_complex *dif_k_box = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*N*(N/2+1));
    fftw_plan dif_r2c  = fftw_plan_dft_r2c_3d(N, N, N, dif_box,   dif_k_box, FFTW_ESTIMATE);

    //Boxes and FFT plans for the Hessian of phi
    double *H_xx_box = (double*) fftw_malloc(sizeof(double)*N*N*N);
    double *H_yy_box = (double*) fftw_malloc(sizeof(double)*N*N*N);
    double *H_zz_box = (double*) fftw_malloc(sizeof(double)*N*N*N);
    double *H_xy_box = (double*) fftw_malloc(sizeof(double)*N*N*N);
    double *H_xz_box = (double*) fftw_malloc(sizeof(double)*N*N*N);
    double *H_yz_box = (double*) fftw_malloc(sizeof(double)*N*N*N);
    fftw_complex *H_xx_k_box = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*N*(N/2+1));
    fftw_complex *H_yy_k_box = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*N*(N/2+1));
    fftw_complex *H_zz_k_box = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*N*(N/2+1));
    fftw_complex *H_xy_k_box = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*N*(N/2+1));
    fftw_complex *H_xz_k_box = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*N*(N/2+1));
    fftw_complex *H_yz_k_box = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*N*(N/2+1));
    fftw_plan H_xx_c2r  = fftw_plan_dft_c2r_3d(N, N, N, H_xx_k_box, H_xx_box, FFTW_ESTIMATE);
    fftw_plan H_yy_c2r  = fftw_plan_dft_c2r_3d(N, N, N, H_yy_k_box, H_yy_box, FFTW_ESTIMATE);
    fftw_plan H_zz_c2r  = fftw_plan_dft_c2r_3d(N, N, N, H_zz_k_box, H_zz_box, FFTW_ESTIMATE);
    fftw_plan H_xy_c2r  = fftw_plan_dft_c2r_3d(N, N, N, H_xy_k_box, H_xy_box, FFTW_ESTIMATE);
    fftw_plan H_xz_c2r  = fftw_plan_dft_c2r_3d(N, N, N, H_xz_k_box, H_xz_box, FFTW_ESTIMATE);
    fftw_plan H_yz_c2r  = fftw_plan_dft_c2r_3d(N, N, N, H_yz_k_box, H_yz_box, FFTW_ESTIMATE);

    //Initialize the difference box using phi = 0 initially
    for (int x=0; x<N; x++) {
        for (int y=0; y<N; y++) {
            for (int z=0; z<N; z++) {
                dif_box[box_idx(N, x, y, z)] = src_box[box_idx(N, x, y, z)];
            }
        }
    }

    //Initialize the phi_k_box with phi = 0
    for (int x=0; x<N; x++) {
        for (int y=0; y<N; y++) {
            for (int z=0; z<=N/2; z++) {
                phi_k_box[half_box_idx(N, x, y, z)][0] = 0;
                phi_k_box[half_box_idx(N, x, y, z)][1] = 0;
            }
        }
    }

    double err1 = 0, err2 = 0;
    int ITER = 0;
    bool done = false;

    while(!done) {
        fftw_execute(dif_r2c);

        //Normalization
        for (int x=0; x<N; x++) {
            for (int y=0; y<N; y++) {
                for (int z=0; z<=N/2; z++) {
                    dif_k_box[half_box_idx(N, x, y, z)][0] *= box_volume/(N*N*N);
                    dif_k_box[half_box_idx(N, x, y, z)][1] *= box_volume/(N*N*N);
                }
            }
        }


        //Multiply with the inverse Poisson kernel and differentiate
        for (int x=0; x<N; x++) {
            for (int y=0; y<N; y++) {
                for (int z=0; z<=N/2; z++) { //note that we stop at the (N/2+1)th entry
                    double k_x = (x > N/2) ? (x - N)*delta_k : x*delta_k; //Mpc^-1
                    double k_y = (y > N/2) ? (y - N)*delta_k : y*delta_k; //Mpc^-1
                    double k_z = (z > N/2) ? (z - N)*delta_k : z*delta_k; //Mpc^-1

                    double k = sqrt(k_x*k_x + k_y*k_y + k_z*k_z);

                    if (k>0) {
                        phi_k_box[half_box_idx(N, x, y, z)][0] += dif_k_box[half_box_idx(N, x, y, z)][0] / (k*k);
                        phi_k_box[half_box_idx(N, x, y, z)][1] += dif_k_box[half_box_idx(N, x, y, z)][1] / (k*k);

                        H_xx_k_box[half_box_idx(N, x, y, z)][0] = phi_k_box[half_box_idx(N, x, y, z)][0] * (k_x * k_x);
                        H_xx_k_box[half_box_idx(N, x, y, z)][1] = phi_k_box[half_box_idx(N, x, y, z)][1] * (k_x * k_x);
                        H_yy_k_box[half_box_idx(N, x, y, z)][0] = phi_k_box[half_box_idx(N, x, y, z)][0] * (k_y * k_y);
                        H_yy_k_box[half_box_idx(N, x, y, z)][1] = phi_k_box[half_box_idx(N, x, y, z)][1] * (k_y * k_y);
                        H_zz_k_box[half_box_idx(N, x, y, z)][0] = phi_k_box[half_box_idx(N, x, y, z)][0] * (k_z * k_z);
                        H_zz_k_box[half_box_idx(N, x, y, z)][1] = phi_k_box[half_box_idx(N, x, y, z)][1] * (k_z * k_z);
                        H_xy_k_box[half_box_idx(N, x, y, z)][0] = phi_k_box[half_box_idx(N, x, y, z)][0] * (k_x * k_y);
                        H_xy_k_box[half_box_idx(N, x, y, z)][1] = phi_k_box[half_box_idx(N, x, y, z)][1] * (k_x * k_y);
                        H_xz_k_box[half_box_idx(N, x, y, z)][0] = phi_k_box[half_box_idx(N, x, y, z)][0] * (k_x * k_z);
                        H_xz_k_box[half_box_idx(N, x, y, z)][1] = phi_k_box[half_box_idx(N, x, y, z)][1] * (k_x * k_z);
                        H_yz_k_box[half_box_idx(N, x, y, z)][0] = phi_k_box[half_box_idx(N, x, y, z)][0] * (k_y * k_z);
                        H_yz_k_box[half_box_idx(N, x, y, z)][1] = phi_k_box[half_box_idx(N, x, y, z)][1] * (k_y * k_z);
                    } else {
                        phi_k_box[half_box_idx(N, x, y, z)][0] = 0;
                        phi_k_box[half_box_idx(N, x, y, z)][1] = 0;

                        H_xx_k_box[half_box_idx(N, x, y, z)][0] = 0;
                        H_xx_k_box[half_box_idx(N, x, y, z)][1] = 0;
                        H_yy_k_box[half_box_idx(N, x, y, z)][0] = 0;
                        H_yy_k_box[half_box_idx(N, x, y, z)][1] = 0;
                        H_zz_k_box[half_box_idx(N, x, y, z)][0] = 0;
                        H_zz_k_box[half_box_idx(N, x, y, z)][1] = 0;
                        H_xy_k_box[half_box_idx(N, x, y, z)][0] = 0;
                        H_xy_k_box[half_box_idx(N, x, y, z)][1] = 0;
                        H_xz_k_box[half_box_idx(N, x, y, z)][0] = 0;
                        H_xz_k_box[half_box_idx(N, x, y, z)][1] = 0;
                        H_yz_k_box[half_box_idx(N, x, y, z)][0] = 0;
                        H_yz_k_box[half_box_idx(N, x, y, z)][1] = 0;
                    }
                }
            }
        }

        fftw_execute(H_xx_c2r);
        fftw_execute(H_yy_c2r);
        fftw_execute(H_zz_c2r);
        fftw_execute(H_xy_c2r);
        fftw_execute(H_xz_c2r);
        fftw_execute(H_yz_c2r);

        //Normalization
        for (int x=0; x<N; x++) {
            for (int y=0; y<N; y++) {
                for (int z=0; z<N; z++) {
                    H_xx_box[box_idx(N, x, y, z)] /= box_volume;
                    H_yy_box[box_idx(N, x, y, z)] /= box_volume;
                    H_zz_box[box_idx(N, x, y, z)] /= box_volume;
                    H_xy_box[box_idx(N, x, y, z)] /= box_volume;
                    H_xz_box[box_idx(N, x, y, z)] /= box_volume;
                    H_yz_box[box_idx(N, x, y, z)] /= box_volume;
                }
            }
        }

        Matrix3d Hessian, Eye = Matrix3d::Identity();
        double mean_dif = 0;

        for (int x=0; x<N; x++) {
            for (int y=0; y<N; y++) {
                for (int z=0; z<N; z++) {
                    int idx = box_idx(N, x, y, z);

                    Hessian << H_xx_box[idx], H_xy_box[idx], H_xz_box[idx],
            				   H_xy_box[idx], H_yy_box[idx], H_yz_box[idx],
                               H_xz_box[idx], H_yz_box[idx], H_zz_box[idx];

                    double Det = (Hessian + Eye).determinant();
                    double rho = 1.0 + src_box[idx];
                    dif_box[idx] = rho - Det;
                    mean_dif += (rho - Det)/(N*N*N);
                }
            }
        }

        //Get rid of the constant mode for diagnostic purposes (doesn't affect
        // phi, but it reveals the real difference)
        double var = 0;
        double rho_var = 0;
        double D_var = 0;
        for (int x=0; x<N; x++) {
            for (int y=0; y<N; y++) {
                for (int z=0; z<N; z++) {
                    dif_box[box_idx(N, x, y, z)] -= mean_dif;
                    var += pow(dif_box[box_idx(N, x, y, z)],2)/(N*N*N-1);

                    int idx = box_idx(N, x, y, z);

                    Hessian << H_xx_box[idx], H_xy_box[idx], H_xz_box[idx],
                               H_xy_box[idx], H_yy_box[idx], H_yz_box[idx],
                               H_xz_box[idx], H_yz_box[idx], H_zz_box[idx];

                    double Det = (Hessian + Eye).determinant();
                    double rho = 1.0 + src_box[idx];

                    rho_var += (rho-1)*(rho-1)/(N*N*N-1);
                    D_var += (Det-1+mean_dif)*(Det-1+mean_dif)/(N*N*N-1);
                }
            }
        }

        err2 = err1;
        err1 = sqrt(var);

        if (MONGE_AMPERE_VERBOSE) {
            std::cout << "" << ITER << ") eps = " << err1 << " \t relative " << (D_var-rho_var)/rho_var << "." << std::endl;
        }

        if (ITER >= MONGE_AMPERE_MAX_ITER || err1 <= tol) {
            done = true;
        } else {
            ITER++;
        }
    }

    if (err1 > tol) {
        throw std::logic_error("Maximum number of iterations exceeded, tolerance not reached.");
    } else if (MONGE_AMPERE_VERBOSE) {
        std::cout << (ITER+1) << ") MAE solver converged." << std::endl;
    }

    //Perform the final IFFT to get phi
    fftw_execute(phi_c2r);

    //Normalization
    for (int x=0; x<N; x++) {
        for (int y=0; y<N; y++) {
            for (int z=0; z<N; z++) {
                phi_box[box_idx(N, x, y, z)] /= box_volume;
            }
        }
    }
}
