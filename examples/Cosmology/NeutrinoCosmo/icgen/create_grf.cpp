/**
 * @create_grf.cpp
 * @author  Willem Elbers <whe@willemelbers.com>
 * @version 1.0
 * @license MIT (https://opensource.org/licenses/MIT)
 *
 * @section DESCRIPTION
 *
 * Generate a 3D Gaussian Random Field.
 */

 #include "create_grf.h"


long int sayhello() {
    return 7369767679;
}

/**
 * Function generate_grf()
 * Generate the Fourier transform of a Gaussian random field with specified power spectrum.
 * A single complex-to-real DFT will produce the desired real-valued random field.
 * Note that k_box has dimension N*N*(N/2+1), exploiting the conjugate symmetry relation.
 *
 * @param oracle Reference to a pre-seeded std::default_random_engine
 * @param k_box Array of complex numbers containing the Fourier transform of the field
 * @param width The dimension of the grid (=N)
 * @param box_len Physical dimension of the box in Mpc
 * @param sigma_func Reference to a function specifying the square root of the power spectrum
 */
void generate_grf(std::default_random_engine& oracle, fftw_complex* k_box, int width, double box_len, double (&sigma_func)(double)) {
	double box_volume = pow(box_len, 3);
	double factor = sqrt(box_volume/2.0); //see fourier.pdf
	double delta_k = 2*M_PI/box_len; //Mpc^-1
	int half_width = width/2;

	if (width%2 == 1) {
		throw std::invalid_argument("We expected an even dimension.");
	}

	std::normal_distribution<double> Gaussian(0.0, 1.0);

 	for (int x=0; x<width; x++) {
		for (int y=0; y<width; y++) {
			for (int z=0; z<=half_width; z++) { //note that we stop at the (N/2+1)th entry
				double k_x = (x > half_width) ? (x - width)*delta_k : x*delta_k; //Mpc^-1
				double k_y = (y > half_width) ? (y - width)*delta_k : y*delta_k; //Mpc^-1
				double k_z = (z > half_width) ? (z - width)*delta_k : z*delta_k; //Mpc^-1

				double k = sqrt(k_x*k_x + k_y*k_y + k_z*k_z);

				//The DC mode needs special treatment
				if (k==0) {
					k_box[half_box_idx(width, x, y, z)][0] = 0;
					k_box[half_box_idx(width, x, y, z)][1] = 0;
				} else {
					double sigma = sigma_func(k);
					double a = Gaussian(oracle) * sigma;
					double b = Gaussian(oracle) * sigma;

					k_box[half_box_idx(width, x, y, z)][0] = a * factor;
					k_box[half_box_idx(width, x, y, z)][1] = b * factor;
				}

			}
		}
	}
}


/**
 * Function integrate_sigma_R()
 * Calculate the value of sigma_R (e.g. sigma_8), corresponding to a given power spectrum,
 * by integrating over the whole Fourier domain.
 *
 * @param width The dimension of the grid (=N) --- not used
 * @param box_len Physical dimension of the box in Mpc --- not used
 * @param R_filter Filtering scale in Mpc
 * @param sigma_func Reference to a function specifying the square root of the power spectrum
 */
double integrate_sigma_R(int width, double box_len, double R_filter, double (&sigma_func)(double)) {
	double total = 0;

	//Integrate sigma^2(k) between k_min and k_max
	double k_max = 10.0; // 1/Mpc
	double k_min = 0.0001; // 1/Mpc

	double k = k_max; // 1/Mpc
	double k_scaling = 1.001;

	//Output debug file
	std::ofstream of(std::string(OUTPUT_DIR) + "sigma_8_integration.txt");
	of << "k(1/Mpc);W(k);P(k);partial_sum(sigma^2)\n";

 	while(k>k_min) {
		double W = 3.0 / pow(k*R_filter, 3) *(sin(k*R_filter) - k*R_filter*cos(k*R_filter)); //top-hat filter
		double Delta_k = k-k/k_scaling;
		total += Delta_k * (W*W) * pow(sigma_func(k),2) * pow(k,2) / (2*pow(M_PI,2));
		of << k << ";" << W << ";" << pow(sigma_func(k),2) << ";" << total << "\n";
		k /= k_scaling;
	}

	//Close file
	of.close();

	return sqrt(total);
}
