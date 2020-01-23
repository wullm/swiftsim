/**
 * @calc_powerspec.cpp
 * @author  Willem Elbers <whe@willemelbers.com>
 * @version 1.0
 *
 * @section DESCRIPTION
 *
 * Calculates the 3D power spectrum from a grid using DFT.
 */


double sinc(double x) {
    return (x==0) ? 0 : x/sin(x);
}

/**
 * Function calc_sigma_R()
 * Calculates fluctuations over a given (top-hat) filter scale, assuming the power
 * spectrum has just been calculated. See fourier.pdf. Most useful for sigma_8.
 *
 * @param R_filter The filtering scale in Mpc
 * @param bins The number of logarithmic k bins
 * @param k_in_bins Reference to an array containing the average value of k in each bin
 * @param power_in_bins Reference to an array containing the average power of k in each bin
 * @param obs_in_bins Reference to an array containing the number of observations in each bin
 */
double calc_sigma_R(double R_filter, int bins, double *k_in_bins, double *power_in_bins, int *obs_in_bins) {
	double sigma_8_squared = 0;
	for (int i=0; i<bins-1; i++) {
		if (obs_in_bins[i]>0 && obs_in_bins[i+1]) {
			double k = k_in_bins[i];
			double delta_k = k_in_bins[i+1] - k_in_bins[i];
			double W_R = 3.0 / pow(k*R_filter, 3) *(sin(k*R_filter) - k*R_filter*cos(k*R_filter)); //top-hat filter
			sigma_8_squared += delta_k * abs(W_R)*abs(W_R) * power_in_bins[i] * pow(k,2) / (2*pow(M_PI,2));
		}
	}

	return sqrt(sigma_8_squared);
}


/**
 * Function calc_cross_powerspec()
 * Calculate the cross power spectrum of two fields A and B, defined by
 * P_AB(k) * (2pi)^3 * delta^3(k+k') = <A(k)B(k')^*>. See fourier.pdf.
 *
 * @param width The dimension of the grid (=N)
 * @param box_len Physical dimension of the box in Mpc
 * @param k_box1 (Field A) Reference to an array of N*N*(N/2+1) complex numbers
 * @param k_box2 (Field B) Reference to an array of N*N*(N/2+1) complex numbers
 * @param bins The number of logarithmic k bins
 * @param k_in_bins Reference to an array containing the average value of k in each bin
 * @param power_in_bins Reference to an array containing the average power of k in each bin
 * @param obs_in_bins Reference to an array containing the number of observations in each bin
 * @param deconvolve (optional) The window function that should be decolved (0=no deconvolution; default)
 */
void calc_cross_powerspec(int width, double box_len, fftw_complex *k_box1, fftw_complex *k_box2, int bins, double *k_in_bins, double *power_in_bins, int *obs_in_bins, int deconvolve=0) {
	int N = width, half_width = N/2;

	double box_volume = pow(box_len, 3); //Mpc^3
	double delta_k = 2*M_PI/box_len; //Mpc^-1

	double max_k = sqrt(3)*delta_k*half_width;
	double min_k = delta_k;

	for (int x=0; x<width; x++) {
		for (int y=0; y<width; y++) {
			for (int z=0; z<=half_width; z++) {
				double k_x = (x > half_width) ? (x - width)*delta_k : x*delta_k;
				double k_y = (y > half_width) ? (y - width)*delta_k : y*delta_k;
				double k_z = (z > half_width) ? (z - width)*delta_k : z*delta_k;

				double k = sqrt(k_x*k_x + k_y*k_y + k_z*k_z);

				if (k>0) {
					//What bin are we in? <That rhymes ##>
					int bin = floor((bins - 1) * (log(k) - log(min_k)) / (log(max_k) - log(min_k)));
					if (bin >= 0 && bin < bins) {
						double a1 = k_box1[half_box_wrap_idx(N, x, y, z)][0];
						double b1 = k_box1[half_box_wrap_idx(N, x, y, z)][1];

						double a2 = k_box2[half_box_wrap_idx(N, x, y, z)][0];
						double b2 = k_box2[half_box_wrap_idx(N, x, y, z)][1];

						double power = a1*a2 + b1*b2;

						//All except the z=0 and (if N is even) the z=N/2 planes count for double
						int multiplicity;
						if (z==0 || z==half_width) {
							multiplicity = 1;
						} else {
							multiplicity = 2;
						}

						//Deconvolution is optional and depends on the mass assignment method
						if (deconvolve == CLOUD_IN_CELL) {
							//The CIC Window function in Fourier space
							double W_x = (k_x == 0) ? 1 : pow(sinc(0.5*k_x*box_len/N),2);
							double W_y = (k_y == 0) ? 1 : pow(sinc(0.5*k_y*box_len/N),2);
							double W_z = (k_z == 0) ? 1 : pow(sinc(0.5*k_z*box_len/N),2);
							double W = W_x*W_y*W_z;
							double WW = abs(W)*abs(W);

							power /= WW;
						}

						/**
							double snc_x = 1 - (2./3.)*pow(sin(0.5*k_z*box_len/N),2);
							double snc_y = 1 - (2./3.)*pow(sin(0.5*k_y*box_len/N),2);
							double snc_z = 1 - (2./3.)*pow(sin(0.5*k_z*box_len/N),2);
							double shot_noise_correction = (snc_x*snc_y*snc_z)/pow(N,3);
							power -= shot_noise_correction;
						**/


						k_in_bins[bin] += multiplicity * k;
						power_in_bins[bin] += multiplicity * power;
						obs_in_bins[bin] += multiplicity;
					} else {
						//Should not happen
						throw std::logic_error("Incorrect bins in power spectrum calculation.");
					}
				}
			}
		}
	}

	for (int i=0; i<bins; i++) {
		k_in_bins[i] /= obs_in_bins[i];
		power_in_bins[i] /= obs_in_bins[i];
		power_in_bins[i] /= box_volume;
	}

}


/**
 * Function calc_powerspec()
 * Calculate the power spectrum of a field A, defined by
 * P_A(k) * (2pi)^3 * delta^3(k+k') = <A(k)A(k')^*>. See fourier.pdf.
 *
 * @param width The dimension of the grid (=N)
 * @param box_len Physical dimension of the box in Mpc
 * @param k_box (Field A) Reference to an array of N*N*(N/2+1) complex numbers
 * @param bins The number of logarithmic k bins
 * @param k_in_bins Reference to an array containing the average value of k in each bin
 * @param power_in_bins Reference to an array containing the average power of k in each bin
 * @param obs_in_bins Reference to an array containing the number of observations in each bin
 * @param deconvolve (optional) The window function that should be decolved (0=no deconvolution; default)
 */
void calc_powerspec(int width, double box_len, fftw_complex *k_box, int bins, double *k_in_bins, double *power_in_bins, int *obs_in_bins, int deconvolve=0) {
	//We just calculate the cross power spectrum of the field with itself
	calc_cross_powerspec(width, box_len, k_box, k_box, bins, k_in_bins, power_in_bins, obs_in_bins, deconvolve);
}
