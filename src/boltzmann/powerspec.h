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
 * powerspec.h  -  allows you to calculate the power spectrum from a grid.
 */

#ifndef POWERSPEC_H
#define POWERSPEC_H

/* Config parameters. */
#include "../config.h"

/* error() */
#include "../../logger/logger_header.h"

#ifdef HAVE_FFTW
#include <fftw3.h>
#endif

//Possible convolution windows
#define NO_DECONVOLUTION		0
#define CLOUD_IN_CELL			1
#define TRIANGULAR_CLOUD		2

/**
 * @brief Returns index of a half-complex N*N*(N/2+1) array using row-major style.
 *
 * Wraps around in the corresponding dimension if any of the 3 indices are out of
 * bounds.
 *
 * @param i Index along x.
 * @param j Index along y.
 * @param k Index along z.
 * @param N Size of the array along one of the larger axes.
 */
__attribute__((always_inline)) INLINE static int row_major_id_half_periodic(int i,
                                                                       int j,
                                                                       int k,
                                                                       int N) {
    int x = (i >= 0) ? i%N : i + ceil(-1.*i/N) * N;
	int y = (j >= 0) ? j%N : j + ceil(-1.*j/N) * N;
	int z = (k >= 0) ? k%(N/2+1) : k + ceil(-1.*k/(N/2+1)) * (N/2+1);
	return z + (N/2+1) * (y + N*x);
}

/**
 * @brief Inline sinc(x) = sin(x)/x function for convenience
 *
 * @param x Argument
 */
__attribute__((always_inline)) INLINE static double sinc(double x) {
    return (x == 0.0) ? 1.0 : sin(x)/x;
}

void calc_cross_powerspec(int width, double box_len, fftw_complex *k_box1,
                        fftw_complex *k_box2, int bins, double *k_in_bins,
                        double *power_in_bins, int *obs_in_bins, int deconvolve);


#endif /* POWERSPEC_H */
