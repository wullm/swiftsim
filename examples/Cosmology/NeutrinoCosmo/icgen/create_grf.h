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

#ifndef CREATE_GRF_H
#define CREATE_GRF_H

#include <random>
#include <stdexcept>
#include <fftw3.h>

#include "field_io.h"
#include "config.h"


void generate_grf(std::default_random_engine&, fftw_complex*, int, double, double (&f)(double));
double integrate_sigma_R(int, double, double, double (&f)(double));
long int sayhello();


#endif
