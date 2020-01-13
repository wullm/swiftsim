/**
 * @infini_lpt.h
 * @author  Willem Elbers <whe@willemelbers.com>
 * @version 1.0
 * @license MIT (https://opensource.org/licenses/MIT)
 *
 * @section DESCRIPTION
 *
 * Infinite Lagrangian perturbation theory.
 */

#ifndef INFINI_LPT_H
#define INFINI_LPT_H

#define MONGE_AMPERE_MAX_ITER 100
#define MONGE_AMPERE_VERBOSE 1

#include <fftw3.h>
#include <random>
#include <stdexcept>

#include "config.h"
#include "field_io.h"

//Linear algebra
#include <Eigen/Dense>
typedef Eigen::Matrix<double, 3, 3> Matrix3d;

void do_infini_lpt(double*, double*, int, double, double tol = 1e-14);

#endif
