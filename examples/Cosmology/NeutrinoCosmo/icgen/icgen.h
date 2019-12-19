/**
 * @icgen.h
 * @author  Willem Elbers <whe@willemelbers.com>
 * @version 1.0
 * @license MIT (https://opensource.org/licenses/MIT)
 *
 * @section DESCRIPTION
 *
 * Generate initial conditions for a cosmological N-body simulation.
 */

#ifndef ICGEN_H
#define ICGEN_H

struct corpuscle {
  long int id;
  double X, Y, Z;
  double v_X, v_Y, v_Z;
  double mass;
  double smoothing_length;
  double internal_energy;
};

#include <math.h>
#include <chrono>
#include <iostream>
#include <random>
#include <vector>

// #include <sstream>
// #include <fstream>
// #include <algorithm>
// #include <cstring>
//
// //Linear algebra
// #include <Eigen/Dense>
//
#include <fftw3.h>

// Tables of transfer functions
std::vector<double> TF_ks;
std::vector<double> TF_T_rho;
std::vector<double> TF_T_rho_nu;

// Linearly interpolate the Transfer functon
inline double Transfer_interpol(double k, std::vector<double> Transfer) {
  if (k > TF_ks[TF_ks.size() - 1]) {
    return Transfer[TF_ks.size() - 1];
  } else if (k < TF_ks[0]) {
    return Transfer[0];
  } else {
    double k1, k2, T1, T2;

    for (int i = 0; i < TF_ks.size(); i++) {
      if (k >= TF_ks[i] && i < TF_ks.size() - 1) {
        if (k <= TF_ks[i + 1]) {
          k1 = TF_ks[i];
          k2 = TF_ks[i + 1];
          T1 = Transfer[i];
          T2 = Transfer[i + 1];
        }
      }
    }

    double u = (k - k1) / (k2 - k1);
    double T = T1 + u * (T2 - T1);

    return T;
  }
}

// Doesn't have to be normalized yet
// inline double sigma_func_no_transfer(double k) {
//   return sqrt(pow(k, 0.97));
// }

// Doesn't have to be normalized yet
inline double sigma_func_cdm(double k) {
  return sqrt(pow(k, 0.97)) * Transfer_interpol(k, TF_T_rho);
}

// Doesn't have to be normalized yet
inline double sigma_func_neutrino(double k) {
  return sqrt(pow(k, 0.97)) * Transfer_interpol(k, TF_T_rho_nu);
}

void writeGRF_H5(double *box, size_t N, float box_len, std::string fname);

#endif
