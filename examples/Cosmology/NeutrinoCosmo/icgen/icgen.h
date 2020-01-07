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
std::vector<double> TF_T_rho_cdm;
std::vector<double> TF_T_rho_nu;
std::vector<double> TF_T_rho_b;
std::vector<double> TF_T_rho_cb; //cdm and baryon weighted mean

//Velocity transfer functions
std::vector<double> TF_T_theta_nu;
std::vector<double> TF_T_theta_b;
std::vector<double> TF_T_theta_cb; //cdm and baryon weighted mean

// Indexed search table for the transfer function interpolation
const int TF_I_max = 100;
float *TF_index;

// The smallest and largest k-values
float log_k_min, log_k_max;

// Linearly interpolate the Transfer functon
inline double Transfer_interpol(double k, std::vector<double> *Transfer) {
  if (k > TF_ks[TF_ks.size() - 1]) {
    return (*Transfer)[TF_ks.size() - 1];
  } else if (k < TF_ks[0]) {
    return (*Transfer)[0];
  } else {
    double k1, k2, T1, T2;

    // Quickly find a starting index using the indexed seach
    float v = log(k);
    float u = (v - log_k_min) / (log_k_max - log_k_min);
    int I = floor(u * TF_I_max);
    int idx = TF_index[I < TF_I_max ? I : TF_I_max - 1];

    // Search in the TF table
    for (int i = idx; i < TF_ks.size(); i++) {
      if (k >= TF_ks[i] && i < TF_ks.size() - 1) {
        if (k <= TF_ks[i + 1]) {
          k1 = TF_ks[i];
          k2 = TF_ks[i + 1];
          T1 = (*Transfer)[i];
          T2 = (*Transfer)[i + 1];
        }
      }
    }

    // Linear interpolation
    double z = (k - k1) / (k2 - k1);
    double T = T1 + z * (T2 - T1);

    return T;
  }
}

//Density sigma functions
inline double sigma_func_cdm(double k) { //this gives cdm+baryons
  return sqrt(pow(k, N_S)) * Transfer_interpol(k, &TF_T_rho_cb);
}

inline double sigma_func_neutrino(double k) {
  return sqrt(pow(k, N_S)) * Transfer_interpol(k, &TF_T_rho_nu);
}

//Velocity sigma functions
inline double sigma_func_vel_cdm(double k) { //this gives cdm+baryons
  return sqrt(pow(k, N_S)) * Transfer_interpol(k, &TF_T_theta_cb);
}

inline double sigma_func_vel_neutrino(double k) {
  return sqrt(pow(k, N_S)) * Transfer_interpol(k, &TF_T_theta_nu);
}

void writeGRF_H5(double *box, size_t N, float box_len, std::string fname);

#endif
