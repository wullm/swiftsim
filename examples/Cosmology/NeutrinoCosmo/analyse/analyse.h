/**
 * @analyse.h
 * @author  Willem Elbers <whe@willemelbers.com>
 * @version 1.0
 * @license MIT (https://opensource.org/licenses/MIT)
 *
 * @section DESCRIPTION
 *
 * Analyse the output of cosmological N-body simulation.
 */

#ifndef ANALYSE_H
#define ANALYSE_H

struct corpuscle {
  long int id;
  double X, Y, Z; //position
  double v_X, v_Y, v_Z; //velocity
  double V; //length of velocity vector

  //initial velocities
  double ic_v_X, ic_v_Y, ic_v_Z;
  double ic_V;


  double delta_X, delta_Y, delta_Z; //displacement from grid/glass position
  double mass;
  double smoothing_length;
  double internal_energy;
};

#include <math.h>
#include <chrono>
#include <iostream>
#include <random>
#include <vector>


#include <fftw3.h>

#endif
