/**
 * @sampler.h
 * @author  Willem Elbers <whe@willemelbers.com>
 * @version 1.0
 * @license MIT (https://opensource.org/licenses/MIT)
 *
 * @section DESCRIPTION
 *
 * Generates random numbers according to a specified probability density
 * function.
 */

#ifndef SAMPLER_H
#define SAMPLER_H

#include <random>
#include <vector>

// Linked list of intervals
struct interval {
  int id;
  double l, r;            // endpoints left and right
  double Fl, Fr;          // cdf evaluations at endpoints
  double a0, a1, a2, a3;  // cubic Hermite coefficients
  double error;           // error at midpoint
  int nid;                // the next interval
};

inline bool compareByLeft(const interval &a, const interval &b) {
  return a.Fl < b.Fl;
}
inline bool compareByRight(const interval &a, const interval &b) {
  return a.Fr < b.Fr;
}

struct sampler {
  // The normalization of the pdf (will be updated later)
  double norm = 1;

  // The random number generator
  std::default_random_engine oracle;
  std::uniform_real_distribution<double> Uniform;

  // The intervals used in the interpolation
  std::vector<interval> intervals;

  // The number of intervals
  int intervalNum;

  // Make an indexed search table
  const int I_max = 5000;
  float *index;
};

void seed_rng(struct sampler *s, int seed);

void prepare_intervals(struct sampler *s, double T, double mu);

double sampler_draw(struct sampler *s);

#endif
