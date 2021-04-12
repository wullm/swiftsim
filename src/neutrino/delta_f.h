/**
 *  Copyright (c) 2021 Willem Elbers (whe@willemelbers.com)
 *
 *  @file phase_space.h
 *  @brief Computes weights for neutrino particles according to
 *  the delta-f method (arXiv:2010.07321).
 */

/* Standard headers */
#include <math.h>
#include <stdint.h>

/* Lookup table for Fermi-Dirac transform */
#include "fermi_dirac.h"

static inline double hypot32(double x, double y, double z) {
    return hypot(x, hypot(y, z));
}

/* Pseudo-random number generator */
static inline uint64_t splitmix64(uint64_t *state) {
    uint64_t result = *state;

    *state = result + 0x9E3779B97f4A7C15;
    result = (result ^ (result >> 30)) * 0xBF58476D1CE4E5B9;
    result = (result ^ (result >> 27)) * 0x94D049BB133111EB;
    return result ^ (result >> 31);
}

/* Generate a uniform variable on the open unit interval */
static inline double sampleUniform(uint64_t *state) {
    const uint64_t A = splitmix64(state);
    const double RM = (double)UINT64_MAX + 1;
    return ((double)A + 0.5) / RM;
}

/**
 * @brief Calculate the present-day momentum in electronvolts. We can use
 * energy units because E = a*sqrt(p^2 + m^2) = a*p, since p/m >> 1 at
 * decoupling. Note that this quantity is constant in a homogeneous Universe.
 *
 * @param v Array of 3 velocity components (a^2 dx/dt, where x is comoving)
 * @param m_eV Neutrino mass in electronvolts
 * @param c Speed of light
 */
static inline double fermi_dirac_momentum2(double *v, double m_eV, double c) {
    const double u = hypot32(v[0], v[1], v[2]);
    const double p_eV = u * m_eV / c;

    return p_eV;
}

/**
 * @brief Calculate the neutrino density at the particle's location in phase
 * space, according to the 0th order background model: f_0(x,p,t).
 *
 * @param p_eV Present-day momentum in electronvolts
 * @param T_eV Present-day temperature in electronvolts
 */
static inline double fermi_dirac_density2(double p_eV, double T_eV) {
    return 1.0 / (exp(p_eV / T_eV) + 1.0);
}
