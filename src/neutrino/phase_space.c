/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Willem Elbers (whe@willemelbers.com)
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
 *  @file phase_space.c
 *  @brief Functions relating to the phase space density function, for use
 *  by pseudo-particles (neutrinos).
 */

/* This object's header. */
#include "phase_space.h"

/* The general cosmology header */
#include "../cosmology.h"

/* Some standard headers */
#include <math.h>

/* Calculation of sqrt(x^2 + y^2 + z^2), without undue overflow or underflow. */
__attribute__((always_inline, const)) INLINE static double hypot3(double x,
                                                                  double y,
                                                                  double z) {

  return hypot(x, hypot(y, z));
}

double fermi_dirac_density(const struct engine *e, float *v, double m_eV,
                           double T_factor) {
  const struct phys_const *physical_constants = e->physical_constants;

/* Retrieve the neutrino temperature today */
#ifdef NEUTRINO_BACKGROUND
  const struct cosmology *cosmo = e->cosmology;
  const double T_nu = cosmo->T_nu;
#else
  /* No neutrino cosmology module, use the fiducial value */
  const struct unit_system *internal_units = e->internal_units;
  const double T_nu_K = NEUTRINO_FIDUCIAL_TEMPERATURE_KELVIN;
  const double T_nu = T_nu_K * internal_units->UnitTemperature_in_cgs;
#endif

  /* Convert temperature to eV (to prevent overflows)*/
  const double k_b = physical_constants->const_boltzmann_k;
  const double eV = physical_constants->const_electron_volt;
  const double T_eV = k_b * T_nu / eV;  // temperature in eV

  /* Calculate the momentum in eV */
  double p_eV = fermi_dirac_momentum(e, v, m_eV);

  return 1.0 / (exp(p_eV / (T_eV * T_factor)) + 1.0);
}

/* Calculate the momentum in energy units, using E = a*sqrt(p^2 + m^2) ~ ap.
 * Note that this is the present-day momentum, i.e. p0 = ap, which is
 * constant in a homogenous Universe.
 */
double fermi_dirac_momentum(const struct engine *e, float *v, double m_eV) {
  const struct cosmology *cosmo = e->cosmology;
  const struct phys_const *physical_constants = e->physical_constants;
  const double c = physical_constants->const_speed_light_c;
  const double a = cosmo->a;

  // The internal velocity V = a^2*(dx/dt), where x=r/a is comoving
  double V = hypot3(v[0], v[1], v[2]);

  // Calculate the length of the physical 3-velocity u=a*|dx/dt|
  double u = V / a;

  // If we use relativistic equations of motion, then V already
  // includes the gamma factor. If we don't use relativistic equations
  // of motion, gamma is assumed to be 1.
  double gamma = 1.0;

  // Convert to present-day momentum in eV
  double p_eV = u * gamma * m_eV / c;  // The physical 3-momentum in eV
  double p0_eV = p_eV * a;             // present-day momentum in eV

  return p0_eV;
}

/* Compute the energy E=sqrt(p^2+m^2) of a micrscopic neutrino in eV */
double fermi_dirac_energy(const struct engine *e, float vi, double m_eV) {
  const struct cosmology *cosmo = e->cosmology;
  const struct phys_const *physical_constants = e->physical_constants;
  const double c = physical_constants->const_speed_light_c;
  const double a = cosmo->a;

  /* Compute the momentum and energy at scale-factor a in eV */
  double p_eV = vi * m_eV / (a * c);
  double energy_eV = hypot(m_eV, p_eV);

  return energy_eV;
}

/* Compute the ratio of macro particle mass in internal mass units to
 * the mass of one microscopic neutrino in eV.
 */
double neutrino_mass_factor(const struct engine *e) {
  const struct phys_const *physical_constants = e->physical_constants;
  const struct space *s = e->s;

  /* Some constants */
  const double k_b = physical_constants->const_boltzmann_k;
  const double hbar = physical_constants->const_planck_hbar;
  const double c = physical_constants->const_speed_light_c;
  const double eV = physical_constants->const_electron_volt;
  const double eV_mass = eV / (c * c);  // 1 eV/c^2 in internal mass units
  const double prefactor = (1.5 * M_ZETA_3) / (M_PI * M_PI);

/* Retrieve the neutrino temperature today & number of flavours */
#ifdef NEUTRINO_BACKGROUND
  const struct cosmology *cosmo = e->cosmology;
  const double T_nu = cosmo->T_nu;
  const double flavours = cosmo->N_nu;
#else
  /* No neutrino cosmology module, use the fiducial values */
  const struct unit_system *internal_units = e->internal_units;
  const double T_nu_K = NEUTRINO_FIDUCIAL_TEMPERATURE_KELVIN;
  const double T_nu = T_nu_K * internal_units->UnitTemperature_in_cgs;
  const double flavours = NEUTRINO_FIDUCIAL_FLAVOURS_NUMBER;
#endif

  /* Compute the comoving number density per flavour */
  const double n = prefactor * pow(k_b * T_nu / (hbar * c), 3);
  const double volume = s->dim[0] * s->dim[1] * s->dim[2];

  /* The total number of neutrino macropaticles present */
  const long long nuparts = e->total_nr_nuparts;

  /* Compute the conversion factor */
  const double mass_factor = nuparts / (flavours * n * volume);

  /* Convert to eV */
  const double mass_factor_eV = mass_factor / eV_mass;

  return mass_factor_eV;
}
