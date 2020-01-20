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


double fermi_dirac_density(const struct engine *engine, double* x, float* v) {
    const struct cosmology *cosmo = engine->cosmology;
    const struct phys_const *physical_constants = engine->physical_constants;

    const double T_nu = cosmo->T_nu;
    const double k_b = physical_constants->const_boltzmann_k;
    const double eV = physical_constants->const_electron_volt;
    const double T_eV = k_b*T_nu/eV; // temperature in eV

    //Calculate momentum in eV
    double p = fermi_dirac_momentum(engine, v);

    double norm = 1.16748e+11;

    return norm * p*p / (exp(p / T_eV) + 1.0);
}

double sample_density(const struct engine *engine, double* x, float* v) {
    const struct cosmology *cosmo = engine->cosmology;
    const struct phys_const *physical_constants = engine->physical_constants;

    const double T_nu = cosmo->T_nu;
    const double k_b = physical_constants->const_boltzmann_k;
    const double eV = physical_constants->const_electron_volt;
    const double T_eV = k_b*T_nu/eV; // temperature in eV

    //Calculate momentum in eV
    double p = fermi_dirac_momentum(engine, v);

    double norm = 8573.24;

    return norm * 1.0 / (exp(p / T_eV) + 1.0);
}

/* Calculate the momentum in eV, using E = a*sqrt(p^2 + m^2) ~ ap. */
double fermi_dirac_momentum(const struct engine *engine, float* v) {
    const struct cosmology *cosmo = engine->cosmology;
    const struct phys_const *physical_constants = engine->physical_constants;

    //Some constants
    const double c = physical_constants->const_speed_light_c;
    const double cc = c*c;
    const double a = cosmo->a;
    const double aa = a*a;
    const double eV = physical_constants->const_electron_volt;
    const double eV_mass = eV/cc; // 1 eV/c^2 in internal mass units

    //Calculate the neutrino mass in internal units
    const double M_nu = cosmo->M_nu[0] * eV_mass; // just select the first species for now

    //The internal velocity V = a^2*(dx/dt), where x=r/a is comoving
    double VV = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    double V = sqrt(VV);

    //Calculate the length of the physical 3-velocity u=a*|dx/dt|
    double u = a*c*V/sqrt(aa*cc + aa*VV - VV);
    double gamma = 1.0/sqrt(1.0 - u*u/cc); //Lorentz factor
    double p_ph = u*gamma*M_nu; //The physical 3-momentum
    double p_eV = p_ph * c/eV; //in eV

    return p_eV;
}

/* Calculate the energy in units of M_nu */
double fermi_dirac_energy(const struct engine *engine, float* v) {
    const struct cosmology *cosmo = engine->cosmology;

    //Calculate the energy in eV
    double M_nu = cosmo->M_nu[0]; // just select the first species for now
    double p_eV = fermi_dirac_momentum(engine, v);
    double E_eV = sqrt(p_eV*p_eV + M_nu*M_nu);

    return E_eV/M_nu; //energy in units of M_nu
}
