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

double fermi_dirac_density(const struct cosmology *cosmo, double* x, float* v) {
    double T_nu = 1.952784; //K
    double k_b = 8.617333262145e-5; //eV/K

    //Convert to eV
    double p = fermi_dirac_momentum(cosmo, x, v);

    return p*p/(exp(p/(k_b*T_nu)) + 1.0);
}

double fermi_dirac_momentum(const struct cosmology *cosmo, double* x, float* v) {
    double a = cosmo->a;
    double aa = a*a;

    double VV = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    double V = sqrt(VV);

    //Convert to Mpc/Gyr
    V *= 1.0;
    VV = V*V;

    //Some constants
    double c = 306.4;
    double cc = c*c;
    double M_nu = 0.2; //eV
    double eV_mass = 1.782662e-36;  // kg (this is eV/c^2)
    double M_nu_kg = M_nu*eV_mass; //kg
    double Mpc = 3.086e22;  // m
    double Gyr = 3.154e16;  // s
    double eV = 1.602176634e-19;    // J

    //Relativistic correction, u=(dx/dt) with x=r/a comoving coordinates and t coordinate time
    double u = a*c*V/sqrt(aa*cc + aa*VV - VV);

    //The Lorentz factor
    double gamma = 1.0/sqrt(1.0 - u*u/cc);

    //Get the physical momentum
    double p_ph = u*gamma*M_nu_kg;

    //Convert to eV
    double p = p_ph * c/eV * (Mpc/Gyr)*(Mpc/Gyr);

    return p;
}
