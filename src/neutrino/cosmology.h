/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 * Copyright (c) 2019 Willem Elbers (whe@willemelbers.com)
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
#ifndef SWIFT_NEUTRINO_COSMOLOGY_H
#define SWIFT_NEUTRINO_COSMOLOGY_H

/**
 * @brief Cosmological parameters
 */
struct cosmology {

  /*! Current expansion factor of the Universe */
  double a;

  /*! Inverse of the current expansion factor of the Universe */
  double a_inv;

  /*! Inverse square of the current expansion factor of the Universe */
  double a2_inv;

  /*! Inverse cube of the current expansion factor of the Universe */
  double a3_inv;

  /*! Power of the scale-factor used for internal energy conversion to physical
   */
  double a_factor_internal_energy;

  /*! Power of the scale-factor used for pressure conversion to physical */
  double a_factor_pressure;

  /*! Power of the scale-factor used for sound-speed conversion to physical */
  double a_factor_sound_speed;

  /*! Power of the scale-factor used for relative velocities in visc. terms */
  double a_factor_mu;

  /*! Power of the scale-factor used for epsilon term in the Balsara switch */
  double a_factor_Balsara_eps;

  /*! Power of the scale-factor used for gravity accelerations */
  double a_factor_grav_accel;

  /*! Power of the scale-factor used for hydro accelerations */
  double a_factor_hydro_accel;

  /*! Current redshift */
  double z;

  /*! Hubble constant at the current redshift (in internal units) */
  double H;

  /*! The critical density at the current redshift (in internal physical units)
   */
  double critical_density;

  /*! The critical density at redshift 0 (in internal physical units) */
  double critical_density_0;

  /*! The mean density at the current redshift (in internal physical units) */
  double mean_density;

  /*! The mean baryonic density at the current redshift (in internal physical
   * units) */
  double mean_density_Omega_b;

  /*! Conversion factor from internal time-step size to cosmological step */
  double time_step_factor;

  /*! Expansion rate at the current redshift (in internal units) */
  double a_dot;

  /*! Time (in internal units) since the Big Bang */
  double time;

  /*! Conformal time (in internal units) since the Big Bang */
  double conformal_time;

  /*! Lookback time (in internal units) */
  double lookback_time;

  /*! Dark-energy equation of state at the current time */
  double w;

  /*! Scale-factor at the previous time-step */
  double a_old;

  /*! Redshit at the previous time-step */
  double z_old;

  /*------------------------------------------------------------------ */

  /*! Starting expansion factor */
  double a_begin;

  /*! Final expansion factor */
  double a_end;

  /*! Time (in internal units) since the Big Bang at the start */
  double time_begin;

  /*! Time (in internal units) since the Big Bang at the end */
  double time_end;

  /*! Conversion factor from integer time-line to \f$ d\log{a} \f$ */
  double time_base;

  /*! Inverse of conversion factor from integer time-line to \f$ d\log{a} \f$ */
  double time_base_inv;

  /*! Reduced Hubble constant (H0 / (100km/s/Mpc)) */
  double h;

  /*! Hubble constant at z = 0 (in internal units) */
  double H0;

  /*! Hubble time 1/H0 */
  double Hubble_time;

  /*! Matter density parameter */
  double Omega_m;

  /*! Baryon density parameter */
  double Omega_b;

  /*! Cosmological constant density parameter */
  double Omega_lambda;

  /*! Radiation density parameter */
  double Omega_r;

  /*! Photon density parameter */
  double Omega_g;

  /*! Curvature density parameter */
  double Omega_k;

  /*! Number of neutrino species */
  size_t N_nu;

  /*! Neutrino masses in eV */
  double *M_nu;

  /*! The total neutrino mass in eV */
  double M_nu_tot;

  /*! Total neutrino density parameter */
  double Omega_nu;

  /*! Effective number of relativistic species at early times */
  double N_eff;

  /*! CMB temperature today */
  double T_CMB;

  /*! Neutrino temperature today */
  double T_nu;

  /*! Dark-energy equation of state at z=0 */
  double w_0;

  /*! Dark-energy evolution parameter */
  double w_a;

  /*! Log of starting expansion factor */
  double log_a_begin;

  /*! Log of final expansion factor */
  double log_a_end;

  /*! Log of starting expansion factor of the interpolation tables */
  double log_a_table_begin;

  /*! Log of final expansion factor of the interpolation tables */
  double log_a_table_end;

  /*! Scale-factor interpolation table */
  double *log_a_interp_table;

  /*! Drift factor interpolation table */
  double *drift_fac_interp_table;

  /*! Kick factor (gravity) interpolation table */
  double *grav_kick_fac_interp_table;

  /*! Kick factor (hydro) interpolation table */
  double *hydro_kick_fac_interp_table;

  /*! Kick factor (hydro correction) interpolation table (GIZMO-MFV only) */
  double *hydro_kick_corr_interp_table;

  /*! Time since Big Bang interpolation table */
  double *time_interp_table;

  /*! Conformal time since Big Bang interpolation table */
  double *conformal_time_interp_table;

  /*! Neutrino density interpolation table */
  double *neutrino_density_interp_table;

  /*! Log of scale factor at which the neutrino table begins */
  double log_a_nutab_begin;

  /*! Log of scale factor at which the neutrino table ends */
  double log_a_nutab_end;

  /*! Time at the present-day (a=1) */
  double universe_age_at_present_day;
};

double cosmology_get_neutrino_density_param(const struct cosmology *c,
                                            double a);

void cosmology_neutrino_init(struct swift_params *params,
                             const struct unit_system *us,
                             const struct phys_const *phys_const,
                             struct cosmology *c);

#endif /* SWIFT_NEUTRINO_COSMOLOGY_H */
