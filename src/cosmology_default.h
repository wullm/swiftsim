/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_COSMOLOGY_DEFAULT_H
#define SWIFT_COSMOLOGY_DEFAULT_H

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

  /*! Neutrino density parameter */
  double Omega_nu;

  /*! Photon density parameter */
  double Omega_g;

  /*! Baryon density parameter */
  double Omega_b;

  /*! Cosmological constant density parameter */
  double Omega_lambda;

  /*! Radiation density parameter */
  double Omega_r;

  /*! Curvature density parameter */
  double Omega_k;

  /*! Dark-energy equation of state at z=0 */
  double w_0;

  /*! Dark-energy evolution parameter */
  double w_a;

  /*! Log of starting expansion factor */
  double log_a_begin;

  /*! Log of final expansion factor */
  double log_a_end;

  /*! Drift factor interpolation table */
  double *drift_fac_interp_table;

  /*! Kick factor (gravity) interpolation table */
  double *grav_kick_fac_interp_table;

  /*! Kick factor (hydro) interpolation table */
  double *hydro_kick_fac_interp_table;

  /*! Kick factor (hydro correction) interpolation table (GIZMO-MFV only) */
  double *hydro_kick_corr_interp_table;

  /*! Time interpolation table */
  double *time_interp_table;

  /*! Scale factor interpolation table */
  double *scale_factor_interp_table;

  /*! Time between Big Bang and first entry in the table */
  double time_interp_table_offset;

  /*! Time at the present-day (a=1) */
  double universe_age_at_present_day;
};

#endif /* SWIFT_COSMOLOGY_DEFAULT_H */
