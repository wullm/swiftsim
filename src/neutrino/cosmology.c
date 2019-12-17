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

/**
 *  @file cosmology.c
 *  @brief Functions relating cosmological parameters, specific to the neutrino
 *  model
 */

/* The general cosmology header (this object's header is included there). */
#include "../cosmology.h"

/* Some standard headers */
#include <math.h>

/* Local headers */
#include "../adiabatic_index.h"
#include "../align.h"
#include "../common_io.h"
#include "../inline.h"
#include "../memuse.h"
#include "../minmax.h"
#include "../restart.h"

#ifdef HAVE_LIBGSL
#include <../gsl/gsl_integration.h>
#endif

/*! Number of values stored in the cosmological interpolation tables */
const int cosmology_table_length = 10000;

/*! Number of values stored in longer cosmological tables */
const int cosmology_long_table_length = 100000;

#ifdef HAVE_LIBGSL
/*! Size of the GSL workspace */
const size_t GSL_workspace_size = 100000;
#endif

/**
 * @brief Returns the interpolated value from a table.
 *
 * Uses linear interpolation.
 *
 * @brief table The table of value to interpolate from (should be of length
 * cosmology_table_length).
 * @brief x The value to interpolate at.
 * @brief x_min The mininum of the range of x.
 * @brief x_max The maximum of the range of x.
 */
static INLINE double interp_table(const double *table, const double x,
                                  const double x_min, const double x_max) {

  const double xx =
      ((x - x_min) / (x_max - x_min)) * ((double)cosmology_table_length);

  const int i = (int)xx;
  const int ii = min(cosmology_table_length - 1, i);

  /* Indicate that the whole array is aligned on boundaries */
  swift_align_information(double, table, SWIFT_STRUCT_ALIGNMENT);

  if (ii <= 1)
    return table[0] * xx;
  else
    return table[ii - 1] + (table[ii] - table[ii - 1]) * (xx - ii);
}

/**
 * @brief Returns the interpolated value from a longer cosmological table.
 *
 * Uses linear interpolation. Slightly different from interp_table in the
 * way it handles the lower bound of the table.
 *
 * @brief table The table of value to interpolate from (should be of length
 * cosmology_long_table_length).
 * @brief x The value to interpolate at.
 * @brief x_min The mininum of the range of x.
 * @brief x_max The maximum of the range of x.
 */
static INLINE double interp_long_table(const double *table, const double x,
                                       const double x_min, const double x_max) {

  const double xx =
      ((x - x_min) / (x_max - x_min)) * ((double)cosmology_long_table_length);

  const int i = (int)xx;
  const int ii = min(cosmology_long_table_length - 1, i);

  /* Indicate that the whole array is aligned on boundaries */
  swift_align_information(double, table, SWIFT_STRUCT_ALIGNMENT);

  if (ii < 1)
    return table[0];
  else
    return table[ii - 1] + (table[ii] - table[ii - 1]) * (xx - ii);
}

/**
 * @brief Computes the dark-energy equation of state at a given scale-factor a.
 *
 * We follow the convention of Linder & Jenkins, MNRAS, 346, 573, 2003
 *
 * @param a The current scale-factor
 * @param w_0 The equation of state parameter at z=0
 * @param w_a The equation of state evolution parameter
 */
__attribute__((const)) static INLINE double cosmology_dark_energy_EoS(
    const double a, const double w_0, const double w_a) {

  return w_0 + w_a * (1. - a);
}

/**
 * @brief Computes the integral of the dark-energy equation of state
 * up to a scale-factor a.
 *
 * We follow the convention of Linder & Jenkins, MNRAS, 346, 573, 2003
 * and compute \f$ \tilde{w}(a) = \int_0^a\frac{1 + w(z)}{1+z}dz \f$.
 *
 * @param a The current scale-factor.
 * @param w0 The equation of state parameter at z=0.
 * @param wa The equation of state evolution parameter.
 */
__attribute__((const)) static INLINE double w_tilde(const double a,
                                                    const double w0,
                                                    const double wa) {
  return (a - 1.) * wa - (1. + w0 + wa) * log(a);
}

/**
 * @brief Compute \f$ E(z) \f$.
 *
 * @param Omega_r The radiation density parameter \f$ \Omega_r \f$.
 * @param Omega_m The matter density parameter \f$ \Omega_m \f$.
 * @param Omega_k The curvature density parameter \f$ \Omega_k \f$.
 * @param Omega_l The cosmological constant density parameter \f$ \Omega_\Lambda
 * \f$.
 * @param w0 The equation of state parameter at z=0.
 * @param wa The equation of state evolution parameter.
 * @param a The current scale-factor.
 */
__attribute__((const)) static INLINE double E(
    const double Omega_r, const double Omega_m, const double Omega_k,
    const double Omega_l, const double w0, const double wa, const double a) {

  const double a_inv = 1. / a;

  return sqrt(Omega_r * a_inv * a_inv * a_inv * a_inv + /* Radiation */
              Omega_m * a_inv * a_inv * a_inv +         /* Matter */
              Omega_k * a_inv * a_inv +                 /* Curvature */
              Omega_l * exp(3. * w_tilde(a, w0, wa)));  /* Lambda */
}

/**
 * @brief Returns the time (in internal units) since Big Bang at a given
 * scale-factor.
 *
 * @param c The current #cosmology.
 * @param a Scale-factor of interest.
 */
double cosmology_get_time_since_big_bang(const struct cosmology *c, double a) {

#ifdef SWIFT_DEBUG_CHECKS
  if (a < c->a_begin) error("Error a can't be smaller than a_begin");
#endif

  /* Time between a_begin and a */
  const double delta_t =
      interp_table(c->time_interp_table, log(a), c->log_a_begin, c->log_a_end);

  return c->time_interp_table_offset + delta_t;
}

/**
 * @brief Update the cosmological parameters to the current simulation time.
 *
 * @param c The #cosmology struct.
 * @param phys_const The physical constants in the internal units.
 * @param ti_current The current (integer) time.
 */
void cosmology_update(struct cosmology *c, const struct phys_const *phys_const,
                      integertime_t ti_current) {

  /* Save the previous state */
  c->z_old = c->z;
  c->a_old = c->a;

  /* Get scale factor and powers of it */
  const double a = c->a_begin * exp(ti_current * c->time_base);
  const double a_inv = 1. / a;
  c->a = a;
  c->a_inv = a_inv;
  c->a2_inv = a_inv * a_inv;
  c->a3_inv = a_inv * a_inv * a_inv;
  c->a_factor_internal_energy =
      pow(a, -3. * hydro_gamma_minus_one);          /* a^{3*(1-gamma)} */
  c->a_factor_pressure = pow(a, -3. * hydro_gamma); /* a^{-3*gamma} */
  c->a_factor_sound_speed =
      pow(a, -1.5 * hydro_gamma_minus_one); /* a^{3*(1-gamma)/2} */
  c->a_factor_grav_accel = a_inv * a_inv;   /* 1 / a^2 */
  c->a_factor_hydro_accel =
      pow(a, -3. * hydro_gamma + 2.); /* 1 / a^(3*gamma - 2) */
  c->a_factor_mu =
      pow(a, 0.5 * (3. * hydro_gamma - 5.)); /* a^{(3*gamma - 5) / 2} */
  c->a_factor_Balsara_eps =
      pow(a, 0.5 * (1. - 3. * hydro_gamma)); /* a^{(1 - 3*gamma) / 2} */

  /* Redshift */
  c->z = a_inv - 1.;

  /* Dark-energy equation of state */
  c->w = cosmology_dark_energy_EoS(a, c->w_0, c->w_a);

  /* E(z) */
  const double Omega_g = c->Omega_g;
  const double Omega_nu = cosmology_get_neutrino_density_param(c, a);
  const double Omega_r = (Omega_g > 0) ? Omega_g + Omega_nu : c->Omega_r;
  const double Omega_m = c->Omega_m;
  const double Omega_k = c->Omega_k;
  const double Omega_l = c->Omega_lambda;
  const double w_0 = c->w_0;
  const double w_a = c->w_a;
  const double E_z = E(Omega_r, Omega_m, Omega_k, Omega_l, w_0, w_a, a);

  /* H(z) */
  c->H = c->H0 * E_z;

  /* Expansion rate */
  c->a_dot = c->H * c->a;

  /* Critical density */
  c->critical_density =
      3. * c->H * c->H / (8. * M_PI * phys_const->const_newton_G);

  /* Mean density */
  c->mean_density = c->critical_density_0 * c->a3_inv;

  /* Mean baryonic density */
  c->mean_density_Omega_b = c->mean_density * c->Omega_b;

  /* Time-step conversion factor */
  c->time_step_factor = c->H;

  /* Time */
  c->time = cosmology_get_time_since_big_bang(c, a);
  c->lookback_time = c->universe_age_at_present_day - c->time;
}

/**
 * @brief Computes \f$ dt / a^2 \f$ for the current cosmology.
 *
 * @param a The scale-factor of interest.
 * @param param The current #cosmology.
 */
double drift_integrand(double a, void *param) {

  const struct cosmology *c = (const struct cosmology *)param;
  const double Omega_g = c->Omega_g;
  const double Omega_nu = cosmology_get_neutrino_density_param(c, a);
  const double Omega_r = (Omega_g > 0) ? Omega_g + Omega_nu : c->Omega_r;
  const double Omega_m = c->Omega_m;
  const double Omega_k = c->Omega_k;
  const double Omega_l = c->Omega_lambda;
  const double w_0 = c->w_0;
  const double w_a = c->w_a;
  const double H0 = c->H0;

  const double a_inv = 1. / a;
  const double E_z = E(Omega_r, Omega_m, Omega_k, Omega_l, w_0, w_a, a);
  const double H = H0 * E_z;

  return (1. / H) * a_inv * a_inv * a_inv;
}

/**
 * @brief Computes \f$ dt / a \f$ for the current cosmology.
 *
 * @param a The scale-factor of interest.
 * @param param The current #cosmology.
 */
double gravity_kick_integrand(double a, void *param) {

  const struct cosmology *c = (const struct cosmology *)param;
  const double Omega_g = c->Omega_g;
  const double Omega_nu = cosmology_get_neutrino_density_param(c, a);
  const double Omega_r = (Omega_g > 0) ? Omega_g + Omega_nu : c->Omega_r;
  const double Omega_m = c->Omega_m;
  const double Omega_k = c->Omega_k;
  const double Omega_l = c->Omega_lambda;
  const double w_0 = c->w_0;
  const double w_a = c->w_a;
  const double H0 = c->H0;

  const double a_inv = 1. / a;
  const double E_z = E(Omega_r, Omega_m, Omega_k, Omega_l, w_0, w_a, a);
  const double H = H0 * E_z;

  return (1. / H) * a_inv * a_inv;
}

/**
 * @brief Computes \f$ dt / a^{3(\gamma - 1) + 1} \f$ for the current cosmology.
 *
 * @param a The scale-factor of interest.
 * @param param The current #cosmology.
 */
double hydro_kick_integrand(double a, void *param) {

  const struct cosmology *c = (const struct cosmology *)param;
  const double Omega_g = c->Omega_g;
  const double Omega_nu = cosmology_get_neutrino_density_param(c, a);
  const double Omega_r = (Omega_g > 0) ? Omega_g + Omega_nu : c->Omega_r;
  const double Omega_m = c->Omega_m;
  const double Omega_k = c->Omega_k;
  const double Omega_l = c->Omega_lambda;
  const double w_0 = c->w_0;
  const double w_a = c->w_a;
  const double H0 = c->H0;

  const double a_inv = 1. / a;
  const double E_z = E(Omega_r, Omega_m, Omega_k, Omega_l, w_0, w_a, a);
  const double H = H0 * E_z;

  /* Note: we can't use the pre-defined pow_gamma_xxx() function as
     as we need double precision accuracy for the GSL routine. */
  return (1. / H) * pow(a_inv, 3. * hydro_gamma_minus_one) * a_inv;
}

/**
 * @brief Computes \f$a dt\f$ for the current cosmology.
 *
 * @param a The scale-factor of interest.
 * @param param The current #cosmology.
 */
double hydro_kick_corr_integrand(double a, void *param) {

  const struct cosmology *c = (const struct cosmology *)param;
  const double Omega_g = c->Omega_g;
  const double Omega_nu = cosmology_get_neutrino_density_param(c, a);
  const double Omega_r = (Omega_g > 0) ? Omega_g + Omega_nu : c->Omega_r;
  const double Omega_m = c->Omega_m;
  const double Omega_k = c->Omega_k;
  const double Omega_l = c->Omega_lambda;
  const double w_0 = c->w_0;
  const double w_a = c->w_a;
  const double H0 = c->H0;

  const double E_z = E(Omega_r, Omega_m, Omega_k, Omega_l, w_0, w_a, a);
  const double H = H0 * E_z;

  return 1. / H;
}

/**
 * @brief Computes \f$ dt \f$ for the current cosmology.
 *
 * @param a The scale-factor of interest.
 * @param param The current #cosmology.
 */
double time_integrand(double a, void *param) {

  const struct cosmology *c = (const struct cosmology *)param;
  const double Omega_g = c->Omega_g;
  const double Omega_nu = cosmology_get_neutrino_density_param(c, a);
  const double Omega_r = (Omega_g > 0) ? Omega_g + Omega_nu : c->Omega_r;
  const double Omega_m = c->Omega_m;
  const double Omega_k = c->Omega_k;
  const double Omega_l = c->Omega_lambda;
  const double w_0 = c->w_0;
  const double w_a = c->w_a;
  const double H0 = c->H0;

  const double a_inv = 1. / a;
  const double E_z = E(Omega_r, Omega_m, Omega_k, Omega_l, w_0, w_a, a);
  const double H = H0 * E_z;

  return (1. / H) * a_inv;
}

/**
 * @brief Evaluates the neutrino density momentum integrand
 * \f$ x^2 \sqrt{x^2 + y^2} / (1+e^x) \f$, where
 * \f$ y = a*M_nu/(kb*T) \f$ is the neutrino mass scaled with redshift.
 *
 * This is used to evaluate the integral on (0, 1).
 *
 * @param x The momentum integration variable
 * @param param Neutrino mass y scaled by temperature at redshift of interest.
 */
double neutrino_density_integrand(double x, void *param) {
  double y = *(double *)param;
  return x * x * hypot(x, y) / (1 + exp(x));
}

/**
 * @brief Evaluates the transformed neutrino density momentum integrand
 * \f$ w^{-4} \sqrt{w^{-2} + y^2} / (1+e^{-w}) \f$, where
 * \f$ y = a*M_nu/(kb*T) \f$ is the neutrino mass scaled with redshift.
 *
 * This is used to evaluate the integral on (1, infinity).
 *
 * @param w The transformed momentum integration variable w=1/x
 * @param param Neutrino mass y scaled by temperature at redshift of interest.
 */
double neutrino_density_integrand_transformed(double w, void *param) {
  return neutrino_density_integrand(1. / w, param) / (w * w);
}

/**
 * @brief Performs the neutrino density momentum integrand
 * \f$ \int_0^\infty x^2 \sqrt{x^2 + y^2} / (1+e^x) dx \f$, where
 * \f$ y = a*M_nu/(kb*T) \f$ is the neutrino mass scaled with redshift,
 * without pre-factors.
 *
 * @param space The GSL working space
 * @param y Neutrino mass y scaled by temperature at redshift of interest.
 */
double neutrino_density_integrate(gsl_integration_workspace *space, double y) {
  double intermediate, abserr;

  double result = 0;

  gsl_function F1 = {&neutrino_density_integrand, &y};
  gsl_function F2 = {&neutrino_density_integrand_transformed, &y};
  /* Integrate between 0 and 1 */
  gsl_integration_qag(&F1, 0.0, 1.0, 0, 1.0e-12, GSL_workspace_size,
                      GSL_INTEG_GAUSS61, space, &intermediate, &abserr);
  result += intermediate;
  /* Integrate between 1 and infinity */
  gsl_integration_qag(&F2, 0.0, 1.0, 0, 1.0e-12, GSL_workspace_size,
                      GSL_INTEG_GAUSS61, space, &intermediate, &abserr);
  result += intermediate;

  return result;
}

/**
 * @brief Initialise the interpolation tables for the neutrino integrals.
 */
void cosmology_init_neutrino_tables(struct cosmology *c,
                                    const struct phys_const *phys_const) {

#ifdef HAVE_LIBGSL

  /* Initalise the GSL workspace */
  gsl_integration_workspace *space =
      gsl_integration_workspace_alloc(GSL_workspace_size);

  /* Create Omega_nu(a) interpolation table for massive neutrinos */
  if (swift_memalign("cosmo.table", (void **)&c->neutrino_density_interp_table,
                     SWIFT_STRUCT_ALIGNMENT,
                     cosmology_long_table_length * sizeof(double)) != 0)
    error("Failed to allocate long cosmology interpolation table");

  const double kb = phys_const->const_boltzmann_k;
  const double eV = phys_const->const_electron_volt;
  const size_t N_nu = c->N_nu;
  const double *M_nu = c->M_nu;
  const double pre_factor = 15. * pow(c->T_nu / c->T_CMB, 4) / pow(M_PI, 4);
  const double fermi_factor = 7. / 8. * pow(4. / 11., 4. / 3.);

  /* Iteratively find a time when the neutrinos are still relativistic */
  double a_start = c->a_begin;
  double N_eff, N_eff_prev = 0, err = 1;
  int iters = 0, max_iter = 1000;
  while (err > 1e-7 && iters < max_iter) {
    N_eff = 0;
    for (size_t j = 0; j < N_nu; j++) {
      /* Massless neutrino case */
      if (M_nu[j] == 0) {
        N_eff += c->N_eff / c->N_nu;
      } else {
        /* Integrate the FD distribtuion */
        double y = a_start * M_nu[j] * eV / (kb * c->T_nu);
        double result = neutrino_density_integrate(space, y);
        N_eff += result * pre_factor / fermi_factor;
      }
    }
    err = fabs(N_eff - N_eff_prev) / N_eff;
    a_start /= 1.1;
    N_eff_prev = N_eff;
    iters++;
  }

  double abs_err = fabs(N_eff - c->N_eff) / c->N_eff;

  if (iters == max_iter) {
    error("Could not find time when neutrinos were relativistic (max iter).");
  } else if (abs_err > 1e-5) {
    error("N_eff and neutrino density do not agree (err=%.10e)", abs_err);
  }

  c->log_a_nutab_begin = log(a_start);
  c->log_a_nutab_end = 0;  // scale-factor a=1

  /* Prepare a longer table of scale factors for the integral bounds */
  const double delta_fine_a =
      (c->log_a_nutab_end - c->log_a_nutab_begin) / cosmology_long_table_length;
  double *long_a_table = (double *)swift_malloc(
      "cosmo.table", cosmology_long_table_length * sizeof(double));
  for (int i = 0; i < cosmology_long_table_length; i++)
    long_a_table[i] = exp(c->log_a_nutab_begin + delta_fine_a * (i + 1));

  for (int i = 0; i < cosmology_long_table_length; i++) {
    double Onu_d_g = 0;  // Omega_nu / Omega_g
    for (size_t j = 0; j < N_nu; j++) {
      /* Massless neutrino case */
      if (M_nu[j] == 0) {
        Onu_d_g += (c->N_eff / c->N_nu) * fermi_factor;
      } else {
        /* Integrate the FD distribtuion */
        double y = long_a_table[i] * M_nu[j] * eV / (kb * c->T_nu);
        double result = neutrino_density_integrate(space, y);
        Onu_d_g += result * pre_factor;
      }
    }
    c->neutrino_density_interp_table[i] = Onu_d_g * c->Omega_g;
  }

  /* Free the temp array */
  swift_free("cosmo.table", long_a_table);

  /* Free the workspace and temp array */
  gsl_integration_workspace_free(space);

#else

  error("Code not compiled with GSL. Can't compute cosmology integrals.");

#endif
}

/**
 * @brief Initialise the interpolation tables for the integrals.
 */
void cosmology_init_tables(struct cosmology *c) {

#ifdef HAVE_LIBGSL

  /* Retrieve some constants */
  const double a_begin = c->a_begin;

  /* Allocate memory for the interpolation tables */
  if (swift_memalign("cosmo.table", (void **)&c->drift_fac_interp_table,
                     SWIFT_STRUCT_ALIGNMENT,
                     cosmology_table_length * sizeof(double)) != 0)
    error("Failed to allocate cosmology interpolation table");
  if (swift_memalign("cosmo.table", (void **)&c->grav_kick_fac_interp_table,
                     SWIFT_STRUCT_ALIGNMENT,
                     cosmology_table_length * sizeof(double)) != 0)
    error("Failed to allocate cosmology interpolation table");
  if (swift_memalign("cosmo.table", (void **)&c->hydro_kick_fac_interp_table,
                     SWIFT_STRUCT_ALIGNMENT,
                     cosmology_table_length * sizeof(double)) != 0)
    error("Failed to allocate cosmology interpolation table");
  if (swift_memalign("cosmo.table", (void **)&c->hydro_kick_corr_interp_table,
                     SWIFT_STRUCT_ALIGNMENT,
                     cosmology_table_length * sizeof(double)) != 0)
    error("Failed to allocate cosmology interpolation table");
  if (swift_memalign("cosmo.table", (void **)&c->time_interp_table,
                     SWIFT_STRUCT_ALIGNMENT,
                     cosmology_table_length * sizeof(double)) != 0)
    error("Failed to allocate cosmology interpolation table");
  if (swift_memalign("cosmo.table", (void **)&c->scale_factor_interp_table,
                     SWIFT_STRUCT_ALIGNMENT,
                     cosmology_table_length * sizeof(double)) != 0)
    error("Failed to allocate cosmology interpolation table");

  /* Prepare a table of scale factors for the integral bounds */
  const double delta_a =
      (c->log_a_end - c->log_a_begin) / cosmology_table_length;
  double *a_table = (double *)swift_malloc(
      "cosmo.table", cosmology_table_length * sizeof(double));
  for (int i = 0; i < cosmology_table_length; i++)
    a_table[i] = exp(c->log_a_begin + delta_a * (i + 1));

  /* Initalise the GSL workspace */
  gsl_integration_workspace *space =
      gsl_integration_workspace_alloc(GSL_workspace_size);

  double result, abserr;

  /* Integrate the drift factor \int_{a_begin}^{a_table[i]} dt/a^2 */
  gsl_function F = {&drift_integrand, c};
  for (int i = 0; i < cosmology_table_length; i++) {
    gsl_integration_qag(&F, a_begin, a_table[i], 0, 1.0e-10, GSL_workspace_size,
                        GSL_INTEG_GAUSS61, space, &result, &abserr);

    /* Store result */
    c->drift_fac_interp_table[i] = result;
  }

  /* Integrate the kick factor \int_{a_begin}^{a_table[i]} dt/a */
  F.function = &gravity_kick_integrand;
  for (int i = 0; i < cosmology_table_length; i++) {
    gsl_integration_qag(&F, a_begin, a_table[i], 0, 1.0e-10, GSL_workspace_size,
                        GSL_INTEG_GAUSS61, space, &result, &abserr);

    /* Store result */
    c->grav_kick_fac_interp_table[i] = result;
  }

  /* Integrate the kick factor \int_{a_begin}^{a_table[i]} dt/a^(3(g-1)+1) */
  F.function = &hydro_kick_integrand;
  for (int i = 0; i < cosmology_table_length; i++) {
    gsl_integration_qag(&F, a_begin, a_table[i], 0, 1.0e-10, GSL_workspace_size,
                        GSL_INTEG_GAUSS61, space, &result, &abserr);

    /* Store result */
    c->hydro_kick_fac_interp_table[i] = result;
  }

  /* Integrate the kick correction factor \int_{a_begin}^{a_table[i]} a dt */
  F.function = &hydro_kick_corr_integrand;
  for (int i = 0; i < cosmology_table_length; i++) {
    gsl_integration_qag(&F, a_begin, a_table[i], 0, 1.0e-10, GSL_workspace_size,
                        GSL_INTEG_GAUSS61, space, &result, &abserr);

    /* Store result */
    c->hydro_kick_corr_interp_table[i] = result;
  }

  /* Integrate the time \int_{a_begin}^{a_table[i]} dt */
  F.function = &time_integrand;
  for (int i = 0; i < cosmology_table_length; i++) {
    gsl_integration_qag(&F, a_begin, a_table[i], 0, 1.0e-10, GSL_workspace_size,
                        GSL_INTEG_GAUSS61, space, &result, &abserr);

    /* Store result */
    c->time_interp_table[i] = result;
  }

  /* Integrate the time \int_{0}^{a_begin} dt */
  gsl_integration_qag(&F, 0., a_begin, 0, 1.0e-10, GSL_workspace_size,
                      GSL_INTEG_GAUSS61, space, &result, &abserr);
  c->time_interp_table_offset = result;

  /* Integrate the time \int_{0}^{1} dt */
  gsl_integration_qag(&F, 0., 1, 0, 1.0e-13, GSL_workspace_size,
                      GSL_INTEG_GAUSS61, space, &result, &abserr);
  c->universe_age_at_present_day = result;

  /* Update the times */
  c->time_begin = cosmology_get_time_since_big_bang(c, c->a_begin);
  c->time_end = cosmology_get_time_since_big_bang(c, c->a_end);

  /*
   * Inverse t(a)
   */

  const double delta_t = (c->time_end - c->time_begin) / cosmology_table_length;

  /* index in the time_interp_table */
  int i_a = 0;

  for (int i_time = 0; i_time < cosmology_table_length; i_time++) {
    /* Current time
     * time_interp_table = \int_a_begin^a => no need of time_begin */
    double time_interp = delta_t * (i_time + 1);

    /* Find next time in time_interp_table */
    while (i_a < cosmology_table_length &&
           c->time_interp_table[i_a] <= time_interp) {
      i_a++;
    }

    /* Find linear interpolation scaling */
    double scale = 0;
    if (i_a != cosmology_table_length) {
      scale = time_interp - c->time_interp_table[i_a - 1];
      scale /= c->time_interp_table[i_a] - c->time_interp_table[i_a - 1];
    }

    scale += i_a;

    /* Compute interpolated scale factor */
    double log_a = c->log_a_begin + scale * (c->log_a_end - c->log_a_begin) /
                                        cosmology_table_length;
    c->scale_factor_interp_table[i_time] = exp(log_a) - c->a_begin;
  }

  /* Free the workspace and temp array */
  gsl_integration_workspace_free(space);
  swift_free("cosmo.table", a_table);

#else

  error("Code not compiled with GSL. Can't compute cosmology integrals.");

#endif
}

/**
 * @brief Initialises the #cosmology from the values read in the parameter file.
 *
 * @param params The parsed values.
 * @param us The current internal system of units.
 * @param phys_const The physical constants in the current system of units.
 * @param c The #cosmology to initialise.
 */
void cosmology_init(struct swift_params *params, const struct unit_system *us,
                    const struct phys_const *phys_const, struct cosmology *c) {

  /* Read in the cosmological parameters */
  c->Omega_m = parser_get_param_double(params, "Cosmology:Omega_m");
  c->Omega_r = parser_get_opt_param_double(params, "Cosmology:Omega_r", 0.);
  c->Omega_lambda = parser_get_param_double(params, "Cosmology:Omega_lambda");
  c->Omega_b = parser_get_param_double(params, "Cosmology:Omega_b");
  c->w_0 = parser_get_opt_param_double(params, "Cosmology:w_0", -1.);
  c->w_a = parser_get_opt_param_double(params, "Cosmology:w_a", 0.);
  c->h = parser_get_param_double(params, "Cosmology:h");

  /* Read in neutrino related quantities */
  c->Omega_g = parser_get_opt_param_double(params, "Cosmology:Omega_g", 0.);
  c->T_CMB = parser_get_opt_param_double(params, "Cosmology:T_CMB", 0);
  c->N_eff = parser_get_opt_param_double(params, "Cosmology:N_eff", 0);
  c->N_nu = parser_get_opt_param_int(params, "Cosmology:N_nu", 0);
  c->T_nu = parser_get_opt_param_double(params, "Cosmology:T_nu", 0);

  /* If there are neutrinos, load the optional mass array */
  if (c->N_nu > 0) {
    c->M_nu = calloc(c->N_nu, sizeof(double *));
    parser_get_opt_param_double_array(params, "Cosmology:M_nu", c->N_nu,
                                      c->M_nu);
  }

  /* Read the start and end of the simulation */
  c->a_begin = parser_get_param_double(params, "Cosmology:a_begin");
  c->a_end = parser_get_param_double(params, "Cosmology:a_end");
  c->log_a_begin = log(c->a_begin);
  c->log_a_end = log(c->a_end);
  c->time_base = (c->log_a_end - c->log_a_begin) / max_nr_timesteps;
  c->time_base_inv = 1. / c->time_base;

  /* If a_begin == a_end we hang */

  if (c->a_begin >= c->a_end)
    error("a_begin must be strictly before (and not equal to) a_end");

  /* Construct derived quantities */

  /* Dark-energy equation of state */
  c->w = cosmology_dark_energy_EoS(c->a_begin, c->w_0, c->w_a);

  /* Hubble constant in internal units */
  const double km = 1.e5 / units_cgs_conversion_factor(us, UNIT_CONV_LENGTH);
  const double H0_cgs =
      100. * c->h * (km / (1.e6 * phys_const->const_parsec)); /* s^-1 */
  c->H0 = H0_cgs * units_cgs_conversion_factor(us, UNIT_CONV_TIME);
  c->Hubble_time = 1. / c->H0;

  /* Critical density at present day */
  c->critical_density_0 =
      3. * c->H0 * c->H0 / (8. * M_PI * phys_const->const_newton_G);

  /* Initialise quantities related to neutrinos and relativisitc species */
  cosmology_neutrino_init(params, us, phys_const, c);

  /* If necessary, initialize the neutrino density interpolation table */
  if (c->M_nu_tot > 0) {
    c->neutrino_density_interp_table = NULL;
    cosmology_init_neutrino_tables(c, phys_const);
    c->Omega_nu = cosmology_get_neutrino_density_param(c, 1);  // scale-factor 1
  } else {
    /* All massless case */
    const double fermi_factor = 7. / 8. * pow(4. / 11., 4. / 3.);
    c->Omega_nu = c->Omega_g * c->N_eff * fermi_factor;
  }

  /* Curvature density (for closure) */
  if (c->Omega_g > 0) {
    c->Omega_k = 1. - (c->Omega_m + c->Omega_nu + c->Omega_g + c->Omega_lambda);
  } else {
    c->Omega_k = 1. - (c->Omega_m + c->Omega_r + c->Omega_lambda);
  }

  /* Initialise the interpolation tables */
  c->drift_fac_interp_table = NULL;
  c->grav_kick_fac_interp_table = NULL;
  c->hydro_kick_fac_interp_table = NULL;
  c->time_interp_table = NULL;
  c->time_interp_table_offset = 0.;
  cosmology_init_tables(c);

  /* Set remaining variables to alid values */
  cosmology_update(c, phys_const, 0);

  /* Update the times */
  c->time_begin = cosmology_get_time_since_big_bang(c, c->a_begin);
  c->time_end = cosmology_get_time_since_big_bang(c, c->a_end);

  /* Initialise the old values to a valid state */
  c->a_old = c->a_begin;
  c->z_old = 1. / c->a_old - 1.;
}

/**
 * @brief Initialises and verifies quantities related to neutrinos and
 * relativistic species
 *
 * @param params The parsed values.
 * @param us The current internal system of units.
 * @param phys_const The physical constants in the current system of units.
 * @param c The #cosmology to initialise.
 */
void cosmology_neutrino_init(struct swift_params *params,
                             const struct unit_system *us,
                             const struct phys_const *phys_const,
                             struct cosmology *c) {

  /* Find the total neutrino mass (eV) */
  c->M_nu_tot = 0;
  for (size_t i = 0; i < c->N_nu; i++) {
    c->M_nu_tot += c->M_nu[i];
  }

  /* The CMB temperature is fixed by Omega_g, but users may specify either */
  const double hbar = phys_const->const_planck_h / (2 * M_PI);
  const double cvel = phys_const->const_speed_light_c;
  const double kb = phys_const->const_boltzmann_k;
  const double crit_energy_density = c->critical_density_0 * cvel * cvel;
  const double e_pre_factor = 15.0 / (M_PI * M_PI) * pow(cvel * hbar, 3);

  /* Ensure that Omega_g and T_CMB are not both specified */
  if (c->T_CMB > 0 && c->Omega_g > 0) {
    error("T_CMB and Omega_g should not both be specified.");
  }

  /* Infer Omega_g from T_CMB or vice versa */
  if (c->T_CMB > 0) {
    c->Omega_g = pow(kb * c->T_CMB, 4) / crit_energy_density / e_pre_factor;
  } else if (c->Omega_g > 0) {
    c->T_CMB = pow(c->Omega_g * crit_energy_density * e_pre_factor, 0.25) / kb;
  }

  /* Ensure that neutrino numbers and Omega_g are consistent */
  if (c->N_nu > 0 && c->Omega_g == 0) {
    error("Specify Omega_g or T_CMB to include neutrinos.");
  } else if (c->N_eff > 0 && c->Omega_g == 0) {
    error("Specify Omega_g or T_CMB to include relativistic species.");
  }

  /* Ensure that T_nu and N_eff are consistent */
  if (c->T_nu > 0 && c->N_eff > 0) {
    error("T_nu and N_eff should not both be specified.");
  } else if (c->T_nu == 0) {
    c->T_nu = c->T_CMB * pow(c->N_eff / c->N_nu, 0.25) * pow(4. / 11., 1. / 3.);
  } else {
    c->N_eff = c->N_nu * pow(c->T_nu / c->T_CMB, 4) * pow(11. / 4., 4. / 3.);
  }

  /* Ensure that neutrino masses and T_nu are consistent */
  if (c->N_nu > 0 && c->T_nu == 0) {
    error("Specify T_nu or N_eff to include neutrinos.");
  }

  /* Ensure that Omega_g and Omega_r are not both specified */
  if (c->Omega_r > 0 && c->Omega_g > 0) {
    error("Omega_r and Omega_g (or T_CMB) should not both be specified.");
  }
}

/**
 * @brief Initialise the #cosmology for non-cosmological time-integration
 *
 * Essentially sets all constants to 1 or 0.
 *
 * @param c The #cosmology to initialise.
 */
void cosmology_init_no_cosmo(struct cosmology *c) {

  c->Omega_m = 0.;
  c->Omega_r = 0.;
  c->Omega_k = 0.;
  c->Omega_lambda = 0.;
  c->Omega_b = 0.;
  c->w_0 = 0.;
  c->w_a = 0.;
  c->h = 1.;
  c->w = -1.;

  c->Omega_g = 0.;
  c->Omega_nu = 0;
  c->T_CMB = 0;
  c->T_nu = 0;
  c->N_nu = 0;
  c->N_eff = 0;
  c->M_nu_tot = 0;

  c->a_begin = 1.;
  c->a_end = 1.;
  c->log_a_begin = 0.;
  c->log_a_end = 0.;

  c->H = 0.;
  c->H0 = 0.;
  c->a = 1.;
  c->z = 0.;
  c->a_inv = 1.;
  c->a2_inv = 1.;
  c->a3_inv = 1.;
  c->a_factor_internal_energy = 1.;
  c->a_factor_pressure = 1.;
  c->a_factor_sound_speed = 1.;
  c->a_factor_mu = 1.;
  c->a_factor_Balsara_eps = 1.;
  c->a_factor_hydro_accel = 1.;
  c->a_factor_grav_accel = 1.;

  c->a_old = 1.;
  c->z_old = 0.;

  c->critical_density = 0.;
  c->critical_density_0 = 0.;
  c->mean_density = 0.;
  c->mean_density_Omega_b = 0;

  c->time_step_factor = 1.;

  c->a_dot = 0.;
  c->time = 0.;
  c->universe_age_at_present_day = 0.;
  c->Hubble_time = 0.;
  c->lookback_time = 0.;

  /* Initialise the interpolation tables */
  c->drift_fac_interp_table = NULL;
  c->grav_kick_fac_interp_table = NULL;
  c->hydro_kick_fac_interp_table = NULL;
  c->hydro_kick_corr_interp_table = NULL;
  c->time_interp_table = NULL;
  c->time_interp_table_offset = 0.;
  c->log_a_nutab_begin = 0.;
  c->log_a_nutab_end = 0.;
  c->scale_factor_interp_table = NULL;

  c->time_begin = 0.;
  c->time_end = 0.;
}

/**
 * @brief Computes the cosmology factor that enters the drift operator.
 *
 * Computes \f$ \int_{a_start}^{a_end} dt/a^2 \f$ using the interpolation table.
 *
 * @param c The current #cosmology.
 * @param ti_start the (integer) time of the start of the drift.
 * @param ti_end the (integer) time of the end of the drift.
 */
double cosmology_get_drift_factor(const struct cosmology *c,
                                  const integertime_t ti_start,
                                  const integertime_t ti_end) {

#ifdef SWIFT_DEBUG_CHECKS
  if (ti_end < ti_start) error("ti_end must be >= ti_start");
#endif

  const double a_start = c->log_a_begin + ti_start * c->time_base;
  const double a_end = c->log_a_begin + ti_end * c->time_base;

  const double int_start = interp_table(c->drift_fac_interp_table, a_start,
                                        c->log_a_begin, c->log_a_end);
  const double int_end = interp_table(c->drift_fac_interp_table, a_end,
                                      c->log_a_begin, c->log_a_end);

  return int_end - int_start;
}

/**
 * @brief Computes the cosmology factor that enters the gravity kick operator.
 *
 * Computes \f$ \int_{a_start}^{a_end} dt/a \f$ using the interpolation table.
 *
 * @param c The current #cosmology.
 * @param ti_start the (integer) time of the start of the drift.
 * @param ti_end the (integer) time of the end of the drift.
 */
double cosmology_get_grav_kick_factor(const struct cosmology *c,
                                      const integertime_t ti_start,
                                      const integertime_t ti_end) {

#ifdef SWIFT_DEBUG_CHECKS
  if (ti_end < ti_start) error("ti_end must be >= ti_start");
#endif

  const double a_start = c->log_a_begin + ti_start * c->time_base;
  const double a_end = c->log_a_begin + ti_end * c->time_base;

  const double int_start = interp_table(c->grav_kick_fac_interp_table, a_start,
                                        c->log_a_begin, c->log_a_end);
  const double int_end = interp_table(c->grav_kick_fac_interp_table, a_end,
                                      c->log_a_begin, c->log_a_end);

  return int_end - int_start;
}

/**
 * @brief Computes the cosmology factor that enters the hydro kick operator.
 *
 * Computes \f$ \int_{a_start}^{a_end} dt/a^{3(gamma - 1)} \f$ using the
 * interpolation table.
 *
 * @param c The current #cosmology.
 * @param ti_start the (integer) time of the start of the drift.
 * @param ti_end the (integer) time of the end of the drift.
 */
double cosmology_get_hydro_kick_factor(const struct cosmology *c,
                                       const integertime_t ti_start,
                                       const integertime_t ti_end) {

#ifdef SWIFT_DEBUG_CHECKS
  if (ti_end < ti_start) error("ti_end must be >= ti_start");
#endif

  const double a_start = c->log_a_begin + ti_start * c->time_base;
  const double a_end = c->log_a_begin + ti_end * c->time_base;

  const double int_start = interp_table(c->hydro_kick_fac_interp_table, a_start,
                                        c->log_a_begin, c->log_a_end);
  const double int_end = interp_table(c->hydro_kick_fac_interp_table, a_end,
                                      c->log_a_begin, c->log_a_end);

  return int_end - int_start;
}

/**
 * @brief Computes the cosmology factor that enters the hydro kick correction
 * operator for the meshless schemes (GIZMO-MFV).
 *
 * Computes \f$ \int_{a_start}^{a_end} a dt \f$ using the interpolation table.
 *
 * @param c The current #cosmology.
 * @param ti_start the (integer) time of the start of the drift.
 * @param ti_end the (integer) time of the end of the drift.
 */
double cosmology_get_corr_kick_factor(const struct cosmology *c,
                                      const integertime_t ti_start,
                                      const integertime_t ti_end) {

#ifdef SWIFT_DEBUG_CHECKS
  if (ti_end < ti_start) error("ti_end must be >= ti_start");
#endif

  const double a_start = c->log_a_begin + ti_start * c->time_base;
  const double a_end = c->log_a_begin + ti_end * c->time_base;

  const double int_start = interp_table(c->hydro_kick_corr_interp_table,
                                        a_start, c->log_a_begin, c->log_a_end);
  const double int_end = interp_table(c->hydro_kick_corr_interp_table, a_end,
                                      c->log_a_begin, c->log_a_end);

  return int_end - int_start;
}

/**
 * @brief Computes the cosmology factor that enters the thermal variable kick
 * operator.
 *
 * Computes \f$ \int_{a_start}^{a_end} dt/a^2 \f$ using the interpolation table.
 *
 * @param c The current #cosmology.
 * @param ti_start the (integer) time of the start of the drift.
 * @param ti_end the (integer) time of the end of the drift.
 */
double cosmology_get_therm_kick_factor(const struct cosmology *c,
                                       const integertime_t ti_start,
                                       const integertime_t ti_end) {

#ifdef SWIFT_DEBUG_CHECKS
  if (ti_end < ti_start) error("ti_end must be >= ti_start");
#endif

  const double a_start = c->log_a_begin + ti_start * c->time_base;
  const double a_end = c->log_a_begin + ti_end * c->time_base;

  const double int_start = interp_table(c->drift_fac_interp_table, a_start,
                                        c->log_a_begin, c->log_a_end);
  const double int_end = interp_table(c->drift_fac_interp_table, a_end,
                                      c->log_a_begin, c->log_a_end);

  return int_end - int_start;
}

/**
 * @brief Compute the cosmic time (in internal units) between two points
 * on the integer time line.
 *
 * @param c The current #cosmology.
 * @param ti_start the (integer) time of the start.
 * @param ti_end the (integer) time of the end.
 */
double cosmology_get_delta_time(const struct cosmology *c,
                                const integertime_t ti_start,
                                const integertime_t ti_end) {

#ifdef SWIFT_DEBUG_CHECKS
  if (ti_end < ti_start) error("ti_end must be >= ti_start");
#endif

  const double log_a_start = c->log_a_begin + ti_start * c->time_base;
  const double log_a_end = c->log_a_begin + ti_end * c->time_base;

  /* Time between a_begin and a_start */
  const double t1 = interp_table(c->time_interp_table, log_a_start,
                                 c->log_a_begin, c->log_a_end);

  /* Time between a_begin and a_end */
  const double t2 = interp_table(c->time_interp_table, log_a_end,
                                 c->log_a_begin, c->log_a_end);

  return t2 - t1;
}

/**
 * @brief Compute the cosmic time (in internal units) between two scale factors
 *
 * @param c The current #cosmology.
 * @param a_start the starting scale factor
 * @param a_end the ending scale factor
 */
double cosmology_get_delta_time_from_scale_factors(const struct cosmology *c,
                                                   const double a_start,
                                                   const double a_end) {

#ifdef SWIFT_DEBUG_CHECKS
  if (a_end < a_start) error("a_end must be >= a_start");
#endif

  const double log_a_start = log(a_start);
  const double log_a_end = log(a_end);

  /* Time between a_begin and a_start */
  const double t1 = interp_table(c->time_interp_table, log_a_start,
                                 c->log_a_begin, c->log_a_end);

  /* Time between a_begin and a_end */
  const double t2 = interp_table(c->time_interp_table, log_a_end,
                                 c->log_a_begin, c->log_a_end);

  return t2 - t1;
}

/**
 * @brief Compute scale factor from time since big bang (in internal units).
 *
 * @param c The current #cosmology.
 * @param t time since the big bang
 * @return The scale factor.
 */
double cosmology_get_scale_factor(const struct cosmology *c, double t) {
  /* scale factor between time_begin and t */
  const double a =
      interp_table(c->scale_factor_interp_table, t, c->time_interp_table_offset,
                   c->universe_age_at_present_day);
  return a + c->a_begin;
}

/**
 * @brief Compute neutrino density parameter Omega_nu at the given scale-factor
 * This is the effective present day value, i.e. must be multiplied by (1+z)^4
 *
 * @param c The current #cosmology.
 * @param a The scale factor
 * @return The density parameter
 */
double cosmology_get_neutrino_density_param(const struct cosmology *c,
                                            double a) {

  /* All massless or no neutrinos at all */
  if (c->M_nu_tot == 0) {
    return c->Omega_nu;
  } else {
    /* Accounting for the R/NR transition of massive neutrinos */
    return interp_long_table(c->neutrino_density_interp_table, log(a),
                             c->log_a_nutab_begin, c->log_a_nutab_end);
  }
}

/**
 * @brief Prints the #cosmology model to stdout.
 */
void cosmology_print(const struct cosmology *c) {

  message(
      "Density parameters: [O_m, O_l, O_b, O_nu, O_k, O_r, O_g] = [%f, %f, %f, "
      "%f, "
      "%f, %f %f]",
      c->Omega_m, c->Omega_lambda, c->Omega_b, c->Omega_nu, c->Omega_k,
      c->Omega_r, c->Omega_g);
  message("Dark energy equation of state: w_0=%f w_a=%f", c->w_0, c->w_a);
  message("Hubble constant: h = %f, H_0 = %e U_t^(-1)", c->h, c->H0);
  message("Hubble time: 1/H0 = %e U_t", c->Hubble_time);
  message("Universe age at present day: %e U_t",
          c->universe_age_at_present_day);
}

void cosmology_clean(struct cosmology *c) {

  swift_free("cosmo.table", c->drift_fac_interp_table);
  swift_free("cosmo.table", c->grav_kick_fac_interp_table);
  swift_free("cosmo.table", c->hydro_kick_fac_interp_table);
  swift_free("cosmo.table", c->hydro_kick_corr_interp_table);
  swift_free("cosmo.table", c->time_interp_table);
  swift_free("cosmo.table", c->scale_factor_interp_table);
  if (c->M_nu_tot > 0) {
    swift_free("cosmo.table", c->neutrino_density_interp_table);
  }
}

#ifdef HAVE_HDF5
void cosmology_write_model(hid_t h_grp, const struct cosmology *c) {

  io_write_attribute_d(h_grp, "a_beg", c->a_begin);
  io_write_attribute_d(h_grp, "a_end", c->a_end);
  io_write_attribute_d(h_grp, "time_beg [internal units]", c->time_begin);
  io_write_attribute_d(h_grp, "time_end [internal units]", c->time_end);
  io_write_attribute_d(h_grp, "Universe age [internal units]", c->time);
  io_write_attribute_d(h_grp, "Lookback time [internal units]",
                       c->lookback_time);
  io_write_attribute_d(h_grp, "h", c->h);
  io_write_attribute_d(h_grp, "H0 [internal units]", c->H0);
  io_write_attribute_d(h_grp, "H [internal units]", c->H);
  io_write_attribute_d(h_grp, "Hubble time [internal units]", c->Hubble_time);
  io_write_attribute_d(h_grp, "Omega_m", c->Omega_m);
  io_write_attribute_d(h_grp, "Omega_r", c->Omega_r);
  io_write_attribute_d(h_grp, "Omega_b", c->Omega_b);
  io_write_attribute_d(h_grp, "Omega_k", c->Omega_k);
  io_write_attribute_d(h_grp, "Omega_lambda", c->Omega_lambda);

  io_write_attribute_d(h_grp, "Omega_g", c->Omega_g);
  io_write_attribute_d(h_grp, "Omega_nu", c->Omega_nu);
  io_write_attribute_d(h_grp, "T_CMB", c->T_CMB);
  io_write_attribute_d(h_grp, "T_nu", c->T_nu);
  io_write_attribute_d(h_grp, "N_nu", c->N_nu);
  io_write_attribute_d(h_grp, "N_eff", c->N_eff);
  io_write_attribute_d(h_grp, "M_nu_tot", c->M_nu_tot);
  io_write_attribute(h_grp, "M_nu", DOUBLE, c->M_nu, c->N_nu);

  io_write_attribute_d(h_grp, "w_0", c->w_0);
  io_write_attribute_d(h_grp, "w_a", c->w_a);
  io_write_attribute_d(h_grp, "w", c->w);
  io_write_attribute_d(h_grp, "Redshift", c->z);
  io_write_attribute_d(h_grp, "Scale-factor", c->a);
  io_write_attribute_d(h_grp, "Critical density [internal units]",
                       c->critical_density);
}
#endif

/**
 * @brief Write a cosmology struct to the given FILE as a stream of bytes.
 *
 * @param cosmology the struct
 * @param stream the file stream
 */
void cosmology_struct_dump(const struct cosmology *cosmology, FILE *stream) {
  restart_write_blocks((void *)cosmology, sizeof(struct cosmology), 1, stream,
                       "cosmology", "cosmology function");
}

/**
 * @brief Restore a cosmology struct from the given FILE as a stream of
 * bytes.
 *
 * @param enabled whether cosmology is enabled.
 * @param cosmology the struct
 * @param stream the file stream
 */
void cosmology_struct_restore(int enabled, struct cosmology *cosmology,
                              const struct phys_const *phys_const,
                              FILE *stream) {
  restart_read_blocks((void *)cosmology, sizeof(struct cosmology), 1, stream,
                      NULL, "cosmology function");

  /* Re-initialise the tables if using a cosmology. */
  if (enabled) {
    cosmology_init_neutrino_tables(cosmology, phys_const);
    cosmology_init_tables(cosmology);
  }
}
