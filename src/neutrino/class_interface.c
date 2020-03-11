/*******************************************************************************
 * This file is part of SWIFT.
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

/* Config parameters. */
#include "../../config.h"

/* This object's header. */
#include "class_interface.h"

/* Run CLASS and compute the neutrino perturbation */
void rend_perturb_from_class(struct renderer *rend, struct swift_params *params,
                             const struct engine *e) {
#ifdef WITH_CLASS_INTERFACE

  /* Load internal units and physical constants */
  const struct unit_system *us = e->internal_units;
  const struct phys_const *pc = e->physical_constants;
  const struct cosmology *cosmo = e->cosmology;

  /* CLASS to internal units conversion factor */
  const double Mpc_to_cm = _Mpc_over_m_ * 100;  // CLASS uses hard-coded Mpc's
  const double unit_length_factor = Mpc_to_cm / us->UnitLength_in_cgs;
  const double unit_time_factor = unit_length_factor / pc->const_speed_light_c;

  message("Converting CLASS Units:");
  message("(CLASS) Unit system: U_L = \t %.6e cm", Mpc_to_cm);
  message("(CLASS) Unit system: U_T = \t %.6e s",
          unit_time_factor * us->UnitTime_in_cgs);
  message("to:");
  message("(internal) Unit system: U_L = \t %.6e cm", us->UnitLength_in_cgs);
  message("(internal) Unit system: U_T = \t %.6e s", us->UnitTime_in_cgs);

  /* Define the CLASS structures */
  struct precision pr;  /* for precision parameters */
  struct background ba; /* for cosmological background */
  struct thermo th;     /* for thermodynamics */
  struct perturbs pt;   /* for source functions */
  struct transfers tr;  /* for transfer functions */
  struct primordial pm; /* for primordial spectra */
  struct spectra sp;    /* for output spectra */
  struct nonlinear nl;  /* for non-linear spectra */
  struct lensing le;    /* for lensed spectra */
  struct output op;     /* for output files */
  ErrorMsg errmsg;      /* for CLASS-specific error messages */

  /* If no class .ini file was specified, infer parameters from the cosmology */
  if (rend->class_ini_fname[0] == '\0') {
    message("Inferring CLASS parameters from the cosmology.");

    /* Infer CLASS parameters from the cosmology module */
    struct file_content fc;
    rend_infer_class_parameters(rend, e, &fc);

    if (input_init(&fc, &pr, &ba, &th, &pt, &tr, &pm, &sp, &nl, &le, &op,
                   errmsg) == _FAILURE_) {
      error("Error running input_init_from_arguments \n=>%s\n", errmsg);
    }
  } else {
    /* Otherwise, initialize CLASS with the parameter files */
    int class_argc = 2;
    char *class_argv[] = {"", rend->class_ini_fname, rend->class_pre_fname};

    message("Reading CLASS parameters from '%s'.", rend->class_ini_fname);
    message("Reading CLASS precision parameters from '%s'.",
            rend->class_pre_fname);

    if (input_init_from_arguments(class_argc, class_argv, &pr, &ba, &th, &pt,
                                  &tr, &pm, &sp, &nl, &le, &op,
                                  errmsg) == _FAILURE_) {
      error("Error running input_init_from_arguments \n=>%s\n", errmsg);
    }
  }

  message("Running CLASS.");

  if (background_init(&pr, &ba) == _FAILURE_) {
    error("Error running background_init \n%s\n", ba.error_message);
  }

  /* Compare cosmological parameters between CLASS & SWIFT */
  rend_print_cosmology(&ba, cosmo);

  if (thermodynamics_init(&pr, &ba, &th) == _FAILURE_) {
    error("Error in thermodynamics_init \n%s\n", th.error_message);
  }

  if (perturb_init(&pr, &ba, &th, &pt) == _FAILURE_) {
    error("Error in perturb_init \n%s\n", pt.error_message);
  }

  /* Try getting a source */
  int index_md = pt.index_md_scalars;  // scalar mode
  int index_ic = 0;                    // index of the initial condition
  int index_tp;                        // type of source function

  /* Size of the perturbations */
  size_t k_size = pt.k_size[index_md];
  size_t tau_size = pt.tau_size;

  /* The number of transfer functions to be read */
  const size_t n_functions = 2;

  /* Little h, which CLASS uses but Swift doesn't */
  const double h = ba.h;

  /* Vector of the wavenumbers */
  rend->transfer.k_size = k_size;
  rend->transfer.k = (double *)calloc(k_size, sizeof(double));

  /* Vector of the conformal times at which the perturbation is sampled */
  rend->transfer.tau_size = tau_size;
  rend->transfer.log_tau = (double *)calloc(tau_size, sizeof(double));

  /* The number of transfer functions to be read */
  rend->transfer.n_functions = n_functions;

  /* What functions should be read */
  int *functions = malloc(n_functions * sizeof(double));
  functions[0] = pt.index_tp_delta_ncdm1;
  functions[1] = pt.index_tp_delta_cdm;

  /* Vector with the transfer functions T(tau, k) */
  rend->transfer.delta =
      (double *)calloc(n_functions * k_size * tau_size, sizeof(double));

  /* Read out the log conformal times */
  for (size_t index_tau = 0; index_tau < tau_size; index_tau++) {
    /* Convert tau from Mpc to U_T */
    double tau = pt.tau_sampling[index_tau] * unit_time_factor;
    rend->transfer.log_tau[index_tau] = log(tau);
  }

  /* Read out the wavenumbers */
  for (size_t index_k = 0; index_k < k_size; index_k++) {
    /* Convert k from h/Mpc to 1/U_L */
    double k = pt.k[index_md][index_k] * h / unit_length_factor;
    rend->transfer.k[index_k] = k;
  }

  /* Convert and store the transfer functions */
  for (size_t index_tau = 0; index_tau < tau_size; index_tau++) {
    for (size_t index_k = 0; index_k < k_size; index_k++) {
      for (size_t index_func = 0; index_func < n_functions; index_func++) {
        index_tp = functions[index_func];  // type of source function
        double p = pt.sources[index_md][index_ic * pt.tp_size[index_md] +
                                        index_tp][index_tau * k_size + index_k];

        /* Convert transfer functions from CLASS format to CAMB/HeWon/icgen/
         *  Eisenstein-Hu format by multiplying by -1/k^2.
         */
        double k = pt.k[index_md][index_k] * h / unit_length_factor;
        double T = -p / k / k;
        rend->transfer.delta[tau_size * k_size * index_func +
                             k_size * index_tau + index_k] = T;
      }
    }
  }

  message("The neutrino density is sampled at %zu * %zu points.", k_size,
          tau_size);

  /* Pre-empt segfault in CLASS if there is no interacting dark radiation */
  if (ba.has_idr == _FALSE_) {
    pt.alpha_idm_dr = (double *)malloc(0);
    pt.beta_idr = (double *)malloc(0);
  }

  /* Close CLASS again */
  if (perturb_free(&pt) == _FAILURE_) {
    error("Error in freeing class memory \n%s\n", pt.error_message);
  }

  if (thermodynamics_free(&th) == _FAILURE_) {
    error("Error in thermodynamics_free \n%s\n", th.error_message);
  }

  if (background_free(&ba) == _FAILURE_) {
    error("Error in background_free \n%s\n", ba.error_message);
  }
#else
  error("No CLASS library found. Cannot compute transfer functions.");
#endif
}

/* Print cosmological information from SWIFT & CLASS */
void rend_print_cosmology(struct background *ba,
                          const struct cosmology *cosmo) {
#ifdef WITH_CLASS_INTERFACE
  /* For diagnostics, find the total ncdm pressure from CLASS at z=0 */
  double rho_ncdm_tot = 0;
  double p_ncdm_tot = 0;
  for (int index_ncdm = 0; index_ncdm < ba->N_ncdm; index_ncdm++) {
    double rho_ncdm;
    double p_ncdm;
    background_ncdm_momenta(
        ba->q_ncdm_bg[index_ncdm], ba->w_ncdm_bg[index_ncdm],
        ba->q_size_ncdm_bg[index_ncdm], ba->M_ncdm[index_ncdm],
        ba->factor_ncdm[index_ncdm], 0.0, NULL, &rho_ncdm, &p_ncdm, NULL, NULL);

    rho_ncdm_tot += rho_ncdm;
    p_ncdm_tot += p_ncdm;
  }

  /* The relativistic part of the neutrino density */
  double rho_ncdm_rel_tot = 3. * p_ncdm_tot;
  double Omega_nu_mat = (rho_ncdm_tot - rho_ncdm_rel_tot) / pow(ba->H0, 2);
  double Omega_nu_rel = rho_ncdm_rel_tot / pow(ba->H0, 2);

#ifdef NEUTRINO_BACKGROUND
  double Og = cosmo->Omega_g;
  double Or = cosmo->Omega_g;
#else
  double Or = cosmo->Omega_r;
  double Og = cosmo->Omega_r;
#endif

  /*
   * WITH NEUTRINO_BACKGROUND enabled:
   * In our implementation, neutrinos are treated separately from
   * matter and radiation. In CLASS, the neutrino density is split
   * into a radiation part ~3*p_nu and matter part rho_nu - 3*p_nu.
   * This difference is taken into account in the comparison below.
   *
   * WITHOUT NEUTRINO_BACKGROUND enabled:
   * Neutrinos are treated as pure matter.
   */

  message("Comparison of density parameters:");
  message(
      "(internal) \t[O_m, O_l, O_b, O_nu, O_k, O_r, O_g] = "
      "[%f, %f, %f, %f, %f, %f, %f]",
      cosmo->Omega_m, cosmo->Omega_lambda, cosmo->Omega_b, cosmo->Omega_nu,
      cosmo->Omega_k, Or, Og);
  message(
      "(CLASS) \t[O_m, O_l, O_b, O_nu, O_k, O_r, O_g] = "
      "[%f, %f, %f, %f, %f, %f, %f]",
      ba->Omega0_m - Omega_nu_mat, ba->Omega0_lambda, ba->Omega0_b,
      ba->Omega0_ncdm_tot, ba->Omega0_k,
      ba->Omega0_r - ba->Omega0_ur - Omega_nu_rel, ba->Omega0_g);
  message("Effective matter & radiation contributions of the neutrinos:");
  message("(CLASS) Massive neutrinos: \t[O_nu_m, O_nu_r] = [%f, %f]",
          Omega_nu_mat, Omega_nu_rel);
  message("(CLASS) Extra rel. species: \t[O_ur] = [%f]", ba->Omega0_ur);
  message("(CLASS) Extra fluid species: \t[O_fld] = [%f]", ba->Omega0_fld);

#else
  error("No CLASS library found. Cannot compute transfer functions.");
#endif
}

/* Infer CLASS parameters from the cosmology module */
void rend_infer_class_parameters(struct renderer *rend, const struct engine *e,
                                 struct file_content *fc) {
#ifdef WITH_CLASS_INTERFACE
  const struct cosmology *cosmo = e->cosmology;
  const struct unit_system *us = e->internal_units;

  /* CLASS to internal units conversion factor */
  const double Mpc_to_cm = _Mpc_over_m_ * 100;  // CLASS uses hard-coded Mpc's
  const double unit_length_factor = Mpc_to_cm / us->UnitLength_in_cgs;

  /* For CLASS-specific error messages */
  ErrorMsg errmsg;

  /* Extract the main cosmological parameters */
  double h = cosmo->h;
  double Omega_m = cosmo->Omega_m;
  double Omega_b = cosmo->Omega_b;
  double Omega_k = cosmo->Omega_k;
  double Omega_c = Omega_m - Omega_b;
  double Omega_l = cosmo->Omega_lambda;
#ifndef NEUTRINO_BACKGROUND
  double Omega_r = cosmo->Omega_r;
#endif

  /* Start and end of the simulation */
  double z_start = 1. / cosmo->a_begin - 1.;
  double z_end = 1. / cosmo->a_end - 1.;

  /* Compute redshifts of interest around which the precision is highest */
  double z_0 = z_start;
  double z_1 = z_start + 0.4 * (z_end - z_start);
  double z_2 = z_start + 0.8 * (z_end - z_start);
  double z_3 = z_start + 0.9 * (z_end - z_start);
  double z_4 = z_start + 0.95 * (z_end - z_start);
  double z_5 = z_end;

  /* The maximum and minimum wavenumbers necessary, given the box size */
  double box_len = e->s->dim[0];  // U_L
  double box_len_Mpc = box_len / unit_length_factor;
  double k_max = sqrt(3.) * M_PI * e->mesh->N / box_len_Mpc;
  double k_min = 2 * M_PI / box_len;  // not used currently

  /* To account for the unequal CLASS k steps, we multiply by two */
  double requested_k_max = k_max * 2.0;

  /* Default precision values */
  int fluid_approx = 3;         // method 3
  int quadrature_strategy = 3;  // strategy 3
  int nbins = 5;                // number of momentum bins
  int lmax_ncdm = 20;           // highest multipole for ncdm
  int k_per_decade = 50;

  message("Precision required at z = %.1f, %.1f, %.1f, %.1f, %.1f, %.1f", z_0,
          z_1, z_2, z_3, z_4, z_5);
  message("Box size means we need wavenumbers in range [%f,%f] 1/Mpc", k_min,
          k_max);

  /* Counter for the parameters */
  int num = 0;

  /* We will put the parameters into a new file: internal.ini */
  class_parser_init(fc, 50, "internal.ini", errmsg);

  /* Insert the parameters */
  strcpy(fc->name[num], "output");
  strcpy(fc->value[num++], "dTk,vTk");  // we want density & velocity transfers
  strcpy(fc->name[num], "z_pk");
  sprintf(fc->value[num++], "%f, %f, %f, %f, %f, %f", z_0, z_1, z_2, z_3, z_4,
          z_5);
  strcpy(fc->name[num], "h");
  sprintf(fc->value[num++], "%e", h);
  strcpy(fc->name[num], "P_k_max_1/Mpc");
  sprintf(fc->value[num++], "%e", requested_k_max);
  strcpy(fc->name[num], "k_per_decade_for_pk");
  sprintf(fc->value[num++], "%i", k_per_decade);
  strcpy(fc->name[num], "extra metric transfer functions");
  strcpy(fc->value[num++], "yes");
  strcpy(fc->name[num], "Nbody gauge transfer functions");
  strcpy(fc->value[num++], "yes");
  strcpy(fc->name[num], "Omega_Lambda");
  sprintf(fc->value[num++], "%e", Omega_l);
  strcpy(fc->name[num], "Omega_b");
  sprintf(fc->value[num++], "%e", Omega_b);
  strcpy(fc->name[num], "Omega_cdm");
  sprintf(fc->value[num++], "%e", Omega_c);
  strcpy(fc->name[num], "Omega_k");
  sprintf(fc->value[num++], "%e", Omega_k);
  strcpy(fc->name[num], "l_max_ncdm");
  sprintf(fc->value[num++], "%i", lmax_ncdm);
  strcpy(fc->name[num], "ncdm_fluid_approximation");
  sprintf(fc->value[num++], "%i", fluid_approx);
  // strcpy(fc->name[num], "input_verbose");
  // sprintf(fc->value[num++], "%i", 1);
  // strcpy(fc->name[num], "background_verbose");
  // sprintf(fc->value[num++], "%i", 1);

/* Extra available parameters if using the neutrino cosmology module */
#ifdef NEUTRINO_BACKGROUND
  /* Use the internal parameters */
  strcpy(fc->name[num], "T_cmb");
  sprintf(fc->value[num++], "%e", cosmo->T_CMB);
  strcpy(fc->name[num], "N_ncdm");
  sprintf(fc->value[num++], "%zu", cosmo->N_nu);

  float N_ur = 0.00707;  // effectively extra ultra-relativistic species

  /* This is debatable */
  message("We assumed that N_ur = %.5f to get N_eff = 3.046.", N_ur);
  strcpy(fc->name[num], "N_ur");
  sprintf(fc->value[num++], "%e", N_ur);

  /* Insert the neutrino masses for each species*/
  strcpy(fc->name[num], "m_ncdm");
  for (size_t i = 0; i < cosmo->N_nu; i++) {
    sprintf(fc->value[num] + strlen(fc->value[num]), "%e", cosmo->M_nu[i]);
    /* Add commas until the last species */
    if (i < e->cosmology->N_nu - 1) {
      sprintf(fc->value[num] + strlen(fc->value[num]), ",");
    }
  }
  num++;

  /* Insert number of momentum bins for each species */
  strcpy(fc->name[num], "Number of momentum bins");
  for (size_t i = 0; i < cosmo->N_nu; i++) {
    sprintf(fc->value[num] + strlen(fc->value[num]), "%i", nbins);
    /* Add commas until the last species */
    if (i < e->cosmology->N_nu - 1) {
      sprintf(fc->value[num] + strlen(fc->value[num]), ",");
    }
  }
  num++;

  /* Insert the quadrature strategy */
  strcpy(fc->name[num], "Quadrature strategy");
  for (size_t i = 0; i < cosmo->N_nu; i++) {
    sprintf(fc->value[num] + strlen(fc->value[num]), "%i", quadrature_strategy);
    /* Add commas until the last species */
    if (i < e->cosmology->N_nu - 1) {
      sprintf(fc->value[num] + strlen(fc->value[num]), ",");
    }
  }
  num++;
#else
  /* Infer neutrino mass from the simulation particles */
  float neutrino_mass_min = e->neutrino_mass_min;
  float neutrino_mass_max = e->neutrino_mass_max;

  if (neutrino_mass_min != neutrino_mass_max) {
    error(
        "Running without CLASS parameter file or neutrino background, \
        and with non-degenerate simulation particles. Cannot consistently \
        infer CLASS parameters.");
  } else {
    float m_ncdm = neutrino_mass_min;
    float N_ur = 0.00707;  // effectively extra ultra-relativistic species

    message("WARNING: no CLASS parameter file or neutrino background.");
    message("We inferred that M_nu = %.4f eV from the particles.", m_ncdm);
    message("We assumed that there are three degenerate neutrinos.");

    if (Omega_r != 0) {
      message("We inferred that Omega_g = %.3e (photon density).", Omega_r);
      strcpy(fc->name[num], "Omega_g");
      sprintf(fc->value[num++], "%e", Omega_r);
    }

    /* This is debatable */
    message("We assumed that N_ur = %.5f to get N_eff = 3.046.", N_ur);
    strcpy(fc->name[num], "N_ur");
    sprintf(fc->value[num++], "%e", N_ur);

    /* This is reasonable */
    strcpy(fc->name[num], "N_ncdm");
    sprintf(fc->value[num++], "%i", 3);
    strcpy(fc->name[num], "m_ncdm");
    sprintf(fc->value[num++], "%e,%e,%e", m_ncdm, m_ncdm, m_ncdm);
    strcpy(fc->name[num], "Number of momentum bins");
    sprintf(fc->value[num++], "%i,%i,%i", nbins, nbins, nbins);
    strcpy(fc->name[num], "Quadrature strategy");
    sprintf(fc->value[num++], "%i,%i,%i", quadrature_strategy,
            quadrature_strategy, quadrature_strategy);
  }
#endif
#else
  error("No CLASS library found. Cannot compute transfer functions.");
#endif
}
