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

/* Some standard headers. */
#include "../config.h"

/* Some standard headers */
#include <fenv.h>
#include <math.h>

/* Includes. */
#include "swift.h"

#if defined(NEUTRINO_BACKGROUND)

#define N_CHECK 20
#define TOLERANCE 1e-6

void test_params_init(struct swift_params *params, int testnr) {
  switch (testnr) {
    /* Three massless neutrinos */
    case 0:
      parser_init("", params);
      parser_set_param(params, "Cosmology:Omega_m:0.3075");
      parser_set_param(params, "Cosmology:Omega_lambda:0.6910");
      parser_set_param(params, "Cosmology:Omega_b:0.0486");
      parser_set_param(params, "Cosmology:h:0.6774");
      parser_set_param(params, "Cosmology:a_begin:0.1");
      parser_set_param(params, "Cosmology:a_end:1.0");
      parser_set_param(params, "Cosmology:Omega_g:5.2e-5");
      parser_set_param(params, "Cosmology:N_nu:3");
      parser_set_param(params, "Cosmology:N_eff:3.046");
      break;
    /* Three massive neutrinos */
    case 1:
      parser_init("", params);
      parser_set_param(params, "Cosmology:Omega_m:0.3075");
      parser_set_param(params, "Cosmology:Omega_lambda:0.6910");
      parser_set_param(params, "Cosmology:Omega_b:0.0486");
      parser_set_param(params, "Cosmology:h:0.6774");
      parser_set_param(params, "Cosmology:a_begin:0.1");
      parser_set_param(params, "Cosmology:a_end:1.0");
      parser_set_param(params, "Cosmology:Omega_g:5.2e-5");
      parser_set_param(params, "Cosmology:N_nu:3");
      parser_set_param(params, "Cosmology:N_eff:3.046");
      parser_set_param(params, "Cosmology:M_nu:0.02,0.03,0.04");
      break;
    /* Two massive neutrinos + one massless */
    case 2:
      parser_init("", params);
      parser_set_param(params, "Cosmology:Omega_m:0.3075");
      parser_set_param(params, "Cosmology:Omega_lambda:0.6910");
      parser_set_param(params, "Cosmology:Omega_b:0.0486");
      parser_set_param(params, "Cosmology:h:0.6774");
      parser_set_param(params, "Cosmology:a_begin:0.1");
      parser_set_param(params, "Cosmology:a_end:1.0");
      parser_set_param(params, "Cosmology:Omega_g:5.2e-5");
      parser_set_param(params, "Cosmology:N_nu:3");
      parser_set_param(params, "Cosmology:N_eff:3.046");
      parser_set_param(params, "Cosmology:M_nu:0,0.02,0.04");
      break;
    /* Three massless neutrinos - same as case 0 but with different
     * parametrisation */
    case 3:
      parser_init("", params);
      parser_set_param(params, "Cosmology:Omega_m:0.3075");
      parser_set_param(params, "Cosmology:Omega_lambda:0.6910");
      parser_set_param(params, "Cosmology:Omega_b:0.0486");
      parser_set_param(params, "Cosmology:Omega_r:8.7971982e-5");
      parser_set_param(params, "Cosmology:h:0.6774");
      parser_set_param(params, "Cosmology:a_begin:0.1");
      parser_set_param(params, "Cosmology:a_end:1.0");
      parser_set_param(params, "Cosmology:N_nu:0");
      parser_set_param(params, "Cosmology:N_eff:0");
      break;
    /* Three very massive neutrinos from very early times to today */
    case 4:
      parser_init("", params);
      parser_set_param(params, "Cosmology:Omega_m:0.3075");
      parser_set_param(params, "Cosmology:Omega_lambda:0.6910");
      parser_set_param(params, "Cosmology:Omega_b:0.0486");
      parser_set_param(params, "Cosmology:h:0.6774");
      parser_set_param(params, "Cosmology:a_begin:1e-2");
      parser_set_param(params, "Cosmology:a_end:1.0");
      parser_set_param(params, "Cosmology:Omega_g:5.2e-5");
      parser_set_param(params, "Cosmology:N_nu:3");
      parser_set_param(params, "Cosmology:N_eff:3.046");
      parser_set_param(params, "Cosmology:M_nu:0.1,0.2,0.3");
      break;
  }
}

int main(int argc, char *argv[]) {

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

  /* We expect an overflow on the indefinite integral a la 1 / exp(bignum) ,
   * which is fine, so do not hang on FP-exceptions.
   */
#ifdef HAVE_FE_ENABLE_EXCEPT
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  /* Test a number of different cosmologies */
  int N_cosmos = 4;
  double times1[N_CHECK];  // compare cosmology 0
  double times2[N_CHECK];  // with cosmology 3
  for (int testnr = 0; testnr < N_cosmos; testnr++) {
    message("Initialization of cosmology %i", testnr);

    /* pseudo initialization of params */
    struct swift_params params;
    test_params_init(&params, testnr);

    /* initialization of unit system */
    struct unit_system us;
    units_init_cgs(&us);

    /* initialization of phys_const */
    struct phys_const phys_const;
    phys_const_init(&us, &params, &phys_const);

    /* initialization of cosmo */
    struct cosmology cosmo;
    cosmology_init(&params, &us, &phys_const, &cosmo);

    message("Start checking computation...");

    for (int i = 0; i < N_CHECK; i++) {
      double a = 0.1 + 0.9 * i / (N_CHECK - 1.);
      /* Compute a(t(a)) and check if same results */
      double time = cosmology_get_time_since_big_bang(&cosmo, a);

      /* Store the value for cosmologies 0 and 3 to compare later */
      if (testnr == 0) {
        times1[i] = time;
      } else if (testnr == 3) {
        times2[i] = time;
      }

      // double my_a = cosmology_get_scale_factor_from_time(&cosmo, time); //from the integration branch
      double my_a = cosmology_get_scale_factor(&cosmo, time);

      /* check accuracy */
      double rel_err = (my_a - a) / a;
      message("Accuracy of %g at a=%g", rel_err, a);
      assert(fabs(rel_err) < TOLERANCE);
    }
    message("Everything seems fine with this cosmology.");

    cosmology_clean(&cosmo);
  }

  /* Compare cosmologies 0 and 3 */
  for (int j = 0; j < N_CHECK; j++) {
    /* check accuracy */
    double err = (times1[j] - times2[j]) / times1[j];
    message("Agreement of %g at step %i", err, j);
    assert(fabs(err) < TOLERANCE);
  }

  /* For the massive neutrino cosmology 4, check that it reproduces
   *  the relativistic (early) and non-relativistic (late) limits */

  message("Initialization...");

  /* pseudo initialization of params */
  struct swift_params params;
  test_params_init(&params, 4);

  /* initialization of unit system */
  struct unit_system us;
  units_init_cgs(&us);

  /* initialization of phys_const */
  struct phys_const phys_const;
  phys_const_init(&us, &params, &phys_const);

  /* initialization of cosmo */
  struct cosmology cosmo;
  cosmology_init(&params, &us, &phys_const, &cosmo);

  message("Start checking the fermion integration...");

  /* Relativistic limit */
  double Omega_nu_early = cosmology_get_neutrino_density_param(
      &cosmo, exp(cosmo.log_a_nutab_begin));
  double Omega_nu_r =
      cosmo.Omega_g * cosmo.N_eff * (7. / 8.) * pow(4. / 11., 4. / 3.);
  double err = (Omega_nu_early - Omega_nu_r) / Omega_nu_early;
  message("Accuracy of %g of the relativistic limit", err);
  assert(fabs(err) < TOLERANCE);

  /* Zeta function and other constants */
  const double zeta3 = 1.202056903159594;
  const double kb = phys_const.const_boltzmann_k;
  const double hbar = phys_const.const_planck_h / (2 * M_PI);
  const double cvel = phys_const.const_speed_light_c;
  const double eV = phys_const.const_electron_volt;

  /* Non-relativistic limit (limit not reached, so lower tolerance is fine)*/
  double Omega_nu_late =
      cosmology_get_neutrino_density_param(&cosmo, exp(cosmo.log_a_nutab_end));
  double Omega_nu_nr = 6 * zeta3 / (11 * M_PI * M_PI) *
                       pow(kb * cosmo.T_CMB, 3) / pow(cvel * hbar, 3) *
                       cosmo.M_nu_tot * eV /
                       (cosmo.critical_density_0 * cvel * cvel);
  double err2 = (Omega_nu_late - Omega_nu_nr) / Omega_nu_late;
  message("Accuracy of %g of the non-relativistic limit", err2);
  assert(fabs(err2) < 0.05);

  return 0;
}

#else

int main(int argc, char **argv) { return 0; }

#endif
