/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Loic Hausamman (loic.hausammann@epfl.ch)
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

#define N_CHECK 200
#define TOLERANCE 1e-6

void test_params_init(struct swift_params *params) {
  parser_init("", params);
  parser_set_param(params, "Cosmology:Omega_m:0.3075");
  parser_set_param(params, "Cosmology:Omega_lambda:0.6910");
  parser_set_param(params, "Cosmology:Omega_b:0.0486");
  parser_set_param(params, "Cosmology:h:0.6774");
  parser_set_param(params, "Cosmology:a_begin:0.01");
  parser_set_param(params, "Cosmology:a_end:1.0");
}

int main(int argc, char *argv[]) {

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

/* Choke on FP-exceptions */
#ifdef HAVE_FE_ENABLE_EXCEPT
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  message("Initialization...");

  /* pseudo initialization of params */
  struct swift_params params;
  test_params_init(&params);

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
    double a = 0.01 + 0.99 * i / (N_CHECK - 1.);

    /* Compute a(t(a)) and check if same results */
    double time = cosmology_get_time_since_big_bang(&cosmo, a);
    double my_a = cosmology_get_scale_factor_from_time(&cosmo, time);

    /* check accuracy */
    double rel_err = (my_a - a) / a;
    message("Accuracy of %g at a=%g", rel_err, a);
    assert(fabs(rel_err) < TOLERANCE);
  }

  message("Everything seems fine with cosmology.");

  cosmology_clean(&cosmo);
  return 0;
}
