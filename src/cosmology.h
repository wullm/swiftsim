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
#ifndef SWIFT_COSMOLOGY_H
#define SWIFT_COSMOLOGY_H

/* Config parameters. */
#include "../config.h"

#include "parser.h"
#include "physical_constants.h"
#include "timeline.h"
#include "units.h"

#if defined(NEUTRINO_BACKGROUND)
#include "neutrino/cosmology.h"
#else
#include "cosmology_default.h"
#endif /* defined(NEUTRINO_BACKGROUND) */

void cosmology_update(struct cosmology *c, const struct phys_const *phys_const,
                      integertime_t ti_current);

double cosmology_get_drift_factor(const struct cosmology *cosmo,
                                  const integertime_t ti_start,
                                  const integertime_t ti_end);
double cosmology_get_grav_kick_factor(const struct cosmology *cosmo,
                                      const integertime_t ti_start,
                                      const integertime_t ti_end);
double cosmology_get_hydro_kick_factor(const struct cosmology *cosmo,
                                       const integertime_t ti_start,
                                       const integertime_t ti_end);
double cosmology_get_therm_kick_factor(const struct cosmology *cosmo,
                                       const integertime_t ti_start,
                                       const integertime_t ti_end);
double cosmology_get_corr_kick_factor(const struct cosmology *cosmo,
                                      const integertime_t ti_start,
                                      const integertime_t ti_end);
double cosmology_get_delta_time(const struct cosmology *c,
                                const integertime_t ti_start,
                                const integertime_t ti_end);

double cosmology_get_delta_time_from_scale_factors(const struct cosmology *c,
                                                   const double a_start,
                                                   const double a_end);

double cosmology_get_scale_factor(const struct cosmology *cosmo,
                                  integertime_t ti);
double cosmology_get_scale_factor_from_time(const struct cosmology *cosmo,
                                            double t);

double cosmology_get_time_since_big_bang(const struct cosmology *c, double a);
double cosmology_get_conformal_time(const struct cosmology *c, double a);
double cosmology_get_scale_factor_from_conformal_time(const struct cosmology *c,
                                            double t);
void cosmology_init(struct swift_params *params, const struct unit_system *us,
                    const struct phys_const *phys_const, struct cosmology *c);

void cosmology_init_no_cosmo(struct cosmology *c);

void cosmology_print(const struct cosmology *c);
void cosmology_clean(struct cosmology *c);

#ifdef HAVE_HDF5
void cosmology_write_model(hid_t h_grp, const struct cosmology *c);
#endif

/* Dump/restore. */
void cosmology_struct_dump(const struct cosmology *cosmology, FILE *stream);
void cosmology_struct_restore(int enabled, struct cosmology *cosmology,
                              const struct phys_const *phys_const,
                              FILE *stream);

#endif /* SWIFT_COSMOLOGY_H */
