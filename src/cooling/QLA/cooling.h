/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_COOLING_QLA_H
#define SWIFT_COOLING_QLA_H

/**
 * @file src/cooling/QLA/cooling.h
 * @brief Quick Lyman alpha cooling (COLIBRE with primoridal Z) function
 * declarations
 */

/* Local includes. */
#include "cooling_struct.h"
#include "cooling_tables.h"

struct part;
struct xpart;
struct cosmology;
struct hydro_props;
struct entropy_floor_properties;
struct phys_const;
struct space;

void cooling_update(const struct cosmology *cosmo,
                    struct cooling_function_data *cooling, struct space *s);

void cooling_cool_part(const struct phys_const *phys_const,
                       const struct unit_system *us,
                       const struct cosmology *cosmo,
                       const struct hydro_props *hydro_properties,
                       const struct entropy_floor_properties *floor_props,
                       const struct cooling_function_data *cooling,
                       struct part *p, struct xpart *xp, const float dt,
                       const float dt_therm, const double time);

float cooling_timestep(const struct cooling_function_data *cooling,
                       const struct phys_const *phys_const,
                       const struct cosmology *cosmo,
                       const struct unit_system *us,
                       const struct hydro_props *hydro_props,
                       const struct part *p, const struct xpart *xp);

void cooling_first_init_part(const struct phys_const *phys_const,
                             const struct unit_system *us,
                             const struct hydro_props *hydro_props,
                             const struct cosmology *cosmo,
                             const struct cooling_function_data *cooling,
                             struct part *p, struct xpart *xp);

float cooling_get_temperature_from_gas(
    const struct phys_const *phys_const, const struct cosmology *cosmo,
    const struct cooling_function_data *cooling, const float rho_phys,
    const float XH, const float logZZsol, const float u_phys,
    const int HII_region);

float cooling_get_temperature(const struct phys_const *phys_const,
                              const struct hydro_props *hydro_props,
                              const struct unit_system *us,
                              const struct cosmology *cosmo,
                              const struct cooling_function_data *cooling,
                              const struct part *p, const struct xpart *xp);

float cooling_get_radiated_energy(const struct xpart *xp);

void cooling_split_part(struct part *p, struct xpart *xp, double n);

void cooling_Hydrogen_reionization(const struct cooling_function_data *cooling,
                                   const struct cosmology *cosmo,
                                   struct space *s);

void cooling_init_backend(struct swift_params *parameter_file,
                          const struct unit_system *us,
                          const struct phys_const *phys_const,
                          const struct hydro_props *hydro_props,
                          struct cooling_function_data *cooling);

void cooling_print_backend(const struct cooling_function_data *cooling);

void cooling_clean(struct cooling_function_data *data);

#endif /* SWIFT_COOLING_QLA_H */
