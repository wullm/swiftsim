
/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2017 Stefan Arridge (stefan.arridge@durham.ac.uk)

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
#ifndef SWIFT_STAR_FORMATION_DENS_THRESH_H
#define SWIFT_STAR_FORMATION_NONE_H

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <float.h>

/* Local includes. */
#include "error.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "units.h"
#include "cell.h"
#include "hydro.h"
#include "stragglers.h"

/**
 * @brief Star Formation Properties
 */
struct star_formation_data {
      
  /*! Density threshold in internal units */
  float density_threshold;

  /*! Fraction of gas mass that is Hydrogen. Used to convert density threshold from hydrogen number density
   to internal units*/
  float hydrogen_mass_abundance;
 };

/**
 * @brief Creates a star particle if density of particle is above the density threshold
 *        
 * @param star_formation The proerties of the star formation model.
 * @param phys_const The physical constants in internal units.
 * @param p Pointer to the particle data.
 * @param c Pointer to the cell to which the particle belongs
 */
__attribute__((always_inline)) INLINE static void do_star_formation(
    const struct star_formation_data* restrict star_formation,
    const struct phys_const* restrict phys_const, struct stragglers* restrict stragglers,
    struct part* restrict p, struct cell* restrict c, const float dt) {

  /* Get particle density */
  const float rho = hydro_get_density(p);
  
  if (rho > star_formation->density_threshold){
    
    /* Create gpart with same properties as the gas particle's gpart */
    struct gpart new_gpart = *(p->gpart);

    /* Mark it as a star */
    new_gpart.type = swift_type_star;
    
    /* Add it to the straggler gpart array */
    struct gpart* new_gpart_pointer = stragglers_add_gpart(stragglers,&new_gpart);

    /* Create a link to it from the cell */
    struct g_straggler_link* new_glink = malloc(sizeof(struct g_straggler_link));

    new_glink->gp = new_gpart_pointer;
    new_glink->next = c->g_straggler_next;
    c->g_straggler_next = new_glink;
    c->straggler_gcount++;

    /* Create star with same properties as the gas particle */
    struct spart new_star;
    new_star.id = p->id;
    new_star.mass = p->mass;
    new_star.x[0] = p->x[0];
    new_star.x[1] = p->x[1];
    new_star.x[2] = p->x[2];
    new_star.v[0] = p->v[0];
    new_star.v[1] = p->v[1];
    new_star.v[2] = p->v[2];
    new_star.time_bin = p->time_bin;
    
    /* Make it link to the gpart we just created */
    new_star.gpart =new_gpart_pointer;

    /* Add it to the straggler spart array */
    struct spart* new_star_pointer = stragglers_add_spart(stragglers,&new_star);
  
    /* Create a link to it from the cell */
    struct s_straggler_link* new_slink = malloc(sizeof(struct s_straggler_link));

    new_slink->sp = new_star_pointer;
    new_slink->next = c->s_straggler_next;
    c->s_straggler_next = new_slink;
    c->straggler_scount++;

    /* For now we 'delete' the gas particle by setting  its 
       time bin to be very high so that it no longer interacts. */

    p->time_bin = 128;
    p->gpart->time_bin = 128;
  }
}



/**
 * @brief Initialises the star formation properties in the internal system
 * of units.
 *
 * Nothing to do here.
 *
 * @param parameter_file The parsed parameter file
 * @param phys_const Physical constants in internal units
 * @param us The current internal system of units
 * @param star_formation The star formation properties to initialize
 */
static INLINE void star_formation_init_backend(
    const struct swift_params* parameter_file,
     const struct UnitSystem* us,
    const struct phys_const* phys_const,
    struct star_formation_data* star_formation) {

  const float number_density_threshold_cgs =
    parser_get_param_double(parameter_file, "DensThreshStarFormation:number_density_threshold_cgs");

  star_formation->hydrogen_mass_abundance = parser_get_param_double(
      parameter_file, "DensThreshStarFormation:hydrogen_mass_abundance");

  /* Get the density thresholf in internal units */
  star_formation->density_threshold = number_density_threshold_cgs *
                                      units_cgs_conversion_factor(us, UNIT_CONV_VOLUME) *
                                      phys_const->const_proton_mass /
                                      star_formation->hydrogen_mass_abundance;
}

/**
 * @brief Prints the properties of the star formation model to stdout.
 *
 * @param star_formation The star formation properties.
 */
static INLINE void star_formation_print_backend(
    const struct star_formation_data* star_formation) {

  message("Star formation model is 'Density threshold star formation' with"
          "(density threshold, hydrogen_mass_abundance) "
          " = (%g, %g)",
          star_formation->density_threshold, star_formation->hydrogen_mass_abundance);
}

#endif /* SWIFT_STAR_FORMATION_NONE_H */
