/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_CHEMISTRY_GEAR_H
#define SWIFT_CHEMISTRY_GEAR_H

/**
 * @file src/chemistry/GEAR/chemistry.h
 * @brief Empty infrastructure for the cases without chemistry function
 */

/* Some standard headers. */
#include <float.h>
#include <math.h>
#include <string.h>

/* Local includes. */
#include "chemistry_struct.h"
#include "error.h"
#include "hydro.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "units.h"

/**
 * @brief Compute the metal mass
 *
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 * @param data The global chemistry information.
 */
__attribute__((always_inline)) INLINE static float
chemistry_part_metal_mass(const struct part* restrict p,
			  const struct xpart* restrict xp) {

  return p->chemistry_data.metal_mass[CHEMISTRY_ELEMENT_COUNT - 1];
}

/**
 * @brief Compute the smoothed metal mass fraction
 *
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 * @param data The global chemistry information.
 */
__attribute__((always_inline)) INLINE static float
chemistry_part_smoothed_metal_mass_fraction(const struct part* restrict p,
			  const struct xpart* restrict xp) {

  return p->chemistry_data.smoothed_metal_mass_fraction[CHEMISTRY_ELEMENT_COUNT - 1];
}

/**
 * @brief Compute the metal mass fraction
 *
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 * @param data The global chemistry information.
 */
__attribute__((always_inline)) INLINE static float
chemistry_spart_metal_mass_fraction(const struct spart* restrict sp) {

  return sp->chemistry_data.metal_mass_fraction[CHEMISTRY_ELEMENT_COUNT - 1];
}

/**
 * @brief Copies the chemistry properties of the gas particle over to the
 * star particle.
 *
 * @param p the gas particles.
 * @param xp the additional properties of the gas particles.
 * @param sp the new created star particle with its properties.
 */
INLINE static void chemistry_copy_star_formation_properties(
    const struct part* p, const struct xpart* xp, struct spart* sp) {

  /* Store the chemistry struct in the star particle */
  for(int i = 0; i < CHEMISTRY_ELEMENT_COUNT; i++) {
    sp->chemistry_data.metal_mass_fraction[i] = p->chemistry_data.smoothed_metal_mass_fraction[i];
  }
}


/**
 * @brief Prints the properties of the chemistry model to stdout.
 *
 * @brief The #chemistry_global_data containing information about the current
 * model.
 */
static INLINE void chemistry_print_backend(
    const struct chemistry_global_data* data) {

  message("Chemistry function is 'Gear'.");
}

/**
 * @brief Read the solar abundances and scale with them the initial metallicities.
 *
 * @param parameter_file The parsed parameter file.
 * @param data The properties to initialise.
 */
static INLINE void chemistry_scale_initial_metallicities(struct swift_params* parameter_file,
							 struct chemistry_global_data *data) {
#ifdef HAVE_HDF5

  /* Get the yields table */
  char filename[DESCRIPTION_BUFFER_SIZE];
  parser_get_param_string(parameter_file, "GEARFeedback:YieldsTable", filename);

  /* Open file. */
  hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0) error("unable to open file %s.\n", filename);

  /* Open group. */
  hid_t group_id = H5Gopen(file_id, "Data", H5P_DEFAULT);
  if (group_id < 0) error("unable to open group Data.\n");

  /* Read the data */
  float *sol_ab = (float*) malloc(sizeof(float) * CHEMISTRY_ELEMENT_COUNT);
  io_read_array_attribute(group_id, "SolarMassAbundances", FLOAT, sol_ab, CHEMISTRY_ELEMENT_COUNT);

  /* Close group */
  hid_t status = H5Gclose(group_id);
  if (status < 0) error("error closing group.");

  /* Close file */
  status = H5Fclose(file_id);
  if (status < 0) error("error closing file.");


  /* Scale the initial metallicities */
  char txt[DESCRIPTION_BUFFER_SIZE] = "Scaling initial metallicities by:";
  for(int i = 0; i < CHEMISTRY_ELEMENT_COUNT; i++) {
    data->initial_metallicities[i] *= sol_ab[i];
    char tmp[10];
    sprintf(tmp, " %.2g", sol_ab[i]);
    strcat(txt, tmp);
  }

  message("%s", txt);
#else
  error("Cannot scale the solar abundances without HDF5");
#endif
}

/**
 * @brief Initialises the chemistry properties.
 *
 * Nothing to do here.
 *
 * @param parameter_file The parsed parameter file.
 * @param us The current internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param data The properties to initialise.
 */
static INLINE void chemistry_init_backend(struct swift_params* parameter_file,
                                          const struct unit_system* us,
                                          const struct phys_const* phys_const,
                                          struct chemistry_global_data* data) {

  /* read parameters */
  const float initial_metallicity = parser_get_param_float(
      parameter_file, "GEARChemistry:InitialMetallicity");

  /* Set the initial metallicities */
  for(int i = 0; i < CHEMISTRY_ELEMENT_COUNT; i++) {
      data->initial_metallicities[i] = initial_metallicity;
  }

  /* Check if need to scale the initial metallicity */
  const int scale_metallicity = parser_get_opt_param_int(
      parameter_file, "GEARChemistry:ScaleInitialMetallicity", 0);

  /* Scale the metallicities if required */
  if (scale_metallicity) {
    chemistry_scale_initial_metallicities(parameter_file, data);
  }
  
}

/**
 * @brief Prepares a particle for the smooth metal calculation.
 *
 * Zeroes all the relevant arrays in preparation for the sums taking place in
 * the various smooth metallicity tasks
 *
 * @param p The particle to act upon
 * @param cd #chemistry_global_data containing chemistry informations.
 */
__attribute__((always_inline)) INLINE static void chemistry_init_part(
    struct part* restrict p, const struct chemistry_global_data* cd) {

  struct chemistry_part_data* cpd = &p->chemistry_data;

  for (int i = 0; i < CHEMISTRY_ELEMENT_COUNT; i++) {
    /* Reset the smoothed metallicity */
    cpd->smoothed_metal_mass_fraction[i] = 0.f;

    /* Convert the total mass into mass fraction */
    cpd->metal_mass_fraction[i] = cpd->metal_mass[i] / p->mass;
  }
}

/**
 * @brief Finishes the smooth metal calculation.
 *
 * Multiplies the smoothed metallicity and number of neighbours by the
 * appropiate constants and add the self-contribution term.
 *
 * This function requires the #hydro_end_density to have been called.
 *
 * @param p The particle to act upon.
 * @param cd #chemistry_global_data containing chemistry informations.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void chemistry_end_density(
    struct part* restrict p, const struct chemistry_global_data* cd,
    const struct cosmology* cosmo) {

  /* Some smoothing length multiples. */
  const float h = p->h;
  const float h_inv = 1.0f / h;                       /* 1/h */
  const float factor = pow_dimension(h_inv) / p->rho; /* 1 / h^d * rho */
  const float m = p->mass;

  struct chemistry_part_data* cpd = &p->chemistry_data;

  for (int i = 0; i < CHEMISTRY_ELEMENT_COUNT; i++) {
    /* Final operation on the density (add self-contribution). */
    cpd->smoothed_metal_mass_fraction[i] +=
        m * cpd->metal_mass_fraction[i] * kernel_root;

    /* Finish the calculation by inserting the missing h-factors */
    cpd->smoothed_metal_mass_fraction[i] *= factor;

    /* Convert the mass fraction into a total mass */
    cpd->metal_mass[i] = m * cpd->metal_mass_fraction[i];
  }
}

/**
 * @brief Updates to the chemistry data after the hydro force loop.
 *
 * @param p The particle to act upon.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void chemistry_end_force(
    struct part* restrict p, const struct cosmology* cosmo) {}

/**
 * @brief Sets all particle fields to sensible values when the #part has 0 ngbs.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 * @param cd #chemistry_global_data containing chemistry informations.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void
chemistry_part_has_no_neighbours(struct part* restrict p,
                                 struct xpart* restrict xp,
                                 const struct chemistry_global_data* cd,
                                 const struct cosmology* cosmo) {

  /* Set the smoothed fractions with the non smoothed fractions */
  for(int i = 0; i < CHEMISTRY_ELEMENT_COUNT; i++) {
    p->chemistry_data.smoothed_metal_mass_fraction[i] = p->chemistry_data.metal_mass_fraction[i];
    p->chemistry_data.metal_mass[i] = p->chemistry_data.metal_mass_fraction[i] * p->mass;
  }
}

/**
 * @brief Computes the chemistry-related time-step constraint.
 *
 * No constraints in the GEAR model (no diffusion) --> FLT_MAX
 *
 * @param phys_const The physical constants in internal units.
 * @param cosmo The current cosmological model.
 * @param us The internal system of units.
 * @param hydro_props The properties of the hydro scheme.
 * @param cd The global properties of the chemistry scheme.
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float chemistry_timestep(
    const struct phys_const* restrict phys_const,
    const struct cosmology* restrict cosmo,
    const struct unit_system* restrict us,
    const struct hydro_props* hydro_props,
    const struct chemistry_global_data* cd, const struct part* restrict p) {
  return FLT_MAX;
}

/**
 * @brief Sets the chemistry properties of the (x-)particles to a valid start
 * state.
 *
 * Nothing to do here.
 *
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 * @param data The global chemistry information.
 */
__attribute__((always_inline)) INLINE static void chemistry_first_init_part(
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo,
    const struct chemistry_global_data* data, struct part* restrict p,
    struct xpart* restrict xp) {

  for(int i = 0; i < CHEMISTRY_ELEMENT_COUNT; i++) {
    p->chemistry_data.metal_mass[i] = data->initial_metallicities[i] * p->mass;
  }

  chemistry_init_part(p, data);

}

/**
 * @brief Sets the chemistry properties of the sparticles to a valid start
 * state.
 *
 * @param data The global chemistry information.
 * @param sp Pointer to the sparticle data.
 */
__attribute__((always_inline)) INLINE static void chemistry_first_init_spart(
    const struct chemistry_global_data* data, struct spart* restrict sp) {

  for(int i = 0; i < CHEMISTRY_ELEMENT_COUNT; i++) {
    sp->chemistry_data.metal_mass_fraction[i] = data->initial_metallicities[i];
  }
}

#endif /* SWIFT_CHEMISTRY_GEAR_H */
