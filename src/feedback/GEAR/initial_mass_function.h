/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_INITIAL_MASS_FUNCTION_GEAR_H
#define SWIFT_INITIAL_MASS_FUNCTION_GEAR_H

#include "hdf5_functions.h"
#include "stellar_evolution_struct.h"

/**
 * @brief Get the IMF exponent in between mass_min and mass_max.
 */
__attribute__((always_inline)) INLINE static float initial_mass_function_get_exponent(
    const struct initial_mass_function* imf, float mass_min, float mass_max) {

#ifdef SWIFT_DEBUG_CHECKS
  if (mass_max > imf->mass_max)
    error("Cannot have mass larger than the largest one in the IMF");
  if (mass_min < imf->mass_min)
    error("Cannot have mass smaller than the smallest one in the IMF");
  if (mass_max < mass_min)
    error("Cannot have mass_min larger than mass_max");
#endif
  
  for(int i = 0; i < imf->n_parts; i++) {

    /* Check if in the correct part of the IMF */
    if (mass_min < imf->mass_limits[i+1]) {

      /* Check if in only one segment */
      if (mass_max > imf->mass_limits[i+1]) {
	error("Cannot get a single exponent for the interval [%g, %g]", mass_min, mass_max);
      }

      return imf->exp[i];
    }
  }

  error("Masses outside IMF ranges");

  return -1;
}

/** @brief Print the initial mass function */
__attribute__((always_inline)) INLINE static void initial_mass_function_print(
    const struct initial_mass_function* imf) {

  message("Number of parts: %i", imf->n_parts);
  message("Mass interval: [%g, %g]", imf->mass_min, imf->mass_max);
  for(int i = 0; i < imf->n_parts; i++) {
    message("[%g, %g]: %g * m^{%g}", imf->mass_limits[i], imf->mass_limits[i+1],
	    imf->coef[i], imf->exp[i]);
  }
}

/**
 * @brief Integrate the #interpolation_1d data with the initial mass function.
 *
 * The x are supposed to be linear in log.
 *
 * @param imf The #initial_mass_function.
 * @param interp The #interpolation_1d.
 */
__attribute__((always_inline)) INLINE static void initial_mass_function_integrate(
    const struct initial_mass_function* imf, struct interpolation_1d *interp) {

  /* Index in the data */
  int j = 1;
  const float mass_min = pow(10, interp->xmin);
  const float mass_max = pow(10, interp->xmin + (interp->N - 1) * interp->dx);

  float m = mass_min;

  float *tmp = (float *) malloc(sizeof(float) * interp->N);

  /* Set lower limit */
  tmp[0] = 0;
  for(int i = 0; i < imf->n_parts; i++) {

    /* Check if already in the correct part */
    if (mass_min > imf->mass_limits[i+1]) {
      continue;
    }

    /* Check if already above the maximal mass */
    if (mass_max < imf->mass_limits[i]) {
      break;
    }

    /* Integrate the data */
    while (m < imf->mass_limits[i+1] && j < interp->N) {

      /* Compute the masses */
      const float log_m1 = interp->xmin + (j - 1) * interp->dx;
      const float m1 = pow(10, log_m1);
      const float log_m2 = interp->xmin + j * interp->dx;
      const float m2 = pow(10, log_m2);
      const float dm = m2 - m1;
      const float imf_1 = imf->coef[i] * pow(m1, imf->exp[i]);

      /* Get the imf of the upper limit  */
      float imf_2;
      if (m2 > imf->mass_limits[i+1]) {
      	imf_2 = imf->coef[i+1] * pow(m2, imf->exp[i+1]);
      }
      else {
	imf_2 = imf->coef[i] * pow(m2, imf->exp[i]);
      }

      /* Compute the integral */
      tmp[j] = tmp[j-1] + 0.5 * (imf_1 * interp->data[j-1] + imf_2 * interp->data[j]) * dm;

      /* Update j and m */
      j += 1;
      m = m2;
    }
  }

  /* The rest is extrapolated with 0 */
  for(int k = j; k < interp->N; k++) {
    tmp[k] = tmp[k-1];
  }

  /* Copy temporary array */
  memcpy(interp->data, tmp, interp->N * sizeof(float));

  /* Update the boundary conditions */
  interp->boundary_condition = boundary_condition_zero_const;

  /* clean everything */
  free(tmp);
}


/**
 * @brief Get the IMF coefficient in between mass_min and mass_max.
 *
 * @param imf The #initial_mass_function.
 * @param mass_min The minimal mass of the requested interval.
 * @param mass_max The maximal mass of the requested interval.
 *
 * @return The imf's coefficient of the interval.
 */
__attribute__((always_inline)) INLINE static float initial_mass_function_get_coefficient(
    const struct initial_mass_function* imf, float mass_min, float mass_max) {

  for(int i = 0; i < imf->n_parts; i++) {

    /* Check if in the correct part of the IMF */
    if (mass_min < imf->mass_limits[i+1]) {

      /* Check if in only one segment */
      if (mass_max > imf->mass_limits[i+1]) {
	error("Cannot get a single coefficient for the interval [%g, %g]", mass_min, mass_max);
      }

      return imf->coef[i];
    }
  }
  
  error("Masses outside IMF ranges");

  return -1;
}


/**
 * @brief Compute the mass fraction of the initial mass function.
 *
 * @param imf The #initial_mass_function.
 * @param m The mass to evaluate.
 *
 * @return The mass fraction.
 */
__attribute__((always_inline)) INLINE static float initial_mass_function_get_imf(
    const struct initial_mass_function *imf, float m) {

#ifdef SWIFT_DEBUG_CHECKS
  if (m > imf->mass_max || m < imf->mass_min)
    error("Mass below or above limits expecting %g < %g < %g.",
	  imf->mass_min, m, imf->mass_max);
#endif

  for(int i = 0; i < imf->n_parts; i++) {
    if (m <= imf->mass_limits[i+1]) {
      return imf->coef[i] * pow(m, imf->exp[i]);
    }
  }

  error("Failed to find correct function part: %g larger than mass max %g.",
	m, imf->mass_max);
};

/**
 * @brief Compute the fraction of stars (in number) in a given range of mass.
 *
 * @param imf The #initial_mass_function.
 * @param m1 The lowest mass.
 * @param m2 The largest mass.
 *
 * @return The fraction stars in the interval (in number).
 */
__attribute__((always_inline)) INLINE static float initial_mass_function_get_number(
    const struct initial_mass_function *imf, float m1, float m2) {
  error("This has not been tested. Need to check the units");
#ifdef SWIFT_DEBUG_CHECKS
  if (m1 > imf->mass_max || m1 < imf->mass_min)
    error("Mass 1 below or above limits expecting %g < %g < %g.",
	  imf->mass_min, m1, imf->mass_max);

  if (m2 > imf->mass_max || m2 < imf->mass_min)
    error("Mass 2 below or above limits expecting %g < %g < %g.",
	  imf->mass_min, m2, imf->mass_max);
  
  if (m1 > m2)
    error("Mass 1 larger than mass 2 %g > %g.", m1, m2);
#endif

  float n = 0.;

  for(int i = 0; i < imf->n_parts; i++) {
    /* Are we above the lowest mass? */
    if (m1 > imf->mass_limits[i+1])
      continue;

    /* Are we after the above the largest mass? */
    if (m2 < imf->mass_limits[i+1])
      break;

    /* Get integral limits. */
    const float mass_min = (m1 > imf->mass_limits[i]) ? m1: imf->mass_limits[i];
    const float mass_max = (m2 < imf->mass_limits[i+1]) ? m2: imf->mass_limits[i+1];
    const float exp = imf->exp[i];

    /* Compute the contribution of the current part. */
    n += imf->coef[i] * (pow(mass_max, exp) - pow(mass_min, exp)) / exp;
  }

  return n;
};

/**
 * @brief Compute the fraction of stars (in mass) in a given range of mass.
 *
 * @param imf The #initial_mass_function.
 * @param m1 The lowest mass.
 * @param m2 The largest mass.
 *
 * @return The fraction stars in the interval (in mass).
 */
__attribute__((always_inline)) INLINE static float initial_mass_function_get_imf_mass(
    const struct initial_mass_function *imf, float m1, float m2) {
  error("This has not been tested. Need to check the units");
#ifdef SWIFT_DEBUG_CHECKS
  if (m1 > imf->mass_max || m1 < imf->mass_min)
    error("Mass 1 below or above limits expecting %g < %g < %g.",
	  imf->mass_min, m1, imf->mass_max);

  if (m2 > imf->mass_max || m2 < imf->mass_min)
    error("Mass 2 below or above limits expecting %g < %g < %g.",
	  imf->mass_min, m2, imf->mass_max);
  
  if (m1 > m2)
    error("Mass 1 larger than mass 2 %g > %g.", m1, m2);
#endif

  float mass = 0.;

  for(int i = 0; i < imf->n_parts; i++) {
    /* Are we above the lowest mass? */
    if (m1 > imf->mass_limits[i+1])
      continue;

    /* Are we after the above the largest mass? */
    if (m2 < imf->mass_limits[i+1])
      break;

    /* Get integral limits. */
    const float mass_min = (m1 > imf->mass_limits[i]) ? m1: imf->mass_limits[i];
    const float mass_max = (m2 < imf->mass_limits[i+1]) ? m2: imf->mass_limits[i+1];
    const float exp = imf->exp[i] + 1.;

    /* Compute the contribution of the current part. */
    mass += imf->coef[i] * (pow(mass_max, exp) - pow(mass_min, exp)) / exp;
  }

  return mass;
};



/**
 * @brief Compute the coefficients of the initial mass function.
 *
 * @param imf The #initial_mass_function.
 */
__attribute__((always_inline)) INLINE static void initial_mass_function_compute_coefficients(
    struct initial_mass_function *imf) {

  /* Allocate memory */
  if ((imf->coef = (float *) malloc(sizeof(float) * imf->n_parts)) == NULL)
    error("Failed to allocate the IMF coefficients.");

  /* Suppose that the first coefficients is 1 (will be corrected later) */
  imf->coef[0] = 1.;

  /* Use the criterion of continuity for the IMF */
  for(int i = 1; i < imf->n_parts; i++) {
    float exp = imf->exp[i-1] - imf->exp[i];
    imf->coef[i] = imf->coef[i-1] * pow(imf->mass_limits[i], exp);
  }

  /* Use the criterion on the integral = 1 */
  float integral = 0;
  for(int i = 0; i < imf->n_parts; i++) {
    const float exp = imf->exp[i] + 1.;
    const float m_i = pow(imf->mass_limits[i], exp);
    const float m_i1 = pow(imf->mass_limits[i+1], exp);
    integral += imf->coef[i] * (m_i1 - m_i)/ exp;
  }

  /* Normalize the coefficients (fix initial supposition) */
  for(int i = 0; i < imf->n_parts; i++) {
    imf->coef[i] /= integral;
  }

}

/**
 * @brief Reads the initial mass function parameters from the tables.
 *
 * @param imf The #initial_mass_function.
 * @param params The #swift_params.
 */
__attribute__((always_inline)) INLINE static void initial_mass_function_read_from_table(
    struct initial_mass_function* imf, struct swift_params* params) {

  hid_t file_id, group_id;

  /* Open IMF group */
  h5_open_group(params, "Data/IMF", &file_id, &group_id);

  /* Read number of parts */
  io_read_attribute(group_id, "n", INT, &imf->n_parts);

  /* The tables have a different definition of n */
  imf->n_parts += 1;

  /* Allocate the memory for the exponents */
  if ((imf->exp = (float *)malloc(sizeof(float) * imf->n_parts)) == NULL)
    error("Failed to allocate the IMF exponents.");

  /* Read the exponents */
  io_read_array_attribute(group_id, "as", FLOAT, imf->exp, imf->n_parts);

  /* Allocate the memory for the temporary mass limits */
  if ((imf->mass_limits = (float *)malloc(sizeof(float) * (imf->n_parts + 1))) == NULL)
    error("Failed to allocate the IMF masses.");

  /* Read the mass limits */
  io_read_array_attribute(group_id, "ms", FLOAT, imf->mass_limits, imf->n_parts - 1);

  /* Copy the data (need to shift for mass_min) */
  for(int i = imf->n_parts - 1; i > 0; i--) {
    imf->mass_limits[i] = imf->mass_limits[i-1];
  }

  /* Read the minimal mass limit */
  io_read_attribute(group_id, "Mmin", FLOAT, &imf->mass_limits[0]);
  
  /* Read the maximal mass limit */
  io_read_attribute(group_id, "Mmax", FLOAT, &imf->mass_limits[imf->n_parts]);

  /* Close everything */
  h5_close_group(file_id, group_id);
}

/**
 * @brief Reads the parameters file and if required overwrites the parameters found in the yields table.
 *
 * @param imf The #initial_mass_function.
 * @param params The #swift_params.
 */
__attribute__((always_inline)) INLINE static void initial_mass_function_read_from_params(
    struct initial_mass_function* imf, struct swift_params* params) {

  /* Read the number of elements */
  const int n_parts = parser_get_opt_param_int(params, "GEARInitialMassFunction:number_function_part", imf->n_parts);

  const int n_parts_changed = n_parts != imf->n_parts;
  imf->n_parts = n_parts;

  /* Reallocate the exponent memory */
  if (n_parts_changed) {
    free(imf->exp);
    if ((imf->exp = (float *)malloc(sizeof(float) * imf->n_parts)) == NULL)
      error("Failed to allocate the IMF exponents.");
  }

  /* Read the exponents */
  const char *exponent_name = "GEARInitialMassFunction:exponents";
  if (n_parts_changed) {
    parser_get_param_float_array(params, exponent_name, imf->n_parts, imf->exp);
  }
  else {
    parser_get_opt_param_float_array(params, exponent_name, imf->n_parts, imf->exp);
  }

  /* Reallocate the mass limits memory */
  if (n_parts_changed) {
    free(imf->mass_limits);
    if ((imf->mass_limits = (float *)malloc(sizeof(float) * (imf->n_parts + 1))) ==
	NULL)
      error("Failed to allocate the IMF masses.");
  }

  /* Read the mass limits */
  const char *mass_limits_name = "GEARInitialMassFunction:mass_limits_msun";
  if (n_parts_changed) {
    parser_get_param_float_array(params, mass_limits_name,
				 imf->n_parts + 1, imf->mass_limits);
  }
  else {
    parser_get_opt_param_float_array(params, mass_limits_name,
				     imf->n_parts + 1, imf->mass_limits);
  }
}

/**
 * @brief Initialize the initial mass function.
 *
 * @param imf The #initial_mass_function.
 * @param phys_const The #phys_const.
 * @param us The #unit_system.
 * @param params The #swift_params.
 */
__attribute__((always_inline)) INLINE static void initial_mass_function_init(
    struct initial_mass_function* imf, const struct phys_const* phys_const,
    const struct unit_system* us, struct swift_params* params) {

  /* Read the parameters from the yields table */
  initial_mass_function_read_from_table(imf, params);

  /* Overwrites the parameters if found in the params file */
  initial_mass_function_read_from_params(imf, params);

  /* change the mass limits to internal units */
  for(int i = 0; i < imf->n_parts + 1; i++) {
    imf->mass_limits[i] *= phys_const->const_solar_mass;
  }

  /* Write the masses in the correct attributes */
  imf->mass_min = imf->mass_limits[0];
  imf->mass_max = imf->mass_limits[imf->n_parts];

  /* Compute the coefficients */
  initial_mass_function_compute_coefficients(imf);

}

#endif // SWIFT_INITIAL_MASS_FUNCTION_GEAR_H
