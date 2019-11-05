/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_STELLAR_EVOLUTION_STRUCT_GEAR_H
#define SWIFT_STELLAR_EVOLUTION_STRUCT_GEAR_H

#include "interpolation.h"

/* Number of different type of companion.
   If changed, the IO needs to be updated.
 */
#define GEAR_NUMBER_TYPE_OF_COMPANION 2
#define GEAR_LABELS_SIZE 10

/**
 * @brief Model for the initial mass function.
 *
 * Describe a model such as Kroupa 2001:
 *
 * f(m) = coef[i] * pow(m, exp[i])
 */
struct initial_mass_function {

  /*! Mass limits between IMF parts (n_parts + 1 elements). */
  double *mass_limits;

  /*! Exponent of each IMF parts (n_parts elements). */
  double *exp;

  /*! Coefficient of each IMF parts (n_parts elements). */
  double *coef;

  /*! Number of parts in the function. */
  int n_parts;

  /*! Minimal mass contained in mass_limits, copied for more clarity. */
  double mass_min;

  /*! Maximal mass contained in mass_limits, copied for more clarity. */
  double mass_max;

};

/**
 * @brief Model for the stellar lifetime.
 */
struct lifetime {

  /*! Coefficients for the log10(m)^2 term */
  double quadratic[3];

  /*! Coefficients for the log10(m) term */
  double linear[3];

  /*! Coefficients for the constant term */
  double constant[3];

  /*! Factor for the mass unit change. */
  double log_unit_mass;
};

/**
 * @brief Model for SNIa.
 */
struct supernovae_ia {
  /*! Yields TODO more comment */
  struct {

    /*! Mass of each element ejected by a single supernovae */
    double data[CHEMISTRY_ELEMENT_COUNT];

  } yields;

  /*! White dwarf's mass */
  double mass_white_dwarf;

  /*! Minimal mass of the progenitor */
  double mass_min_progenitor;

  /*! Maximal mass of the progenitor */
  double mass_max_progenitor;

  /*! coefficient of the initial mass function for progenitor divided by progenitor_exponent */
  double progenitor_coef_exp;

  /*! exponent of the initial mass function for progenitor */
  double progenitor_exponent;

  /*! exponent of the initial mass function for binaries */
  double companion_exponent;

  struct {
    /*! Initial mass function's coeffcients */
    double coef;

    /*! Maximal mass of the companion */
    double mass_max;

    /*! Minimal mass of the companion */
    double mass_min;
  } companion[GEAR_NUMBER_TYPE_OF_COMPANION];


};

/**
 * @brief Model for SNII.
 */
struct supernovae_ii {

  /*! Integrated yields TODO more comment: Currently it is the integrated mass fraction */
  struct interpolation_1d integrated_yields[CHEMISTRY_ELEMENT_COUNT];

  /*! Integrated mass ejected (processed): Currently it is the mass fraction */
  struct interpolation_1d integrated_ejected_mass_processed;

  /*! Integrated mass ejected (non processed): Currently it is the mass fraction */
  struct interpolation_1d integrated_ejected_mass;

  /*! Minimal mass for a SNII */
  double mass_min;

  /*! Maximal mass for a SNII */
  double mass_max;

  /*! exponent of the IMF */
  double exponent;

  /*! coefficient of the IMF over the exponent */
  double coef_exp;
  
};

/**
 * @brief The complete stellar model.
 */
struct stellar_model {

  /*! Name of the different elements */
  char elements_name[CHEMISTRY_ELEMENT_COUNT * GEAR_LABELS_SIZE];
  
  /*! The initial mass function */
  struct initial_mass_function imf;

  /*! The stellar lifetime */
  struct lifetime lifetime;

  /*! The supernovae type Ia */
  struct supernovae_ia snia;

  /*! The supernovae type II */
  struct supernovae_ii snii;

};



#endif // SWIFT_STELLAR_EVOLUTION_STRUCT_GEAR_H
