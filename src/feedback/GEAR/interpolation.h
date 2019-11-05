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
#ifndef SWIFT_GEAR_INTERPOLATION_H
#define SWIFT_GEAR_INTERPOLATION_H

enum interpolate_boundary_condition {
  /* No extrapolation => raise errors */
  boundary_condition_error,

  /* Zero as boundary conditions */
  boundary_condition_zero,

  /* Zero (left boundary) and constant (right boundary) boundary conditions */
  boundary_condition_zero_const,
};

struct interpolation_1d {
  /* Data to interpolate */
  double *data;

  /* Minimal x */
  double xmin;

  /* Step size between x points */
  double dx;

  /* Number of element in the data */
  int N;

  /* Type of boundary conditions. */
  enum interpolate_boundary_condition boundary_condition;
};

/**
 * @brief Initialize the #interpolation_1d.
 *
 * Assumes x are linear in log.
 *
 * @params interp The #interpolation_1d.
 * @params xmin Minimal value of x (in log).
 * @params xmax Maximal value of x (in log).
 * @params N Requested number of values.
 * @params log_data_xmin The minimal value of the data (in log).
 * @params step_size The size of the x steps (in log).
 * @params N_data The number of element in the data.
 * @params data The data to interpolate (y).
 * @params N The number of element in data.
 * @params boundary_condition The type of #interpolate_boundary_condition.
 */
__attribute__((always_inline)) static INLINE void interpolate_1d_init(
    struct interpolation_1d *interp, double xmin, double xmax,
    int N, double log_data_xmin, double step_size, int N_data,
    const double *data, enum interpolate_boundary_condition boundary_condition) {

  /* Save the variables */
  interp->N = N;
  interp->xmin = xmin;
  interp->dx = (xmax - xmin) / (N - 1);
  interp->boundary_condition = boundary_condition;

  /* Allocate the memory */
  interp->data = malloc(sizeof(double) * N);
  if (interp->data == NULL)
    error("Failed to allocate memory for the interpolation");

  /* Interpolate the data */
  for(int i = 0; i < N; i++) {
    const double log_x = xmin + i * interp->dx;
    const double x_j = (log_x - log_data_xmin) / step_size;

    /* Check boundaries */
    if (x_j < 0) {
      switch (boundary_condition) {
        case boundary_condition_error:
	  error("Cannot extrapolate");
	  break;
        case boundary_condition_zero:
	  interp->data[i] = 0;
	  break;
        case boundary_condition_zero_const:
	  interp->data[i] = 0;
	  break;
        default:
	  error("Interpolation type not implemented");
      }
      continue;
    }
    else if (x_j >= N_data) {
      switch (boundary_condition) {
        case boundary_condition_error:
	  error("Cannot extrapolate");
	  break;
        case boundary_condition_zero:
	  interp->data[i] = 0;
	  break;
        case boundary_condition_zero_const:
	  interp->data[i] = interp->data[i-1];
	  break;
        default:
	  error("Interpolation type not implemented");
      }
      continue;
    }

    /* Interpolate i */
    const int j = x_j;
    const double f = x_j - j;
    interp->data[i] = (1. - f) * data[j] + f * data[j+1];
  }

}

/**
 * @brief Interpolate the data.
 *
 * @params interp The #interpolation_1d.
 * @params x The x value where to interpolate.
 *
 * @return The interpolated value y.
 */
__attribute__((always_inline)) static INLINE double interpolate_1d(
    const struct interpolation_1d *interp, double x) {

  /* Find indice */
  const double i = (x - interp->xmin) / interp->dx;
  const int idx = i;
  const double dx = i - idx;

  /* Should we extrapolate? */
  if (i < 0) {
    switch (interp->boundary_condition) {
      case boundary_condition_error:
	error("Cannot extrapolate");
	break;
      case boundary_condition_zero:
      case boundary_condition_zero_const:
	return 0;
      default:
	error("Interpolation type not implemented");
    }
  }
  else if (i >= interp->N - 1) {
    switch (interp->boundary_condition) {
      case boundary_condition_error:
	error("Cannot extrapolate");
	break;
      case boundary_condition_zero:
	return 0;
      case boundary_condition_zero_const:
	return interp->data[interp->N-1];
      default:
	error("Interpolation type not implemented");
    }
  }

  /* interpolate */
  return interp->data[idx] * (1. - dx) +
    interp->data[idx+1] * dx;
  
}

/**
 * @brief Print the data.
 *
 * @params interp The #interpolation_1d.
 */
__attribute__((always_inline)) static INLINE void interpolate_1d_print(
    const struct interpolation_1d *interp) {

  message("Interpolation between %g and %g", interp->xmin,
	  interp->xmin + interp->dx * interp->N);

  message("Contains %i values and use the boundary condition %i",
	  interp->N, interp->boundary_condition);

  /* Print values */
  for(int i = 0; i < interp->N; i++) {
    double x = interp->xmin + i * interp->dx;
    message("%.2g: %g", x, interp->data[i]);
  }
}



/**
 * @brief Cleanup the #interpolation_1d structure.
 *
 * @params interp The #interpolation_1d.
 */
__attribute__((always_inline)) static INLINE void interpolate_1d_free(
    struct interpolation_1d *interp) {

  /* Free the allocated memory */
  free(interp->data);
  interp->data = NULL;
}

#endif // SWIFT_GEAR_INTERPOLATION_H
