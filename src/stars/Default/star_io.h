/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_DEFAULT_STAR_IO_H
#define SWIFT_DEFAULT_STAR_IO_H

#include "hydro_properties.h"
#include "io_properties.h"
#

/**
 * @brief return the name of the stellar model
 */
INLINE static char *star_implementation(void) {
  return "Default";
}

/**
 * @brief Specifies which s-particle fields to read from a dataset
 *
 * @param sparts The s-particle array.
 * @param list The list of i/o properties to read.
 * @param num_fields The number of i/o fields to read.
 */
INLINE static void star_read_particles(struct spart* sparts,
                                       struct io_props* list, int* num_fields) {

  /* Say how much we want to read */
  *num_fields = 4;

  /* List what we want to read */
  list[0] = io_make_input_field("Coordinates", DOUBLE, 3, COMPULSORY,
                                UNIT_CONV_LENGTH, sparts, x);
  list[1] = io_make_input_field("Velocities", FLOAT, 3, COMPULSORY,
                                UNIT_CONV_SPEED, sparts, v);
  list[2] = io_make_input_field("Masses", FLOAT, 1, COMPULSORY, UNIT_CONV_MASS,
                                sparts, mass);
  list[3] = io_make_input_field("ParticleIDs", LONGLONG, 1, COMPULSORY,
                                UNIT_CONV_NO_UNITS, sparts, id);
}

/**
 * @brief Specifies which s-particle fields to write to a dataset
 *
 * @param sparts The s-particle array.
 * @param list The list of i/o properties to write.
 * @param num_fields The number of i/o fields to write.
 */
INLINE static void star_write_particles(const struct spart* sparts,
                                        struct io_props* list,
                                        int* num_fields) {

  /* Say how much we want to read */
  *num_fields = 4;

  /* List what we want to read */
  list[0] = io_make_output_field("Coordinates", DOUBLE, 3, UNIT_CONV_LENGTH,
                                 sparts, x);
  list[1] =
      io_make_output_field("Velocities", FLOAT, 3, UNIT_CONV_SPEED, sparts, v);
  list[2] =
      io_make_output_field("Masses", FLOAT, 1, UNIT_CONV_MASS, sparts, mass);
  list[3] = io_make_output_field("ParticleIDs", LONGLONG, 1, UNIT_CONV_NO_UNITS,
                                 sparts, id);
}


/**
 * @brief Initialize the global properties of the stellar model.
 *
 * @param s The #star_props.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 * @param p The #hydro_props (used for default options).
 */
INLINE static void star_props_init(struct star_props *s,
                      const struct phys_const *phys_const,
                      const struct unit_system *us,
		      struct swift_params *params,
		      const struct hydro_props *p) {

  /* Kernel properties */
  s->eta_neighbours = parser_get_opt_param_float(params, "Stars:resolution_eta",
						 p->eta_neighbours);

  /* Tolerance for the smoothing length Newton-Raphson scheme */
  s->h_tolerance = parser_get_opt_param_float(params, "Stars:h_tolerance",
                                              p->h_tolerance);

  /* Get derived properties */
  s->target_neighbours = pow_dimension(s->eta_neighbours) * kernel_norm;
  const float delta_eta = s->eta_neighbours * (1.f + s->h_tolerance);
  s->delta_neighbours =
      (pow_dimension(delta_eta) - pow_dimension(s->eta_neighbours)) *
      kernel_norm;

  /* Maximal smoothing length */
  s->h_max = parser_get_opt_param_float(params, "Stars:h_max",
                                        p->h_max);

  /* Number of iterations to converge h */
  s->max_smoothing_iterations = parser_get_opt_param_int(
      params, "Stars:max_ghost_iterations", p->max_smoothing_iterations);

  /* Max volume change */
  const float max_volume_change = parser_get_opt_param_float(
      params, "Stars:max_volume_change", p->log_max_h_change);
  s->log_max_h_change = logf(powf(max_volume_change, hydro_dimension_inv));

}

/**
 * @brief Print the global properties of the stellar model.
 *
 * @param s The #stars_props.
 */
INLINE static void star_props_print(const struct star_props *s) {

  /* Now star */
  message("Stellar model: %s.", star_implementation());

  message("Stellar kernel: %s with eta=%f (%.2f neighbours).", kernel_name,
          s->eta_neighbours, s->target_neighbours);

  message("Stellar relative tolerance in h: %.5f (+/- %.4f neighbours).",
          s->h_tolerance, s->delta_neighbours);

  message(
      "Hydrodynamic integration: Max change of volume: %.2f "
      "(max|dlog(h)/dt|=%f).",
      pow_dimension(expf(s->log_max_h_change)), s->log_max_h_change);

  message("Maximal smoothing length allowed: %.4f", s->h_max);

  message("Maximal iterations in ghost task set to %d",
	  s->max_smoothing_iterations);

}

#if defined(HAVE_HDF5)
INLINE static void star_props_print_snapshot(hid_t h_grpsph, const struct star_props *s) {

  io_write_attribute_s(h_grpsph, "Scheme", star_implementation());
  io_write_attribute_s(h_grpsph, "Kernel function", kernel_name);
  io_write_attribute_f(h_grpsph, "Kernel target N_ngb", s->target_neighbours);
  io_write_attribute_f(h_grpsph, "Kernel delta N_ngb", s->delta_neighbours);
  io_write_attribute_f(h_grpsph, "Kernel eta", s->eta_neighbours);
  io_write_attribute_f(h_grpsph, "Smoothing length tolerance", s->h_tolerance);
  io_write_attribute_f(h_grpsph, "Maximal smoothing length", s->h_max);
  io_write_attribute_f(h_grpsph, "Volume log(max(delta h))",
                       s->log_max_h_change);
  io_write_attribute_f(h_grpsph, "Volume max change time-step",
                       pow_dimension(expf(s->log_max_h_change)));
  io_write_attribute_i(h_grpsph, "Max ghost iterations",
                       s->max_smoothing_iterations);
}
#endif


#endif /* SWIFT_DEFAULT_STAR_IO_H */
