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
#ifndef SWIFT_STAR_IO_H
#define SWIFT_STAR_IO_H

#include "../config.h"

#include "./const.h"

#ifdef STARS_GEAR
#include "./stars/GEAR/star_io.h"
#else
#include "./stars/Default/star_io.h"
#endif

/**
 * @brief Write a #star_props struct to the given FILE as a stream of bytes.
 *
 * @param s the struct
 * @param stream the file stream
 */
INLINE static void star_props_struct_dump(const struct star_props *s, FILE *stream) {
  restart_write_blocks((void *)s, sizeof(struct star_props), 1, stream,
                       "starprops", "star props");
}

/**
 * @brief Restore a #star_props struct from the given FILE as a stream of
 * bytes.
 *
 * @param s the struct
 * @param stream the file stream
 */
INLINE static void star_props_struct_restore(const struct star_props *s, FILE *stream) {
  restart_read_blocks((void *)s, sizeof(struct star_props), 1, stream, NULL,
                      "star props");
}

#endif /* SWIFT_STAR_IO_H */
