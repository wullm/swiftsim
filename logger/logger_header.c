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
#include "logger_header.h"

#include "logger_loader_io.h"
#include "logger_logfile.h"
#include "logger_reader.h"
#include "logger_tools.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Name of each offset direction. */
const char *logger_offset_name[logger_offset_count] = {
    "Forward",
    "Backward",
    "Corrupted",
};

/**
 * @brief Print the properties of the header to stdout.
 *
 * @param h The #header.
 */
void header_print(const struct header *h) {
#ifdef SWIFT_DEBUG_CHECKS
  message("Debug checks enabled.");
#endif
  message("First Offset:     %lu.", h->offset_first_record);
  message("Offset direction: %s.", logger_offset_name[h->offset_direction]);
  message("Number masks:     %i.", h->number_mask);

  for (size_t i = 0; i < h->number_mask; i++) {
    message("  Mask:  %s.", h->masks[i].name);
    message("  Value: %u.", h->masks[i].mask);
    message("  Size:  %i.", h->masks[i].size);
    message("");
  }
};

/**
 * @brief free the allocated memory.
 *
 * @param h The #header.
 */
void header_free(struct header *h) { free(h->masks); };

/**
 * @brief Check if a field is present in the header.
 *
 * @param h The #header.
 * @param field name of the requested field.
 * @return Index of the field (-1 if not found).
 */
int header_get_field_index(const struct header *h, const char *field) {
  for (size_t i = 0; i < h->number_mask; i++) {
    if (strcmp(h->masks[i].name, field) == 0) {
      return i;
    }
  }

  return -1;
};

/**
 * @brief Update the offset direction in the structure and
 * write it to the logfile.
 *
 * @param h #header file structure.
 * @param new_value The new value to write.
 *
 */
void header_change_offset_direction(struct header *h,
                                    enum logger_offset_direction new_value) {
  h->offset_direction = new_value;
  /* Skip file format and version numbers. */
  size_t offset = LOGGER_VERSION_SIZE + 2 * sizeof(int);

  logger_loader_io_write_data(h->log->log.map + offset, sizeof(unsigned int),
                              &new_value);
}

/**
 * @brief read the logger header.
 *
 * @param h out: The #header.
 * @param log The #logger_logfile.
 */
void header_read(struct header *h, struct logger_logfile *log) {
  void *map = log->log.map;

  /* Set pointer to log. */
  h->log = log;

  /* read the file format. */
  char file_format[STRING_SIZE];
  map = logger_loader_io_read_data(map, LOGGER_VERSION_SIZE, &file_format);
  if (strcmp(file_format, "SWIFT_LOGGER"))
    error("Wrong file format (%s).", file_format);

  /* Read the major version number. */
  map = logger_loader_io_read_data(map, sizeof(int), &h->major_version);

  /* Read the minor version number. */
  map = logger_loader_io_read_data(map, sizeof(int), &h->minor_version);

  struct logger_reader *reader = log->reader;
  if (&reader->log != log) error("Wrong link to the reader.");

  if (reader->verbose > 0)
    message("File version %i.%i.", h->major_version, h->minor_version);

  /* Read the offset directions. */
  map = logger_loader_io_read_data(map, sizeof(int), &h->offset_direction);

  if (!header_is_forward(h) && !header_is_backward(h) &&
      !header_is_corrupted(h))
    error("Wrong offset value in the header (%i).", h->offset_direction);

  /* Read offset to first record. */
  map = logger_loader_io_read_data(map, LOGGER_OFFSET_SIZE,
                                   &h->offset_first_record);

  /* Read the size of the strings. */
  map =
      logger_loader_io_read_data(map, sizeof(unsigned int), &h->string_length);

  /* Check if value defined in this file is large enough. */
  if (STRING_SIZE < h->string_length) {
    error("Name too large in log file %i.", h->string_length);
  }

  /* Read the number of masks. */
  map = logger_loader_io_read_data(map, sizeof(unsigned int), &h->number_mask);

  /* Allocate the masks memory. */
  h->masks = malloc(sizeof(struct mask_data) * h->number_mask);

  /* Loop over all masks. */
  for (size_t i = 0; i < h->number_mask; i++) {
    /* Read the mask name. */
    map = logger_loader_io_read_data(map, h->string_length, h->masks[i].name);

    /* Set the mask value. */
    h->masks[i].mask = 1 << i;

    /* Read the mask data size. */
    map = logger_loader_io_read_data(map, sizeof(unsigned int),
                                     &h->masks[i].size);
  }

  /* Check the logfile header's size. */
  if (map != log->log.map + h->offset_first_record) {
    header_print(h);
    size_t offset = map - log->log.map;
    error("Wrong header size (in header %zi, current %zi).",
          h->offset_first_record, offset);
  }
};

/**
 * @brief Count number of bits in a given mask (without the record header).
 *
 * @param h #header file structure.
 * @param mask Mask to compute.
 *
 * @return number of bits in mask.
 */
size_t header_get_record_size_from_mask(const struct header *h,
                                        const size_t mask) {
  size_t count = 0;
  /* Loop over each masks. */
  for (size_t i = 0; i < h->number_mask; i++) {
    if (mask & h->masks[i].mask) {
      count += h->masks[i].size;
    }
  }
  return count;
}
