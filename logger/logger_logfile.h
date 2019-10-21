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
/**
 * @file logger_logfile.h
 * @brief This file contains the high level function for the log.
 */
#ifndef LOGGER_LOGGER_LOGFILE_H
#define LOGGER_LOGGER_LOGFILE_H

#include "logger_header.h"
#include "logger_time.h"

struct logger_reader;

/**
 * @brief This structure deals with the log file.
 *
 * This structure is initialized by the #logger_reader
 * and deals with the log file.
 * It maps it, reverse the offsets (if required) and unmap it.
 *
 * The structure is initialized with #logger_logfile_init_from_file and
 * freed with #logger_logfile_free.
 */
struct logger_logfile {

  /* Information contained in the file header. */
  struct header header;

  /* The reader that is using this log file. */
  struct logger_reader *reader;

  /* Information about the time records. */
  struct time_array times;

  /* The log's variables. */
  struct {
    /* Mapped data. */
    void *map;

    /* File size. */
    size_t file_size;

  } log;
};

void logger_logfile_init_from_file(struct logger_logfile *log, char *filename,
                                   struct logger_reader *reader,
                                   int only_header);
void logger_logfile_reverse_offset(struct logger_logfile *log, char *filename);
void logger_logfile_free(struct logger_logfile *log);

#endif  // LOGGER_LOGGER_LOGFILE_H
