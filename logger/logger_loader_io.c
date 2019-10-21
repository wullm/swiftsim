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
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include "logger_header.h"
#include "logger_loader_io.h"
#include "logger_tools.h"

/**
 * @brief get the size of a file.
 *
 * @param fd file id.
 *
 * @return file size.
 */
size_t logger_loader_io_get_file_size(int fd) {
  struct stat s;
  int status = fstat(fd, &s);
  if (status != 0) error("Unable to get file size (%s).", strerror(errno));
  return s.st_size;
}

/**
 * @brief Map a file.
 *
 * #logger_loader_io_munmap_file should be called to unmap the file.
 *
 * @param filename file to read.
 * @param file_size (out) size of the file.
 * @param read_only Open the file in read only mode?
 *
 */
void *logger_loader_io_mmap_file(char *filename, size_t *file_size,
                                 int read_only) {
  /* open the file. */
  int fd;

  if (read_only)
    fd = open(filename, O_RDONLY);
  else
    fd = open(filename, O_RDWR);

  if (fd == -1)
    error("Unable to open file %s (%s).", filename, strerror(errno));

  /* get the file size. */
  *file_size = logger_loader_io_get_file_size(fd);

  /* map the memory. */
  int mode = PROT_READ;
  if (!read_only) mode |= PROT_WRITE;

  void *map = mmap(NULL, *file_size, mode, MAP_SHARED, fd, 0);
  if (map == MAP_FAILED)
    error("Failed to allocate map of size %zi bytes (%s).", *file_size,
          strerror(errno));

  /* Close the file. */
  close(fd);

  return map;
}

/**
 * @brief Unmap a file.
 *
 * @param map file mapping.
 * @param file_size The file size.
 *
 */
void logger_loader_io_munmap_file(void *map, size_t file_size) {
  /* unmap the file. */
  if (munmap(map, file_size) != 0) {
    error("Unable to unmap the file (%s).", strerror(errno));
  }
}
