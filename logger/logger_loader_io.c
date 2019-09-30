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

/* SWIFT headers */
#include "io_properties.h"
#include "threadpool.h"

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
void *logger_loader_io_mmap_file(const char *filename, size_t *file_size,
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


/**
 * @brief Mapper function to copy a #logger_particle fields into a buffer.
 */
void logger_io_copy_mapper(void* restrict temp, int N, void* restrict extra_data) {

  const struct io_props props = *((const struct io_props*)extra_data);
  const size_t typeSize = io_sizeof_type(props.type);
  const size_t copySize = typeSize * props.dimension;

  /* How far are we with this chunk? */
  char* restrict temp_c = (char*)temp;
  const ptrdiff_t delta = (temp_c - props.start_temp_c) / copySize;

  for (int k = 0; k < N; k++) {
    memcpy(&temp_c[k * copySize], props.field + (delta + k) * props.partSize,
           copySize);
  }
}

/**
 * @brief Copy the particle data into a temporary buffer ready for i/o.
 *
 * @param temp The buffer to be filled. Must be allocated and aligned properly.
 * @param props The #io_props corresponding to the particle field we are
 * copying.
 * @param N The number of particles to copy
 */
void logger_io_copy_temp_buffer(void* temp, const struct threadpool *threadpool,
                                struct io_props props, size_t N) {

  const size_t typeSize = io_sizeof_type(props.type);
  const size_t copySize = typeSize * props.dimension;

  /* Copy particle data to temporary buffer */
  if (props.conversion == 0) { /* No conversion */

    /* Prepare some parameters */
    char* temp_c = (char*)temp;
    props.start_temp_c = temp_c;

    /* Copy the whole thing into a buffer */
    threadpool_map((struct threadpool*)threadpool, logger_io_copy_mapper,
                   temp_c, N, copySize, 0, (void*)&props);

  } else { /* Converting particle to data */
      error("Missing conversion function");
  }
}


  /**
 * @brief Writes a data array in given HDF5 group.
 *
 * @param grp The group in which to write.
 * @param props The #io_props of the field to read
 * @param N The number of particles to write.
 *
 * @todo A better version using HDF5 hyper-slabs to write the file directly from
 * the part array will be written once the structures have been stabilized.
 */
void logger_writeArray(hid_t grp, const struct threadpool *threadpool,
                       const struct io_props props, size_t N) {

  const size_t typeSize = io_sizeof_type(props.type);
  const size_t num_elements = N * props.dimension;

  /* Allocate temporary buffer */
  void* temp = NULL;
  if (swift_memalign("writebuff", (void**)&temp, IO_BUFFER_ALIGNMENT,
                     num_elements * typeSize) != 0)
    error("Unable to allocate temporary i/o buffer");

  /* Copy the particle data to the temporary buffer */
  logger_io_copy_temp_buffer(temp, threadpool, props, N);

  /* Create data space */
  const hid_t h_space = H5Screate(H5S_SIMPLE);
  if (h_space < 0)
    error("Error while creating data space for field '%s'.", props.name);

  int rank;
  hsize_t shape[2];
  hsize_t chunk_shape[2];

  if (props.dimension > 1) {
    rank = 2;
    shape[0] = N;
    shape[1] = props.dimension;
    chunk_shape[0] = 1 << 20; /* Just a guess...*/
    chunk_shape[1] = props.dimension;
  } else {
    rank = 1;
    shape[0] = N;
    shape[1] = 0;
    chunk_shape[0] = 1 << 20; /* Just a guess...*/
    chunk_shape[1] = 0;
  }

  /* Make sure the chunks are not larger than the dataset */
  if (chunk_shape[0] > N) chunk_shape[0] = N;

  /* Change shape of data space */
  hid_t h_err = H5Sset_extent_simple(h_space, rank, shape, shape);
  if (h_err < 0)
    error("Error while changing data space shape for field '%s'.", props.name);

  /* Dataset properties */
  const hid_t h_prop = H5Pcreate(H5P_DATASET_CREATE);

  /* Set chunk size */
  h_err = H5Pset_chunk(h_prop, rank, chunk_shape);
  if (h_err < 0)
    error("Error while setting chunk size (%llu, %llu) for field '%s'.",
          chunk_shape[0], chunk_shape[1], props.name);

  /* Impose check-sum to verify data corruption */
  h_err = H5Pset_fletcher32(h_prop);
  if (h_err < 0)
    error("Error while setting checksum options for field '%s'.", props.name);

  /* Create dataset */
  const hid_t h_data = H5Dcreate(grp, props.name, io_hdf5_type(props.type),
                                 h_space, H5P_DEFAULT, h_prop, H5P_DEFAULT);
  if (h_data < 0) error("Error while creating dataspace '%s'.", props.name);

  /* Write temporary buffer to HDF5 dataspace */
  h_err = H5Dwrite(h_data, io_hdf5_type(props.type), h_space, H5S_ALL,
                   H5P_DEFAULT, temp);
  if (h_err < 0) error("Error while writing data array '%s'.", props.name);

  /* Write the full description */
  io_write_attribute_s(h_data, "Description", props.description);

  /* Free and close everything */
  swift_free("writebuff", temp);
  H5Pclose(h_prop);
  H5Dclose(h_data);
  H5Sclose(h_space);
}
