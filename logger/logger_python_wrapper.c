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
#include "logger_particle.h"
#include "logger_reader.h"
#include "logger_time.h"

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <errno.h>
#include <numpy/arrayobject.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * @brief load data from the offset without any interpolation
 *
 * <b>offset</b> PyArrayObject list of offset for each particle.
 *
 * <b>filename</b> string filename of the log file.
 *
 * <b>verbose</b> Verbose level.
 *
 * <b>returns</b> dictionnary containing the data read.
 */
static PyObject *loadFromIndex(__attribute__((unused)) PyObject *self,
                               PyObject *args) {

  /* input variables. */
  PyArrayObject *offset = NULL;
  char *filename = NULL;

  /* output variables. */
  PyArrayObject *pos = NULL;
  PyArrayObject *vel = NULL;
  PyArrayObject *acc = NULL;
  PyArrayObject *entropy = NULL;
  PyArrayObject *h_sph = NULL;
  PyArrayObject *rho = NULL;
  PyArrayObject *mass = NULL;
  PyArrayObject *id = NULL;

  size_t time_offset;
  int verbose = 2;

  /* parse arguments. */
  if (!PyArg_ParseTuple(args, "OsL|i", &offset, &filename, &time_offset,
                        &verbose))
    return NULL;

  if (!PyArray_Check(offset)) {
    error("Offset is not a numpy array.");
  }
  if (PyArray_NDIM(offset) != 1) {
    error("Offset is not a 1 dimensional array.");
  }
  if (PyArray_TYPE(offset) != NPY_UINT64) {
    error("Offset does not contain unsigned int.");
  }

  /* initialize the reader. */
  struct logger_reader reader;
  logger_reader_init(&reader, filename, verbose);
  struct header *h = &reader.log.header;

  /* init array. */
  npy_intp dim[2];
  dim[0] = PyArray_DIMS(offset)[0];
  dim[1] = DIM;

  /* Get required time. */
  double time = time_array_get_time(&reader.log.times, time_offset);

  /* init output. */
  if (header_get_field_index(h, "positions") != -1) {
    pos = (PyArrayObject *)PyArray_SimpleNew(2, dim, NPY_DOUBLE);
  }

  if (header_get_field_index(h, "velocities") != -1) {
    vel = (PyArrayObject *)PyArray_SimpleNew(2, dim, NPY_FLOAT);
  }

  if (header_get_field_index(h, "accelerations") != -1) {
    acc = (PyArrayObject *)PyArray_SimpleNew(2, dim, NPY_FLOAT);
  }

  if (header_get_field_index(h, "entropy") != -1) {
    entropy =
        (PyArrayObject *)PyArray_SimpleNew(1, PyArray_DIMS(offset), NPY_FLOAT);
  }

  if (header_get_field_index(h, "smoothing length") != -1) {
    h_sph =
        (PyArrayObject *)PyArray_SimpleNew(1, PyArray_DIMS(offset), NPY_FLOAT);
  }

  if (header_get_field_index(h, "density") != -1) {
    rho =
        (PyArrayObject *)PyArray_SimpleNew(1, PyArray_DIMS(offset), NPY_FLOAT);
  }

  if (header_get_field_index(h, "consts") != -1) {
    mass =
        (PyArrayObject *)PyArray_SimpleNew(1, PyArray_DIMS(offset), NPY_FLOAT);
    id = (PyArrayObject *)PyArray_SimpleNew(1, PyArray_DIMS(offset), NPY_ULONG);
  }

  if (verbose > 1) message("Reading particles.");

  /* loop over all particles. */
  for (npy_intp i = 0; i < PyArray_DIMS(offset)[0]; i++) {
    struct logger_particle part;

    /* Get the offset. */
    size_t offset_particle = *(size_t *)PyArray_GETPTR1(offset, i);

    /* Read the particle. */
    logger_particle_read(&part, &reader, offset_particle, time,
                         logger_reader_lin);

    double *dtmp;
    float *ftmp;
    size_t *stmp;

    /* copy the data. */
    for (size_t k = 0; k < DIM; k++) {
      if (pos) {
        dtmp = PyArray_GETPTR2(pos, i, k);
        *dtmp = part.pos[k];
      }

      if (vel) {
        ftmp = PyArray_GETPTR2(vel, i, k);
        *ftmp = part.vel[k];
      }

      if (acc) {
        ftmp = PyArray_GETPTR2(acc, i, k);
        *ftmp = part.acc[k];
      }
    }

    if (entropy) {
      ftmp = PyArray_GETPTR1(entropy, i);
      *ftmp = part.entropy;
    }

    if (rho) {
      ftmp = PyArray_GETPTR1(rho, i);
      *ftmp = part.density;
    }

    if (h_sph) {
      ftmp = PyArray_GETPTR1(h_sph, i);
      *ftmp = part.h;
    }

    if (mass) {
      ftmp = PyArray_GETPTR1(mass, i);
      *ftmp = part.mass;
    }

    if (id) {
      stmp = PyArray_GETPTR1(id, i);
      *stmp = part.id;
    }
  }

  /* Free the memory. */
  logger_reader_free(&reader);

  /* construct return value. */
  PyObject *dict = PyDict_New();
  PyObject *key = PyUnicode_FromString("positions");
  PyDict_SetItem(dict, key, PyArray_Return(pos));

  if (vel) {
    key = PyUnicode_FromString("velocities");
    PyDict_SetItem(dict, key, PyArray_Return(vel));
  }

  if (acc) {
    key = PyUnicode_FromString("accelerations");
    PyDict_SetItem(dict, key, PyArray_Return(acc));
  }

  if (entropy) {
    key = PyUnicode_FromString("entropy");
    PyDict_SetItem(dict, key, PyArray_Return(entropy));
  }

  if (rho) {
    key = PyUnicode_FromString("rho");
    PyDict_SetItem(dict, key, PyArray_Return(rho));
  }

  if (h_sph) {
    key = PyUnicode_FromString("h_sph");
    PyDict_SetItem(dict, key, PyArray_Return(h_sph));
  }

  if (mass) {
    key = PyUnicode_FromString("mass");
    PyDict_SetItem(dict, key, PyArray_Return(mass));
  }

  if (id) {
    key = PyUnicode_FromString("id");
    PyDict_SetItem(dict, key, PyArray_Return(id));
  }

  return dict;
}

/**
 * @brief Reverse offset in log file
 *
 * <b>filename</b> string filename of the log file
 * <b>verbose</b> Verbose level
 */
static PyObject *pyReverseOffset(__attribute__((unused)) PyObject *self,
                                 PyObject *args) {
  /* input variables. */
  char *filename = NULL;

  int verbose = 0;

  /* parse the arguments. */
  if (!PyArg_ParseTuple(args, "s|i", &filename, &verbose)) return NULL;

  /* initialize the reader which reverse the offset if necessary. */
  struct logger_reader reader;
  logger_reader_init(&reader, filename, verbose);

  /* Free the reader. */
  logger_reader_free(&reader);

  return Py_BuildValue("");
}

/* definition of the method table. */

static PyMethodDef libloggerMethods[] = {
    {"loadFromIndex", loadFromIndex, METH_VARARGS,
     "Load snapshot directly from the offset in an index file."},
    {"reverseOffset", pyReverseOffset, METH_VARARGS,
     "Reverse the offset (from pointing backward to forward)."},

    {NULL, NULL, 0, NULL} /* Sentinel */
};

static struct PyModuleDef libloggermodule = {
    PyModuleDef_HEAD_INIT,
    "liblogger",
    "Module reading a SWIFTsim logger snapshot",
    -1,
    libloggerMethods,
    NULL, /* m_slots */
    NULL, /* m_traverse */
    NULL, /* m_clear */
    NULL  /* m_free */
};

PyMODINIT_FUNC PyInit_liblogger(void) {
  PyObject *m;
  m = PyModule_Create(&libloggermodule);
  if (m == NULL) return NULL;

  import_array();

  return m;
}
