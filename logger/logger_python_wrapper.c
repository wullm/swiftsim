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
#include "logger_python_tools.h"
#include "logger_reader.h"
#include "logger_time.h"

#ifdef HAVE_PYTHON
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct {
  PyObject_HEAD;

  /* Add logger stuff here */
  struct logger_reader reader;

  /* Is the logger ready? */
  int ready;

  /* basename of the logfile */
  char *basename;

  /* Verbose level */
  int verbose;

} PyObjectReader;

/**
 * @brief Deallocate the memory.
 */
static void Reader_dealloc(PyObjectReader *self) {
  if (self->basename != NULL) {
    free(self->basename);
    self->basename = NULL;
  }
  Py_TYPE(self)->tp_free((PyObject *)self);
}

/**
 * @brief Allocate the memory.
 */
static PyObject *Reader_new(PyTypeObject *type, PyObject *args,
                            PyObject *kwds) {
  PyObjectReader *self;

  self = (PyObjectReader *)type->tp_alloc(type, 0);
  self->ready = 0;
  self->basename = NULL;
  return (PyObject *)self;
}
/**
 * @brief Initializer of the Reader.
 *
 * @param basename The basename of the logger file.
 * @param verbose The verbosity level.
 *
 * @return The Reader object.
 */
static int Reader_init(PyObjectReader *self, PyObject *args, PyObject *kwds) {

  /* input variables. */
  char *basename = NULL;
  int verbose = 0;

  /* List of keyword arguments. */
  static char *kwlist[] = {"basename", "verbose", NULL};

  /* parse the arguments. */
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "s|i", kwlist, &basename,
                                   &verbose))
    return -1;

  /* Copy the arguments */
  self->verbose = verbose;

  size_t n_plus_null = strlen(basename) + 1;
  self->basename = (char *)malloc(n_plus_null * sizeof(char));
  if (self->basename == NULL) {
    error_python("Failed to allocate memory");
  }
  strcpy(self->basename, basename);

  return 0;
}

/**
 * @brief Read the minimal and maximal time.
 *
 * <b>returns</b> tuple containing min and max time.
 */
static PyObject *getTimeLimits(PyObject *self, PyObject *Py_UNUSED(ignored)) {
  if (!((PyObjectReader *)self)->ready) {
    error_python(
        "The logger is not ready yet."
        "Did you forget to open it with \"with\"?");
  }

  /* initialize the reader. */
  struct logger_reader *reader = &((PyObjectReader *)self)->reader;

  if (reader->verbose > 1) message("Reading time limits.");

  /* Get the time limits */
  double time_min = logger_reader_get_time_begin(reader);
  double time_max = logger_reader_get_time_end(reader);

  /* Create the output */
  PyObject *out = PyTuple_New(2);
  PyTuple_SetItem(out, 0, PyFloat_FromDouble(time_min));
  PyTuple_SetItem(out, 1, PyFloat_FromDouble(time_max));

  return (PyObject *)out;
}

/**
 * @brief Create a list of numpy array containing the fields.
 *
 * @param output A list of array of fields.
 * @param field_indices The indices (header ordering) of the requested fields.
 * @param n_fields The number of fields requested.
 * @param n_part The number of particles of each type.
 * @param n_tot The total number of particles.
 *
 * @return The python list of numpy array.
 */
__attribute__((always_inline)) INLINE static PyObject *
logger_loader_create_output(void **output, const int *field_indices,
                            const int n_fields, const uint64_t *n_part,
                            uint64_t n_tot) {

  struct logger_python_field python_fields[100];

  /* Create the python list */
  PyObject *list = PyList_New(n_fields);
  struct logger_python_field *current_field;

  /* Get the hydro fields */
  hydro_logger_generate_python(python_fields);
  int total_number_fields = hydro_logger_field_count;
  /* Get the gravity fields */
  gravity_logger_generate_python(python_fields + total_number_fields);
  total_number_fields += gravity_logger_field_count;
  /* Get the stars fields */
  stars_logger_generate_python(python_fields + total_number_fields);
  total_number_fields += stars_logger_field_count;

  /* Get all the requested fields */
  for (int i = 0; i < n_fields; i++) {
    /* Reset the variables. */
    current_field = NULL;
    total_number_fields = 0;

    /* Find in the hydro the field. */
    for (int local = 0; local < hydro_logger_field_count; local++) {
      const int global = hydro_logger_local_to_global[local];
      if (field_indices[i] == global) {
        current_field = &python_fields[local];
      }
    }
    total_number_fields += hydro_logger_field_count;

    /* Find in the gravity the field. */
    for (int local = 0; local < gravity_logger_field_count; local++) {
      const int global = gravity_logger_local_to_global[local];
      const int local_shifted = local + total_number_fields;
      if (field_indices[i] == global) {
        /* Check if we have the same fields for gravity + hydro */
        if (current_field != NULL) {
          if (current_field->dimension !=
                  python_fields[local_shifted].dimension ||
              current_field->typenum != python_fields[local_shifted].typenum) {
            error_python(
                "The python definition of the field %s does not correspond "
                "between"
                " the modules.",
                gravity_logger_field_names[local]);
          }
        }
        current_field = &python_fields[local_shifted];
        break;
      }
    }
    total_number_fields += gravity_logger_field_count;

    /* Find in the stars the field. */
    for (int local = 0; local < stars_logger_field_count; local++) {
      const int global = stars_logger_local_to_global[local];
      const int local_shifted = local + total_number_fields;
      if (field_indices[i] == global) {
        /* Check if we have the same fields for gravity + hydro + stars. */
        if (current_field != NULL) {
          if (current_field->dimension !=
                  python_fields[local_shifted].dimension ||
              current_field->typenum != python_fields[local_shifted].typenum) {
            error_python(
                "The python definition of the field %s does not correspond "
                "between"
                " the modules.",
                stars_logger_field_names[local]);
          }
        }
        current_field = &python_fields[local_shifted];
        break;
      }
    }
    total_number_fields += stars_logger_field_count;

    /* Check if we got a field */
    if (current_field == NULL) {
      error_python("Failed to find the required field");
    }
    PyObject *array = NULL;
    if (current_field->dimension > 1) {
      npy_intp dims[2] = {n_tot, current_field->dimension};
      array =
          PyArray_SimpleNewFromData(2, dims, current_field->typenum, output[i]);
    } else {
      npy_intp dims = n_tot;
      array = PyArray_SimpleNewFromData(1, &dims, current_field->typenum,
                                        output[i]);
    }

    PyList_SetItem(list, i, array);
  }

  return list;
}

static PyObject *pyEnter(__attribute__((unused)) PyObject *self,
                         PyObject *args) {

  PyObjectReader *self_reader = (PyObjectReader *)self;
  logger_reader_init(&self_reader->reader, self_reader->basename,
                     self_reader->verbose);
  self_reader->ready = 1;

  Py_INCREF(self);
  return self;
}

static PyObject *pyExit(__attribute__((unused)) PyObject *self,
                        PyObject *args) {
  PyObjectReader *self_reader = (PyObjectReader *)self;
  if (!self_reader->ready) {
    error_python("It seems that the reader was not initialized");
  }
  logger_reader_free(&self_reader->reader);
  self_reader->ready = 0;

  Py_RETURN_NONE;
}
/**
 * @brief Read some fields at a given time.
 *
 * @param basename The basename of the logger files.
 * @param fields Python list containing the name of the fields (e.g.
 * Coordinates).
 * @param time The time of the fields.
 * @param verbose (Optional) The verbose level of the reader.
 *
 * @return List of numpy array containing the fields requested (in the same
 * order).
 */
static PyObject *pyGetParticleData(__attribute__((unused)) PyObject *self,
                                   PyObject *args) {
  PyObjectReader *self_reader = (PyObjectReader *)self;
  if (!self_reader->ready) {
    error_python(
        "The logger is not ready yet."
        "Did you forget to open it with \"with\"?");
  }

  /* input variables. */
  PyObject *fields = NULL;
  double time = 0;
  /* parse the arguments. */
  if (!PyArg_ParseTuple(args, "Od", &fields, &time)) return NULL;

  /* Check the inputs. */
  if (!PyList_Check(fields)) {
    error_python("Expecting a list of fields");
  }

  /* initialize the reader. */
  struct logger_reader *reader = &self_reader->reader;
  const struct header *h = &reader->log.header;

  /* Get the fields indexes from the header. */
  const int n_fields = PyList_Size(fields);
  int *field_indices = (int *)malloc(n_fields * sizeof(int));
  for (int i = 0; i < n_fields; i++) {
    field_indices[i] = -1;

    /* Get the an item in the list. */
    PyObject *field = PyList_GetItem(fields, i);
    if (!PyUnicode_Check(field)) {
      error_python("Expecting list of string for the fields");
    }

    /* Convert into C string. */
    Py_ssize_t size = 0;
    const char *field_name = PyUnicode_AsUTF8AndSize(field, &size);

    /* Find the equivalent field inside the header. */
    for (int j = 0; j < h->masks_count; j++) {
      if (strcmp(field_name, h->masks[j].name) == 0) {
        field_indices[i] = j;
        break;
      }
    }

    /* Check if we found a field. */
    if (field_indices[i] == -1) {
      error_python("Failed to find the field %s", field_name);
    }
  }

  /* Set the time. */
  logger_reader_set_time(reader, time);

  /* Get the number of particles. */
  int n_type = 0;
  const uint64_t *n_part = logger_reader_get_number_particles(reader, &n_type);
  uint64_t n_tot = 0;
  for (int i = 0; i < n_type; i++) {
    n_tot += n_part[i];
  }

  /* Allocate the output memory. */
  void **output = malloc(n_fields * sizeof(void *));
  for (int i = 0; i < n_fields; i++) {
    output[i] = malloc(n_tot * h->masks[field_indices[i]].size);
  }

  /* Read the particles. */
  logger_reader_read_all_particles(reader, time, logger_reader_lin,
                                   field_indices, n_fields, output, n_part);

  /* Create the python output. */
  PyObject *array = logger_loader_create_output(output, field_indices, n_fields,
                                                n_part, n_tot);

  /* Free the reader. */
  free(field_indices);

  return array;
}

/* definition of the method table. */

static PyMethodDef libloggerMethods[] = {
    {NULL, NULL, 0, NULL} /* Sentinel */
};

/* Definition of the Reader methods */
static PyMethodDef libloggerReaderMethods[] = {
    {"get_time_limits", getTimeLimits, METH_NOARGS,
     "Read the time limits of the simulation.\n\n"
     "Returns\n"
     "-------\n\n"
     "times: tuple\n"
     "  time min, time max\n"},
    {"get_particle_data", pyGetParticleData, METH_VARARGS,
     "Read some fields from the logfile at a given time.\n\n"
     "Parameters\n"
     "----------\n\n"
     "fields: list\n"
     "  The list of fields (e.g. 'Coordinates', 'Entropies', ...)\n\n"
     "time: float\n"
     "  The time at which the fields must be read.\n\n"
     "-------\n\n"
     "list_of_fields: list\n"
     "  Each element is a numpy array containing the corresponding field.\n"},
    {"__enter__", pyEnter, METH_VARARGS, ""},
    {"__exit__", pyExit, METH_VARARGS, ""},
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

static PyTypeObject PyObjectReader_Type = {
    PyVarObject_HEAD_INIT(NULL, 0).tp_name = "liblogger.Reader",
    .tp_basicsize = sizeof(PyObjectReader),
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_doc =
        "Deals with a logger file and provides the data through"
        " its different methods. The reader is really open with the method "
        "__enter__"
        "(called by `with ...``). Outside the area delimited by with, "
        "the logger is not accessible.\n\n"
        "Parameters\n"
        "----------\n\n"
        "basename: str\n"
        "    Basename of the logfile.\n"
        "verbose: int (optional)\n"
        "    Verbose level\n\n"
        "Methods\n"
        "-------\n\n"
        "get_time_limits\n"
        "    Provides the initial and final time of the simulation\n"
        "get_particle_data\n"
        "    Reads some particle's data from the logfile\n\n"
        "Examples\n"
        "--------\n\n"
        ">>> import liblogger as logger\n"
        ">>> with logger.Reader(\"index_0000\") as reader:\n"
        ">>>    t0, t1 = reader.get_time_limits()\n"
        ">>>    pos, ent = reader.get_particle_data([\"Coordinates\", "
        "\"Entropies\"]"
        ", 0.5 * (t0 + t1))\n",
    .tp_init = (initproc)Reader_init,
    .tp_dealloc = (destructor)Reader_dealloc,
    .tp_new = Reader_new,
    .tp_methods = libloggerReaderMethods,
};

PyMODINIT_FUNC PyInit_liblogger(void) {

  /* Create the module. */
  PyObject *m;
  m = PyModule_Create(&libloggermodule);
  if (m == NULL) return NULL;

  /* Deal with SWIFT clock */
  clocks_set_cpufreq(0);

  /* Prepare the classes */
  if (PyType_Ready(&PyObjectReader_Type) < 0) {
    return NULL;
  }

  Py_INCREF(&PyObjectReader_Type);
  PyModule_AddObject(m, "Reader", (PyObject *)&PyObjectReader_Type);

  /* Import numpy. */
  import_array();

  return m;
}

#endif  // HAVE_PYTHON
