/**
 * @field_io.h
 * @author  Willem Elbers <whe@willemelbers.com>
 * @version 1.0
 * @license MIT (https://opensource.org/licenses/MIT)
 *
 * @section DESCRIPTION
 *
 * Loads and saves floating point grid values from binary files.
 */

#ifndef FIELD_IO_H
#define FIELD_IO_H

#include <fstream>
#include <iostream>
#include <vector>

// Array index (this is the row major format)
inline int box_idx(int width, int x, int y, int z) {
  return z + width * (y + width * x);
}

// Array index, but wraps around (i.e. the topology is T^3)
inline int box_wrap_idx(int width, int x, int y, int z) {
  int x2 = (x >= 0) ? x % width : x + ceil(-1. * x / width) * width;
  int y2 = (y >= 0) ? y % width : y + ceil(-1. * y / width) * width;
  int z2 = (z >= 0) ? z % width : z + ceil(-1. * z / width) * width;
  return box_idx(width, x2, y2, z2);
}

// Row major index for a half-complex array with N*N*(N/2+1) complex entries
inline int half_box_idx(int width, int x, int y, int z) {
  return z + (width / 2 + 1) * (y + width * x);
}

// Same as half_box_idx but for a T^3 topology
inline int half_box_wrap_idx(int width, int x, int y, int z) {
  int x2 = (x >= 0) ? x % width : x + ceil(-1. * x / width) * width;
  int y2 = (y >= 0) ? y % width : y + ceil(-1. * y / width) * width;
  int z2 = (z >= 0) ? z % (width / 2 + 1)
                    : z + ceil(-1. * z / (width / 2 + 1)) * (width / 2 + 1);
  return half_box_idx(width, x2, y2, z2);
}

// Read a vector of floats from a binary file
inline std::vector<float> read_floats(std::string fname) {
  // Open the file (in, binary, open at end)
  std::fstream file(
      fname, std::fstream::in | std::fstream::binary | std::fstream::ate);

  if (file.is_open()) {
    int size = file.tellg();         // current pos at end
    int float_size = sizeof(float);  // bytes
    int number_of_floats = size / float_size;
    int width = cbrt(number_of_floats);  // cube root

    if (width < 1 || width * width * width != number_of_floats) {
      std::cerr << "File is not a data cube of floats." << std::endl;
      exit(0);
    }

    std::vector<float> the_floats(number_of_floats);

    // Go to the beginning
    file.seekg(0, std::ios::beg);
    float f;
    int i = 0;
    while (file.read(reinterpret_cast<char*>(&f), float_size)) {
      the_floats[i] = f;
      i++;
    }

    file.close();

    return the_floats;
  } else {
    std::cerr << "File not found." << std::endl;
    exit(0);
  }
}

inline void write_floats(std::string fname, std::vector<float> floats) {
  // Open the file (out, binary)
  std::fstream file(fname, std::fstream::out | std::fstream::binary);

  if (file.is_open()) {
    int float_size = sizeof(float);  // bytes

    for (auto f : floats) {
      file.write(reinterpret_cast<const char*>(&f), float_size);
    }

    file.close();
  } else {
    std::cerr << "Could not write file." << std::endl;
    exit(0);
  }
}

inline void write_array_to_disk(std::string fname, double* x_box, int N) {
  int width = N;

  // Write the box to the disk so we can plot it with python
  std::vector<float> out_box(N * N * N);
  for (int x = 0; x < N; x++) {
    for (int y = 0; y < N; y++) {
      for (int z = 0; z < N; z++) {
        out_box[box_idx(width, x, y, z)] = x_box[box_idx(width, x, y, z)];
      }
    }
  }
  write_floats(fname, out_box);
}

// export_what can be 0 (real), 1 (imag), or magnitude (2)
#define REAL_PART 0
#define IMAG_PART 1
#define MAGNITUDE 2
inline void write_complex_array_to_disk(std::string fname, fftw_complex* x_box,
                                        int N, int export_what = MAGNITUDE) {
  if (!(export_what == REAL_PART || export_what == IMAG_PART ||
        export_what == MAGNITUDE)) {
    throw std::invalid_argument(
        "You must export either the real part (0) or imaginary part (1).");
  }
  int width = N;

  // Write the box to the disk so we can plot it with python
  std::vector<float> out_box(N * N * N);
  for (int x = 0; x < N; x++) {
    for (int y = 0; y < N; y++) {
      for (int z = 0; z < N; z++) {
        if (export_what == REAL_PART) {
          out_box[box_idx(width, x, y, z)] = x_box[box_idx(width, x, y, z)][0];
        } else if (export_what == IMAG_PART) {
          out_box[box_idx(width, x, y, z)] = x_box[box_idx(width, x, y, z)][1];
        } else if (export_what == MAGNITUDE) {
          out_box[box_idx(width, x, y, z)] =
              sqrt(pow(x_box[box_idx(width, x, y, z)][0], 2) +
                   pow(x_box[box_idx(width, x, y, z)][1], 2));
        }
      }
    }
  }
  write_floats(fname, out_box);
}

#endif
