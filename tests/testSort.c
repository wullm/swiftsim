/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2015 Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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

#include <fenv.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "swift.h"

/**
 * @brief Constructs a cell and all of its particle in a valid state prior to
 * a DOPAIR or DOSELF calcuation.
 *
 * @param n The cube root of the number of particles.
 * @param offset The position of the cell offset from (0,0,0).
 * @param size The cell size.
 * @param h The smoothing length of the particles in units of the inter-particle
 *separation.
 * @param density The density of the fluid.
 * @param partId The running counter of IDs.
 * @param pert The perturbation to apply to the particles in the cell in units
 *of the inter-particle separation.
 */
struct cell *make_cell(size_t count, double *offset, double size, double h,
                       double density, long long *partId, double pert) {
  const double volume = size * size * size;
  struct cell *cell = malloc(sizeof(struct cell));
  bzero(cell, sizeof(struct cell));

  if (posix_memalign((void **)&cell->parts, part_align,
                     count * sizeof(struct part)) != 0) {
    error("couldn't allocate particles, no. of particles: %d", (int)count);
  }
  bzero(cell->parts, count * sizeof(struct part));

  /* Cell properties */
  cell->split = 0;
  cell->h_max = h;
  cell->count = count;
  cell->dx_max = 0.;
  cell->h[0] = size;
  cell->h[1] = size;
  cell->h[2] = size;
  cell->loc[0] = offset[0];
  cell->loc[1] = offset[1];
  cell->loc[2] = offset[2];

  cell->ti_end_min = 1;
  cell->ti_end_max = 1;

  cell->sorted = 0;
  cell->sort = NULL;
  cell->sortsize = 0;

  /* Construct the parts */
  generate_random_positions(cell);
  
  struct part *part = cell->parts;
  for (size_t i = 0; i < count; ++i) {
    part->id = ++(*partId);
    part->mass = density * volume / count;
    part->ti_begin = 0;
    part->ti_end = 1;
    ++part;
  }

  return cell;
}

void clean_up(struct cell *ci) {
  free(ci->parts);
  free(ci->sort);
  free(ci);
}

/**
 * @brief Initializes all particles field to be ready for a density calculation
 */
void zero_particle_fields(struct cell *c) {
  for (size_t pid = 0; pid < c->count; pid++) {
    c->parts[pid].rho = 0.f;
    c->parts[pid].rho_dh = 0.f;
    hydro_init_part(&c->parts[pid]);
  }
}

/**
 * @brief Dump all the particles to a file
 */
void dump_particle_fields(char *fileName, struct cell *cj) {

  FILE *file = fopen(fileName, "w");

  /* Write header */
  fprintf(file,
          "# %4s %10s %10s %10s %10s %10s\n",
          "ID", "pos_x", "pos_y", "pos_z", "sort_i", "sort_d");
  
  /* Print each sorted distance along each sort direction. */
  for (size_t sid = 0; sid < 13; sid++) {
    fprintf(file,"# SID: %4zu\n",sid);
    for (size_t pjd = 0; pjd < cj->count; pjd++) {
      fprintf(
          file,
          "%6llu %10f %10f %10f %10d %10f\n",
          cj->parts[pjd].id, cj->parts[pjd].x[0], cj->parts[pjd].x[1],
          cj->parts[pjd].x[2], cj->sort[sid*(cj->count + 1) + pjd].i, cj->sort[sid*(cj->count + 1) + pjd].d);
    }
  }
  fclose(file);
}

/* And go... */
int main(int argc, char *argv[]) {
  size_t runs = 0, particles = 0;
  double h = 1.2348, size = 1., rho = 1.;
  double perturbation = 0.;
  char outputFileNameExtension[200] = "";
  char outputFileName[200] = "";

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

  /* Get some randomness going */
  srand(0);

  char c;
  while ((c = getopt(argc, argv, "m:s:h:p:r:t:d:f:v:")) != -1) {
    switch (c) {
      case 'h':
        sscanf(optarg, "%lf", &h);
        break;
      case 's':
        sscanf(optarg, "%lf", &size);
        break;
      case 'p':
        sscanf(optarg, "%zu", &particles);
        break;
      case 'r':
        sscanf(optarg, "%zu", &runs);
        break;
      case 'd':
        sscanf(optarg, "%lf", &perturbation);
        break;
      case 'm':
        sscanf(optarg, "%lf", &rho);
        break;
      case 'f':
        strcpy(outputFileNameExtension, optarg);
        break;
      case '?':
        error("Unknown option.");
        break;
    }
  }

  if (h < 0 || particles == 0 || runs == 0) {
    printf(
        "\nUsage: %s -p NUMBER_OF_PARTICLES -r NUMBER_OF_RUNS [OPTIONS...]\n"
        "\nGenerates a cell pair, filled with particles on a Cartesian grid."
        "\nThese are then interacted using runner_dopair1_density."
        "\n\nOptions:"
        "\n-h DISTANCE=1.2348 - Smoothing length in units of <x>"
        "\n-m rho             - Physical density in the cell"
        "\n-s size            - Physical size of the cell"
        "\n-d pert            - Perturbation to apply to the particles [0,1["
        "\n-f fileName        - Part of the file name used to save the dumps\n",
        argv[0]);
    exit(1);
  }

  /* Help users... */
  message("Smoothing length: h = %f", h * size);
  message("Kernel:               %s", kernel_name);
  message("Neighbour target: N = %f",
          h * h * h * 4.0 * M_PI * kernel_gamma3 / 3.0);
  message("Density target: rho = %f", rho);
  printf("\n");

  /* Build the infrastructure */
  struct space space;
  space.periodic = 0;
  space.h_max = h;

  struct engine engine;
  engine.s = &space;
  engine.time = 0.1f;
  engine.ti_current = 1;

  struct runner runner;
  runner.e = &engine;

  /* Construct the cell */
  struct cell *cell;
  static long long partId = 0;
  double offset[3] = {0,0,0};
  
  cell = make_cell(particles, offset, size, h, rho,
                    &partId, perturbation);
  
  ticks time = 0;
  for (size_t i = 0; i < runs; ++i) {

    /* Randomise the particle positions */
    generate_random_positions(cell);

    /* Make sure the sort occurs. */
    cell->sorted = 0;
    
    const ticks tic = getticks();

    /* Run the sort */
    runner_do_sort(NULL, cell, 0x1FFF, 0);

    const ticks toc = getticks();
    time += toc - tic;
    
  }

  /* Dump contents of cell */
  sprintf(outputFileName, "swift_dosort_new%s.dat",
              outputFileNameExtension);
  dump_particle_fields(outputFileName,cell);
  
  /* Output timing */
  message("New SWIFT sort took       : %15lli ticks.", time / runs);

  /* Now perform the original version for accuracy tests */
  
  /* Re-seed RNG */
  srand(0);
  
  generate_random_positions(cell);

  time = 0;
  for (size_t i = 0; i < runs; ++i) {
    
    /* Randomise the particle positions */
    generate_random_positions(cell);

    /* Make sure the sort occurs. */
    cell->sorted = 0;

    const ticks tic = getticks();

    /* Run the sort */
    runner_do_sort(NULL, cell, 0x1FFF, 0);

    const ticks toc = getticks();
    time += toc - tic;

  }

  /* Dump contents of cell */
  sprintf(outputFileName, "swift_dosort_old%s.dat",
              outputFileNameExtension);
  dump_particle_fields(outputFileName,cell);

  /* Output timing */
  message("Original SWIFT sort took  : %15lli ticks.", time / runs);

  /* Clean things to make the sanitizer happy ... */
  clean_up(cell);

  return 0;
}
