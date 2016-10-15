/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <string.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "statistics.h"

/* Local headers. */
#include "atomic.h"
#include "cooling.h"
#include "engine.h"
#include "error.h"
#include "hydro.h"
#include "threadpool.h"

/**
 * @brief Information required to compute the statistics in the mapper
 */
struct index_data {
  const struct space *s;
  struct statistics *stats;
};

/**
 * @brief Adds the content of one #statistics aggregator to another one.
 *
 * Performs a += b;
 *
 * @param a The #statistics structure to update.
 * @param b The #statistics structure to add to a.
 */
void stats_add(struct statistics *a, const struct statistics *b) {

  a->E_kin += b->E_kin;
  a->E_int += b->E_int;
  a->E_pot += b->E_pot;
  a->E_rad += b->E_rad;
  a->entropy += b->entropy;
  a->mass += b->mass;
  a->mom[0] += b->mom[0];
  a->mom[1] += b->mom[1];
  a->mom[2] += b->mom[2];
  a->ang_mom[0] += b->ang_mom[0];
  a->ang_mom[1] += b->ang_mom[1];
  a->ang_mom[2] += b->ang_mom[2];
}

/**
 * @brief Initialises a statistics aggregator to a valid state.
 *
 * @param s The #statistics aggregator to initialise
 */
void stats_init(struct statistics *s) {

  /* Zero everything */
  bzero(s, sizeof(struct statistics));

  /* Set the lock */
  lock_init(&s->lock);
}

/**
 * @brief The #threadpool mapper function used to collect statistics
 *
 * @param map_data Pointer to the particles.
 * @param nr_parts The number of particles in this chunk
 * @param extra_data The #statistics aggregator.
 */
void stats_collect_part_mapper(void *map_data, int nr_parts, void *extra_data) {

  /* Unpack the data */
  struct index_data *data = (struct index_data *)extra_data;
  const struct space *s = data->s;
  struct part *restrict parts = (struct part *)map_data;
  struct xpart *restrict xparts = s->xparts + (ptrdiff_t)(parts - s->parts);
  const int ti_current = s->e->ti_current;
  const double timeBase = s->e->timeBase;
  struct statistics *const global_stats = data->stats;

  /* Local accumulator */
  struct statistics stats;
  bzero(&stats, sizeof(struct statistics));

  /* Loop over particles */
  for (int k = 0; k < nr_parts; k++) {

    /* Get the particle */
    struct part *restrict p = &parts[k];
    struct xpart *restrict xp = &xparts[k];

    const float dt = (ti_current - (p->ti_begin + p->ti_end) / 2) * timeBase;
    const double x[3] = {p->x[0], p->x[1], p->x[2]};
    const float v[3] = {xp->v_full[0] + p->a_hydro[0] * dt,
                        xp->v_full[1] + p->a_hydro[1] * dt,
                        xp->v_full[2] + p->a_hydro[2] * dt};

    const float m = hydro_get_mass(p);

    /* Collect mass */
    stats.mass += m;

    /* Collect momentum */
    stats.mom[0] += m * v[0];
    stats.mom[1] += m * v[1];
    stats.mom[2] += m * v[2];

    /* Collect angular momentum */
    stats.ang_mom[0] += m * (x[1] * v[2] - x[2] * v[1]);
    stats.ang_mom[1] += m * (x[2] * v[0] - x[0] * v[2]);
    stats.ang_mom[2] += m * (x[0] * v[1] - x[1] * v[0]);

    /* Collect energies. */
    stats.E_kin += 0.5f * m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    stats.E_pot += 0.;
    stats.E_int += m * hydro_get_internal_energy(p, dt);
    stats.E_rad += cooling_get_radiated_energy(xp);

    /* Collect entropy */
    stats.entropy += m * hydro_get_entropy(p, dt);
  }

  /* Now write back to memory */
  if (lock_lock(&global_stats->lock) == 0) stats_add(global_stats, &stats);
  if (lock_unlock(&global_stats->lock) != 0) error("Failed to unlock stats.");
}

/**
 * @brief Collect physical statistics over all particles in a #space.
 *
 * @param s The #space to collect from.
 * @param stats The #statistics aggregator to fill.
 */
void stats_collect(const struct space *s, struct statistics *stats) {

  const int chunk_size = 1000;

  /* Prepare the data */
  struct index_data extra_data;
  extra_data.s = s;
  extra_data.stats = stats;

  /* Run parallel collection of statistics */
  threadpool_map(&s->e->threadpool, stats_collect_part_mapper, s->parts,
                 s->nr_parts, sizeof(struct part), chunk_size, &extra_data);
}

/**
 * @brief Prints the content of a #statistics aggregator to a file
 *
 * @param file File to write to.
 * @param stats The #statistics object to write to the file
 * @param time The current physical time.
 */
void stats_print_to_file(FILE *file, const struct statistics *stats,
                         double time) {

  const double E_tot = stats->E_kin + stats->E_int + stats->E_pot;

  fprintf(file,
          " %14e %14e %14e %14e %14e %14e %14e %14e %14e %14e %14e %14e %14e "
          "%14e\n",
          time, stats->mass, E_tot, stats->E_kin, stats->E_int, stats->E_pot,
          stats->E_rad, stats->entropy, stats->mom[0], stats->mom[1],
          stats->mom[2], stats->ang_mom[0], stats->ang_mom[1],
          stats->ang_mom[2]);
  fflush(file);
}

#ifdef WITH_MPI

MPI_Datatype statistics_mpi_type;
MPI_Op statistics_mpi_reduce_op;

/**
 * @brief MPI reduce operator for #statistics structures.
 */
void stats_add_MPI(void *in, void *inout, int *len, MPI_Datatype *datatype) {

  for (int i = 0; i < *len; ++i)
    stats_add(&((struct statistics *)inout)[0], &((struct statistics *)in)[i]);
}

/**
 * @brief Registers MPI #statistics type and reduction function.
 */
void stats_create_MPI_type() {

  /* This is not the recommended way of doing this.
     One should define the structure field by field
     But as long as we don't do serialization via MPI-IO
     we don't really care.
     Also we would have to modify this function everytime something
     is added to the statistics structure. */
  if (MPI_Type_contiguous(sizeof(struct statistics) / sizeof(unsigned char),
                          MPI_BYTE, &statistics_mpi_type) != MPI_SUCCESS ||
      MPI_Type_commit(&statistics_mpi_type) != MPI_SUCCESS) {
    error("Failed to create MPI type for statistics.");
  }

  /* Create the reduction operation */
  MPI_Op_create(stats_add_MPI, 1, &statistics_mpi_reduce_op);
}
#endif
