/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
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

/* This object's header. */
#include "engine.h"

/* Local headers. */
#include "active.h"
#include "timeline.h"

/**
 * @brief Data collected from the cells at the end of a time-step
 */
struct end_of_step_data {

  size_t updated, g_updated, s_updated, b_updated;
  size_t inhibited, g_inhibited, s_inhibited, b_inhibited;
  integertime_t ti_hydro_end_min, ti_hydro_end_max, ti_hydro_beg_max;
  integertime_t ti_gravity_end_min, ti_gravity_end_max, ti_gravity_beg_max;
  integertime_t ti_stars_end_min, ti_stars_end_max, ti_stars_beg_max;
  integertime_t ti_black_holes_end_min, ti_black_holes_end_max,
      ti_black_holes_beg_max;
  struct engine *e;
  struct star_formation_history sfh;
};

#ifdef WITH_ENGINEERING

void engine_recurse_fix_timesteps(struct cell *c, struct engine *e){
 
  c->hydro.ti_end_min = e->ti_end_min;
  c->hydro.ti_end_max = e->ti_end_min;
  if(c->split){
    for (int k = 0; k < 8; k++) {
      struct cell *cp = c->progeny[k];
      if(cp != NULL){
        engine_recurse_fix_timesteps(cp, e);
      }
    }

  }else if(c->nodeID == e->nodeID){
    for(int k = 0; k < c->hydro.count; k++){
      struct part *p = &c->hydro.parts[k];
      p->time_bin = e->min_active_bin;
    }
  }
}

void engine_fix_timestep_mapper(void *map_data, int num_elements,
                                        void *extra_data){
  struct engine *e = (struct engine*) extra_data;
  struct cell *current_cells = (struct cell*)map_data;
//  int *cell_indices = e->s->local_cells_top;
  for (int ind = 0; ind < num_elements; ind++) {
    struct cell *c = &current_cells[ind];
    engine_recurse_fix_timesteps(c, e);

  }
}

void engine_fix_timestep(struct engine *e){
  /* Fix the particle timesteps to the minimum computed. */
#ifdef WITH_MPI
  MPI_Allreduce(MPI_IN_PLACE, &e->ti_end_min, 1, MPI_LONG_LONG, MPI_MIN, MPI_COMM_WORLD);
#endif
  threadpool_map(&e->threadpool, engine_fix_timestep_mapper,
                 e->s->cells_top, e->s->nr_cells/*local_cells*/,
                 sizeof(struct cell), 0, e);

}
#endif



/**
 * @brief Recursive function gathering end-of-step data.
 *
 * We recurse until we encounter a timestep or time-step MPI recv task
 * as the values will have been set at that level. We then bring these
 * values upwards.
 *
 * @param c The #cell to recurse into.
 * @param e The #engine.
 */
void engine_collect_end_of_step_recurse_hydro(struct cell *c,
                                              const struct engine *e) {

  /* Skip super-cells (Their values are already set) */
  if (c->timestep != NULL) return;
#ifdef WITH_MPI
  if (cell_get_recv(c, task_subtype_tend_part) != NULL) return;
#endif /* WITH_MPI */

#ifdef SWIFT_DEBUG_CHECKS
    /* if (!c->split) error("Reached a leaf without finding a time-step task!
     * c->depth=%d c->maxdepth=%d c->count=%d c->node=%d", */
    /* 		       c->depth, c->maxdepth, c->hydro.count, c->nodeID); */
#endif

  /* Counters for the different quantities. */
  size_t updated = 0;
  integertime_t ti_hydro_end_min = max_nr_timesteps, ti_hydro_end_max = 0,
                ti_hydro_beg_max = 0;

  /* Local Star formation history properties */
  struct star_formation_history sfh_updated;

  /* Initialize the star formation structs */
  star_formation_logger_init(&sfh_updated);

  /* Collect the values from the progeny. */
  for (int k = 0; k < 8; k++) {
    struct cell *cp = c->progeny[k];
    if (cp != NULL && cp->hydro.count > 0) {

      /* Recurse */
      engine_collect_end_of_step_recurse_hydro(cp, e);

      /* And update */
      ti_hydro_end_min = min(ti_hydro_end_min, cp->hydro.ti_end_min);
      ti_hydro_end_max = max(ti_hydro_end_max, cp->hydro.ti_end_max);
      ti_hydro_beg_max = max(ti_hydro_beg_max, cp->hydro.ti_beg_max);

      updated += cp->hydro.updated;

      /* Check if the cell is inactive and in that case reorder the SFH */
      if (!cell_is_starting_hydro(cp, e)) {
        star_formation_logger_log_inactive_cell(&cp->stars.sfh);
      }

      /* Add the star formation history in this cell to sfh_updated */
      star_formation_logger_add(&sfh_updated, &cp->stars.sfh);

      /* Collected, so clear for next time. */
      cp->hydro.updated = 0;
    }
  }

  /* Store the collected values in the cell. */
  c->hydro.ti_end_min = ti_hydro_end_min;
  c->hydro.ti_end_max = ti_hydro_end_max;
  c->hydro.ti_beg_max = ti_hydro_beg_max;
  c->hydro.updated = updated;
  // c->hydro.inhibited = inhibited;

  /* Store the star formation history in the parent cell */
  star_formation_logger_add(&c->stars.sfh, &sfh_updated);
}

/**
 * @brief Recursive function gathering end-of-step data.
 *
 * We recurse until we encounter a timestep or time-step MPI recv task
 * as the values will have been set at that level. We then bring these
 * values upwards.
 *
 * @param c The #cell to recurse into.
 * @param e The #engine.
 */
void engine_collect_end_of_step_recurse_grav(struct cell *c,
                                             const struct engine *e) {

  /* Skip super-cells (Their values are already set) */
  if (c->timestep != NULL) return;
#ifdef WITH_MPI
  if (cell_get_recv(c, task_subtype_tend_gpart) != NULL) return;
#endif /* WITH_MPI */

#ifdef SWIFT_DEBUG_CHECKS
    //  if (!c->split) error("Reached a leaf without finding a time-step
    //  task!");
#endif

  /* Counters for the different quantities. */
  size_t updated = 0;
  integertime_t ti_grav_end_min = max_nr_timesteps, ti_grav_end_max = 0,
                ti_grav_beg_max = 0;

  /* Collect the values from the progeny. */
  for (int k = 0; k < 8; k++) {
    struct cell *cp = c->progeny[k];
    if (cp != NULL && cp->grav.count > 0) {

      /* Recurse */
      engine_collect_end_of_step_recurse_grav(cp, e);

      /* And update */
      ti_grav_end_min = min(ti_grav_end_min, cp->grav.ti_end_min);
      ti_grav_end_max = max(ti_grav_end_max, cp->grav.ti_end_max);
      ti_grav_beg_max = max(ti_grav_beg_max, cp->grav.ti_beg_max);

      updated += cp->grav.updated;

      /* Collected, so clear for next time. */
      cp->grav.updated = 0;
    }
  }

  /* Store the collected values in the cell. */
  c->grav.ti_end_min = ti_grav_end_min;
  c->grav.ti_end_max = ti_grav_end_max;
  c->grav.ti_beg_max = ti_grav_beg_max;
  c->grav.updated = updated;
}

/**
 * @brief Recursive function gathering end-of-step data.
 *
 * We recurse until we encounter a timestep or time-step MPI recv task
 * as the values will have been set at that level. We then bring these
 * values upwards.
 *
 * @param c The #cell to recurse into.
 * @param e The #engine.
 */
void engine_collect_end_of_step_recurse_stars(struct cell *c,
                                              const struct engine *e) {

  /* Skip super-cells (Their values are already set) */
  if (c->timestep != NULL) return;
#ifdef WITH_MPI
  if (cell_get_recv(c, task_subtype_tend_spart) != NULL) return;
#endif /* WITH_MPI */

#ifdef SWIFT_DEBUG_CHECKS
    // if (!c->split) error("Reached a leaf without finding a time-step task!");
#endif

  /* Counters for the different quantities. */
  size_t updated = 0;
  integertime_t ti_stars_end_min = max_nr_timesteps, ti_stars_end_max = 0,
                ti_stars_beg_max = 0;

  /* Collect the values from the progeny. */
  for (int k = 0; k < 8; k++) {
    struct cell *cp = c->progeny[k];
    if (cp != NULL && cp->stars.count > 0) {

      /* Recurse */
      engine_collect_end_of_step_recurse_stars(cp, e);

      /* And update */
      ti_stars_end_min = min(ti_stars_end_min, cp->stars.ti_end_min);
      ti_stars_end_max = max(ti_stars_end_max, cp->stars.ti_end_max);
      ti_stars_beg_max = max(ti_stars_beg_max, cp->stars.ti_beg_max);

      updated += cp->stars.updated;

      /* Collected, so clear for next time. */
      cp->stars.updated = 0;
    }
  }

  /* Store the collected values in the cell. */
  c->stars.ti_end_min = ti_stars_end_min;
  c->stars.ti_end_max = ti_stars_end_max;
  c->stars.ti_beg_max = ti_stars_beg_max;
  c->stars.updated = updated;
}

/**
 * @brief Recursive function gathering end-of-step data.
 *
 * We recurse until we encounter a timestep or time-step MPI recv task
 * as the values will have been set at that level. We then bring these
 * values upwards.
 *
 * @param c The #cell to recurse into.
 * @param e The #engine.
 */
void engine_collect_end_of_step_recurse_black_holes(struct cell *c,
                                                    const struct engine *e) {

  /* Skip super-cells (Their values are already set) */
  if (c->timestep != NULL) return;
#ifdef WITH_MPI
  if (cell_get_recv(c, task_subtype_tend_bpart) != NULL) return;
#endif /* WITH_MPI */

#ifdef SWIFT_DEBUG_CHECKS
    // if (!c->split) error("Reached a leaf without finding a time-step task!");
#endif

  /* Counters for the different quantities. */
  size_t updated = 0;
  integertime_t ti_black_holes_end_min = max_nr_timesteps,
                ti_black_holes_end_max = 0, ti_black_holes_beg_max = 0;

  /* Collect the values from the progeny. */
  for (int k = 0; k < 8; k++) {
    struct cell *cp = c->progeny[k];
    if (cp != NULL && cp->black_holes.count > 0) {

      /* Recurse */
      engine_collect_end_of_step_recurse_black_holes(cp, e);

      /* And update */
      ti_black_holes_end_min =
          min(ti_black_holes_end_min, cp->black_holes.ti_end_min);
      ti_black_holes_end_max =
          max(ti_black_holes_end_max, cp->black_holes.ti_end_max);
      ti_black_holes_beg_max =
          max(ti_black_holes_beg_max, cp->black_holes.ti_beg_max);

      updated += cp->black_holes.updated;

      /* Collected, so clear for next time. */
      cp->black_holes.updated = 0;
    }
  }

  /* Store the collected values in the cell. */
  c->black_holes.ti_end_min = ti_black_holes_end_min;
  c->black_holes.ti_end_max = ti_black_holes_end_max;
  c->black_holes.ti_beg_max = ti_black_holes_beg_max;
  c->black_holes.updated = updated;
}

/**
 * @brief Mapping function to collect the data from the end of the step
 *
 * This function will call a recursive function on all the top-level cells
 * to collect the information we are after.
 *
 * @param map_data The list of cells with tasks on this node.
 * @param num_elements The number of elements in the list this thread will work
 * on.
 * @param extra_data The #engine.
 */
void engine_collect_end_of_step_mapper(void *map_data, int num_elements,
                                       void *extra_data) {

  struct end_of_step_data *data = (struct end_of_step_data *)extra_data;
  const struct engine *e = data->e;
  const int with_hydro = (e->policy & engine_policy_hydro);
  const int with_self_grav = (e->policy & engine_policy_self_gravity);
  const int with_ext_grav = (e->policy & engine_policy_external_gravity);
  const int with_grav = (with_self_grav || with_ext_grav);
  const int with_stars = (e->policy & engine_policy_stars);
  const int with_black_holes = (e->policy & engine_policy_black_holes);
  struct space *s = e->s;
  int *local_cells = (int *)map_data;
  struct star_formation_history *sfh_top = &data->sfh;

  /* Local collectible */
  size_t updated = 0, g_updated = 0, s_updated = 0, b_updated = 0;
  integertime_t ti_hydro_end_min = max_nr_timesteps, ti_hydro_end_max = 0,
                ti_hydro_beg_max = 0;
  integertime_t ti_gravity_end_min = max_nr_timesteps, ti_gravity_end_max = 0,
                ti_gravity_beg_max = 0;
  integertime_t ti_stars_end_min = max_nr_timesteps, ti_stars_end_max = 0,
                ti_stars_beg_max = 0;
  integertime_t ti_black_holes_end_min = max_nr_timesteps,
                ti_black_holes_end_max = 0, ti_black_holes_beg_max = 0;

  /* Local Star formation history properties */
  struct star_formation_history sfh_updated;

  /* Initialize the star formation structs for this engine to zero */
  star_formation_logger_init(&sfh_updated);

  for (int ind = 0; ind < num_elements; ind++) {
    struct cell *c = &s->cells_top[local_cells[ind]];

    if (c->hydro.count > 0 || c->grav.count > 0 || c->stars.count > 0 ||
        c->black_holes.count > 0) {

      /* Make the top-cells recurse */
      if (with_hydro) {
        engine_collect_end_of_step_recurse_hydro(c, e);
      }
      if (with_grav) {
        engine_collect_end_of_step_recurse_grav(c, e);
      }
      if (with_stars) {
        engine_collect_end_of_step_recurse_stars(c, e);
      }
      if (with_black_holes) {
        engine_collect_end_of_step_recurse_black_holes(c, e);
      }

      /* And aggregate */
      if (c->hydro.ti_end_min > e->ti_current)
        ti_hydro_end_min = min(ti_hydro_end_min, c->hydro.ti_end_min);
      ti_hydro_end_max = max(ti_hydro_end_max, c->hydro.ti_end_max);
      ti_hydro_beg_max = max(ti_hydro_beg_max, c->hydro.ti_beg_max);

      if (c->grav.ti_end_min > e->ti_current)
        ti_gravity_end_min = min(ti_gravity_end_min, c->grav.ti_end_min);
      ti_gravity_end_max = max(ti_gravity_end_max, c->grav.ti_end_max);
      ti_gravity_beg_max = max(ti_gravity_beg_max, c->grav.ti_beg_max);

      if (c->stars.ti_end_min > e->ti_current)
        ti_stars_end_min = min(ti_stars_end_min, c->stars.ti_end_min);
      ti_stars_end_max = max(ti_stars_end_max, c->stars.ti_end_max);
      ti_stars_beg_max = max(ti_stars_beg_max, c->stars.ti_beg_max);

      if (c->black_holes.ti_end_min > e->ti_current)
        ti_black_holes_end_min =
            min(ti_black_holes_end_min, c->black_holes.ti_end_min);
      ti_black_holes_end_max =
          max(ti_black_holes_end_max, c->black_holes.ti_end_max);
      ti_black_holes_beg_max =
          max(ti_black_holes_beg_max, c->black_holes.ti_beg_max);

      updated += c->hydro.updated;
      g_updated += c->grav.updated;
      s_updated += c->stars.updated;
      b_updated += c->black_holes.updated;

      /* Check if the cell is inactive and in that case reorder the SFH */
      if (!cell_is_starting_hydro(c, e)) {
        star_formation_logger_log_inactive_cell(&c->stars.sfh);
      }

      /* Get the star formation history from the current cell and store it in
       * the star formation history struct */
      star_formation_logger_add(&sfh_updated, &c->stars.sfh);

      /* Collected, so clear for next time. */
      c->hydro.updated = 0;
      c->grav.updated = 0;
      c->stars.updated = 0;
      c->black_holes.updated = 0;
    }
  }

  /* Let's write back to the global data.
   * We use the space lock to garanty single access*/
  if (lock_lock(&s->lock) == 0) {
    data->updated += updated;
    data->g_updated += g_updated;
    data->s_updated += s_updated;
    data->b_updated += b_updated;

    /* Add the SFH information from this engine to the global data */
    star_formation_logger_add(sfh_top, &sfh_updated);

    if (ti_hydro_end_min > e->ti_current)
      data->ti_hydro_end_min = min(ti_hydro_end_min, data->ti_hydro_end_min);
    data->ti_hydro_end_max = max(ti_hydro_end_max, data->ti_hydro_end_max);
    data->ti_hydro_beg_max = max(ti_hydro_beg_max, data->ti_hydro_beg_max);

    if (ti_gravity_end_min > e->ti_current)
      data->ti_gravity_end_min =
          min(ti_gravity_end_min, data->ti_gravity_end_min);
    data->ti_gravity_end_max =
        max(ti_gravity_end_max, data->ti_gravity_end_max);
    data->ti_gravity_beg_max =
        max(ti_gravity_beg_max, data->ti_gravity_beg_max);

    if (ti_stars_end_min > e->ti_current)
      data->ti_stars_end_min = min(ti_stars_end_min, data->ti_stars_end_min);
    data->ti_stars_end_max = max(ti_stars_end_max, data->ti_stars_end_max);
    data->ti_stars_beg_max = max(ti_stars_beg_max, data->ti_stars_beg_max);

    if (ti_black_holes_end_min > e->ti_current)
      data->ti_black_holes_end_min =
          min(ti_black_holes_end_min, data->ti_black_holes_end_min);
    data->ti_black_holes_end_max =
        max(ti_black_holes_end_max, data->ti_black_holes_end_max);
    data->ti_black_holes_beg_max =
        max(ti_black_holes_beg_max, data->ti_black_holes_beg_max);
  }

  if (lock_unlock(&s->lock) != 0) error("Failed to unlock the space");
}

/**
 * @brief Collects the next time-step and rebuild flag.
 *
 * The next time-step is determined by making each super-cell recurse to
 * collect the minimal of ti_end and the number of updated particles.  When in
 * MPI mode this routines reduces these across all nodes and also collects the
 * forcerebuild flag -- this is so that we only use a single collective MPI
 * call per step for all these values.
 *
 * Note that the results are stored in e->collect_group1 struct not in the
 * engine fields, unless apply is true. These can be applied field-by-field
 * or all at once using collectgroup1_copy();
 *
 * @param e The #engine.
 * @param apply whether to apply the results to the engine or just keep in the
 *              group1 struct.
 */
void engine_collect_end_of_step(struct engine *e, int apply) {

  const ticks tic = getticks();
  struct space *s = e->s;
  struct end_of_step_data data;
  data.updated = 0, data.g_updated = 0, data.s_updated = 0, data.b_updated = 0;
  data.ti_hydro_end_min = max_nr_timesteps, data.ti_hydro_end_max = 0,
  data.ti_hydro_beg_max = 0;
  data.ti_gravity_end_min = max_nr_timesteps, data.ti_gravity_end_max = 0,
  data.ti_gravity_beg_max = 0;
  data.ti_stars_end_min = max_nr_timesteps, data.ti_stars_end_max = 0,
  data.ti_stars_beg_max = 0;
  data.ti_black_holes_end_min = max_nr_timesteps,
  data.ti_black_holes_end_max = 0, data.ti_black_holes_beg_max = 0;
  data.e = e;

  /* Initialize the total SFH of the simulation to zero */
  star_formation_logger_init(&data.sfh);

  /* Collect information from the local top-level cells */
  threadpool_map(&e->threadpool, engine_collect_end_of_step_mapper,
                 s->local_cells_with_tasks_top, s->nr_local_cells_with_tasks,
                 sizeof(int), 0, &data);

  /* Get the number of inhibited particles from the space-wide counters
   * since these have been updated atomically during the time-steps. */
  data.inhibited = s->nr_inhibited_parts;
  data.g_inhibited = s->nr_inhibited_gparts;
  data.s_inhibited = s->nr_inhibited_sparts;
  data.b_inhibited = s->nr_inhibited_bparts;

  /* Store these in the temporary collection group. */
  collectgroup1_init(
      &e->collect_group1, data.updated, data.g_updated, data.s_updated,
      data.b_updated, data.inhibited, data.g_inhibited, data.s_inhibited,
      data.b_inhibited, data.ti_hydro_end_min, data.ti_hydro_end_max,
      data.ti_hydro_beg_max, data.ti_gravity_end_min, data.ti_gravity_end_max,
      data.ti_gravity_beg_max, data.ti_stars_end_min, data.ti_stars_end_max,
      data.ti_stars_beg_max, data.ti_black_holes_end_min,
      data.ti_black_holes_end_max, data.ti_black_holes_beg_max, e->forcerebuild,
      e->s->tot_cells, e->sched.nr_tasks,
      (float)e->sched.nr_tasks / (float)e->s->tot_cells, data.sfh);

/* Aggregate collective data from the different nodes for this step. */
#ifdef WITH_MPI
  collectgroup1_reduce(&e->collect_group1);

#ifdef SWIFT_DEBUG_CHECKS
  {
    /* Check the above using the original MPI calls. */
    integertime_t in_i[2], out_i[2];
    in_i[0] = 0;
    in_i[1] = 0;
    out_i[0] = data.ti_hydro_end_min;
    out_i[1] = data.ti_gravity_end_min;
    if (MPI_Allreduce(out_i, in_i, 2, MPI_LONG_LONG_INT, MPI_MIN,
                      MPI_COMM_WORLD) != MPI_SUCCESS)
      error("Failed to aggregate ti_end_min.");
    if (in_i[0] != (long long)e->collect_group1.ti_hydro_end_min)
      error("Failed to get same ti_hydro_end_min, is %lld, should be %lld",
            in_i[0], e->collect_group1.ti_hydro_end_min);
    if (in_i[1] != (long long)e->collect_group1.ti_gravity_end_min)
      error("Failed to get same ti_gravity_end_min, is %lld, should be %lld",
            in_i[1], e->collect_group1.ti_gravity_end_min);

    long long in_ll[4], out_ll[4];
    out_ll[0] = data.updated;
    out_ll[1] = data.g_updated;
    out_ll[2] = data.s_updated;
    out_ll[3] = data.b_updated;
    if (MPI_Allreduce(out_ll, in_ll, 4, MPI_LONG_LONG_INT, MPI_SUM,
                      MPI_COMM_WORLD) != MPI_SUCCESS)
      error("Failed to aggregate particle counts.");
    if (in_ll[0] != (long long)e->collect_group1.updated)
      error("Failed to get same updated, is %lld, should be %lld", in_ll[0],
            e->collect_group1.updated);
    if (in_ll[1] != (long long)e->collect_group1.g_updated)
      error("Failed to get same g_updated, is %lld, should be %lld", in_ll[1],
            e->collect_group1.g_updated);
    if (in_ll[2] != (long long)e->collect_group1.s_updated)
      error("Failed to get same s_updated, is %lld, should be %lld", in_ll[2],
            e->collect_group1.s_updated);
    if (in_ll[3] != (long long)e->collect_group1.b_updated)
      error("Failed to get same b_updated, is %lld, should be %lld", in_ll[3],
            e->collect_group1.b_updated);

    out_ll[0] = data.inhibited;
    out_ll[1] = data.g_inhibited;
    out_ll[2] = data.s_inhibited;
    out_ll[3] = data.b_inhibited;
    if (MPI_Allreduce(out_ll, in_ll, 4, MPI_LONG_LONG_INT, MPI_SUM,
                      MPI_COMM_WORLD) != MPI_SUCCESS)
      error("Failed to aggregate particle counts.");
    if (in_ll[0] != (long long)e->collect_group1.inhibited)
      error("Failed to get same inhibited, is %lld, should be %lld", in_ll[0],
            e->collect_group1.inhibited);
    if (in_ll[1] != (long long)e->collect_group1.g_inhibited)
      error("Failed to get same g_inhibited, is %lld, should be %lld", in_ll[1],
            e->collect_group1.g_inhibited);
    if (in_ll[2] != (long long)e->collect_group1.s_inhibited)
      error("Failed to get same s_inhibited, is %lld, should be %lld", in_ll[2],
            e->collect_group1.s_inhibited);
    if (in_ll[3] != (long long)e->collect_group1.b_inhibited)
      error("Failed to get same b_inhibited, is %lld, should be %lld", in_ll[3],
            e->collect_group1.b_inhibited);

    int buff = 0;
    if (MPI_Allreduce(&e->forcerebuild, &buff, 1, MPI_INT, MPI_MAX,
                      MPI_COMM_WORLD) != MPI_SUCCESS)
      error("Failed to aggregate the rebuild flag across nodes.");
    if (!!buff != !!e->collect_group1.forcerebuild)
      error(
          "Failed to get same rebuild flag from all nodes, is %d,"
          "should be %d",
          buff, e->collect_group1.forcerebuild);
  }
#endif
#endif

  /* Apply to the engine, if requested. */
  if (apply) collectgroup1_apply(&e->collect_group1, e);

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}
