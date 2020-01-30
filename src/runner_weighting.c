/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Willem Elbers (whe@willemelbers.com)
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
#include "runner.h"

/* Phase space density functions needed */
#include "neutrino/phase_space.h"

/* Local headers. */
#include "active.h"
#include "cell.h"
#include "engine.h"
#include "feedback.h"
#include "kick.h"
#include "timers.h"
#include "timestep.h"
#include "timestep_limiter.h"
#include "timestep_sync.h"
#include "tracers.h"

/* For convert_gpart_vel() */
#include "gravity_io.h"

/**
 * @brief Weight the active neutrino particles in a cell to satisfy Liouville's
 * equation.
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param timer Are we timing this ?
 */
void runner_do_weighting(struct runner *r, struct cell *c, int timer) {

  const struct engine *e = r->e;
  const int with_cosmology = (e->policy & engine_policy_cosmology);
  struct gpart *restrict gparts = c->grav.parts;
  const int gcount = c->grav.count;

  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_starting_gravity(c, e) && !cell_is_active_gravity(c, e)) return;
  if (!with_cosmology)
    error("Phase space weighting without cosmology not implemented.");

  // const struct phys_const *physical_constants = e->physical_constants;
  // const struct cosmology *cosmo = e->cosmology;
  // const double volume = e->s->dim[0] * e->s->dim[1] * e->s->dim[2];
  // const double H_ratio = cosmo->H0 / cosmo->H;
  // const double rho_crit0 = cosmo->critical_density * H_ratio * H_ratio;
  // const double neutrino_mass = cosmo->Omega_nu * volume * rho_crit0;
  // const double particle_mass = neutrino_mass / e->total_nr_nuparts;
  // const double T_nu = cosmo->T_nu;
  // const double k_b = physical_constants->const_boltzmann_k;
  // const double eV = physical_constants->const_electron_volt;
  // const double T_eV = k_b * T_nu / eV; // temperature in eV

  /* Conversion factor from macro particle mass to neutrino mass in eV */
  const double mult = e->neutrino_mass_conversion_factor;

  /* Recurse? */
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) runner_do_weighting(r, c->progeny[k], 0);
  } else {

    /* Loop over the gparts in this cell. */
    for (int k = 0; k < gcount; k++) {
      /* Get a handle on the part. */
      struct gpart *restrict gp = &gparts[k];

      /* If the g-particle is a neutrino and needs to be weighted */
      if (gp->type == swift_type_neutrino && true) {
        if (gpart_is_active(gp, e)) {
          /* Set up the initial phase space density if necessary */
          if (e->step == 0) {
            gp->mass_i = gp->mass;
            gp->f_phase_i = fermi_dirac_density(e, gp->v_full, gp->mass * mult);
            gp->f_phase = gp->f_phase_i;
            gp->mass = 1e-12;  // dither in the first time step
          } else {
            gp->f_phase = fermi_dirac_density(e, gp->v_full, gp->mass_i * mult);
            gp->mass = gp->mass_i * (1.0 - gp->f_phase / gp->f_phase_i);

            // if (gp->id_or_neg_offset >= 262144 &&
            //     gp->id_or_neg_offset < 262144 + 5) {
            //   // double m = neutrino_mass_factor(e) * gp->mass_i;
            //   // double m2 = neutrino_mass_factor(e) * particle_mass;
            //   double p = fermi_dirac_momentum(e, gp->v_full, gp->mass_i *
            //   mult); message("%.10e %.10e %.10e %.10e %f", p, gp->mass_i *
            //   mult, gp->f_phase, gp->f_phase_i, gp->mass);
            // }
          }
        }
      }
    }

    if (timer) TIMER_TOC(timer_weight);
  }
}
