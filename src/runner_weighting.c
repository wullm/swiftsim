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

#include "neutrino/delta_f.h"

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

#ifdef NEUTRINO_DELTA_F_LINEAR_THEORY
/**
 * @brief Returns 1D index of a 3D NxNxN array using row-major style.
 *
 * Wraps around in the corresponding dimension if any of the 3 indices is >= N
 * or < 0.
 *
 * @param i Index along x.
 * @param j Index along y.
 * @param k Index along z.
 * @param N Size of the array along one axis.
 */
__attribute__((always_inline)) INLINE static int row_major_id_periodic(int i,
                                                                       int j,
                                                                       int k,
                                                                       int N) {
  return (((i + N) % N) * N * N + ((j + N) % N) * N + ((k + N) % N));
}

/**
 * @brief Interpolate values from a the mesh using CIC.
 *
 * @param mesh The mesh to read from.
 * @param i The index of the cell along x
 * @param j The index of the cell along y
 * @param k The index of the cell along z
 * @param tx First CIC coefficient along x
 * @param ty First CIC coefficient along y
 * @param tz First CIC coefficient along z
 * @param dx Second CIC coefficient along x
 * @param dy Second CIC coefficient along y
 * @param dz Second CIC coefficient along z
 */
__attribute__((always_inline)) INLINE static double CIC_get(
    double mesh[6][6][6], int i, int j, int k, double tx, double ty, double tz,
    double dx, double dy, double dz) {

  double temp;
  temp = mesh[i + 0][j + 0][k + 0] * tx * ty * tz;
  temp += mesh[i + 0][j + 0][k + 1] * tx * ty * dz;
  temp += mesh[i + 0][j + 1][k + 0] * tx * dy * tz;
  temp += mesh[i + 0][j + 1][k + 1] * tx * dy * dz;
  temp += mesh[i + 1][j + 0][k + 0] * dx * ty * tz;
  temp += mesh[i + 1][j + 0][k + 1] * dx * ty * dz;
  temp += mesh[i + 1][j + 1][k + 0] * dx * dy * tz;
  temp += mesh[i + 1][j + 1][k + 1] * dx * dy * dz;

  return temp;
}

/**
 * @brief Retrieve value for a gpart from a given mesh using the CIC
 * method.
 *
 * Debugging routine.
 *
 * @param gp The #gpart.
 * @param grid The grid.
 * @param N the size of the mesh along one axis.
 * @param fac width of a mesh cell.
 * @param dim The dimensions of the simulation box.
 */
double grid_to_gparts_CIC(struct gpart *gp, const double *grid, int N,
                          double fac, const double dim[3]) {

  /* Box wrap the gpart's position */
  const double pos_x = box_wrap(gp->x[0], 0., dim[0]);
  const double pos_y = box_wrap(gp->x[1], 0., dim[1]);
  const double pos_z = box_wrap(gp->x[2], 0., dim[2]);

  int i = (int)(fac * pos_x);
  if (i >= N) i = N - 1;
  const double dx = fac * pos_x - i;
  const double tx = 1. - dx;

  int j = (int)(fac * pos_y);
  if (j >= N) j = N - 1;
  const double dy = fac * pos_y - j;
  const double ty = 1. - dy;

  int k = (int)(fac * pos_z);
  if (k >= N) k = N - 1;
  const double dz = fac * pos_z - k;
  const double tz = 1. - dz;

#ifdef SWIFT_DEBUG_CHECKS
  if (i < 0 || i >= N) error("Invalid gpart position in x");
  if (j < 0 || j >= N) error("Invalid gpart position in y");
  if (k < 0 || k >= N) error("Invalid gpart position in z");
#endif

#ifdef SWIFT_GRAVITY_FORCE_CHECKS
  if (gp->a_grav_PM[0] != 0. || gp->potential_PM != 0.)
    error("Particle with non-initalised stuff");
#endif

  /* First, copy the necessary part of the mesh for stencil operations */
  /* This includes box-wrapping in all 3 dimensions. */
  double subgrid[6][6][6];
  for (int iii = -2; iii <= 3; ++iii) {
    for (int jjj = -2; jjj <= 3; ++jjj) {
      for (int kkk = -2; kkk <= 3; ++kkk) {
        subgrid[iii + 2][jjj + 2][kkk + 2] =
            grid[row_major_id_periodic(i + iii, j + jjj, k + kkk, N)];
      }
    }
  }

  /* Some local accumulators */
  double p = 0.;

  /* Indices of (i,j,k) in the local copy of the mesh */
  const int ii = 2, jj = 2, kk = 2;

  /* Simple CIC for the local grid itself */
  p += CIC_get(subgrid, ii, jj, kk, tx, ty, tz, dx, dy, dz);

  /* ---- */
  return p;
}
#endif

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

#ifdef NEUTRINO_DELTA_F_LINEAR_THEORY
#ifdef RENDERER_USED
  /* Locate the linear theory density grid */
  const int N = e->rend->primordial_grid_N;
  const double cell_fac = e->mesh->cell_fac;
  const double *grid = e->rend->density_grid;
  const double dim[3] = {e->s->dim[0], e->s->dim[1], e->s->dim[2]};
#else
  error("Running with linear theory delta-f, but no renderer.");
#endif
#endif

  const struct phys_const *physical_constants = e->physical_constants;
  const struct cosmology *cosmo = e->cosmology;
  // const double volume = e->s->dim[0] * e->s->dim[1] * e->s->dim[2];
  // const double H_ratio = cosmo->H0 / cosmo->H;
  // const double rho_crit0 = cosmo->critical_density * H_ratio * H_ratio;
  // const double neutrino_mass = cosmo->Omega_nu * volume * rho_crit0;
  // const double particle_mass = neutrino_mass / e->total_nr_nuparts;
  const double T_nu = cosmo->T_nu;
  const double k_b = physical_constants->const_boltzmann_k;
  const double eV = physical_constants->const_electron_volt;
  const double T_eV = k_b * T_nu / eV; // temperature in eV

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

      if (e->step == 0) {
        gp->v_i = sqrt(gp->v_full[0] * gp->v_full[0] +
                       gp->v_full[1] * gp->v_full[1] +
                       gp->v_full[2] * gp->v_full[2]);
      }

      /* If the g-particle is a neutrino and needs to be weighted */
      if (gp->type == swift_type_neutrino && true) {
        if (gpart_is_active(gp, e)) {

          double temperature_factor = 1.0;

#ifdef NEUTRINO_DELTA_F_LINEAR_THEORY
#ifdef RENDERER_USED
          double overdensity = grid_to_gparts_CIC(gp, grid, N, cell_fac, dim);
          temperature_factor = cbrt(1.0 + overdensity);
#else
          error("Running with linear theory delta-f, but no renderer.");
#endif
#endif
          /* Is it the first time step? */
          if (e->step == 0) {
            /* The mass of a microscopic neutrino in eV */
            double m_eV = gp->mass * mult;
            double f;
            (void) m_eV;

            /* Use the particle ID as seed */
            uint64_t seed = gp->id_or_neg_offset;

            /* A unique uniform random number for this neutrino */
            const double w = sampleUniform(&seed);
            /* The corresponding initial Fermi-Dirac momentum */
            const double p_eV = fermi_dirac_transform(w) * T_eV;
            /* The corresponding initial density */
            const double f_i = fermi_dirac_density2(p_eV, T_eV);

            /* Store the initial mass & phase space density */
            f = f_i;
            gp->mass_i = gp->mass;
            gp->f_phase = f;
            gp->f_phase_i = f;
            gp->mass = FLT_MIN;  // dither in the first time step
          } else {
            /* The mass of a microscopic neutrino in eV */
            double m_eV = gp->mass_i * mult;
            double f;

            /* Compute the phase space density */
            f = fermi_dirac_density(e, gp->v_full, m_eV, temperature_factor);
            gp->f_phase = f;

            /* We use the energy instead of the mass: M -> sqrt(M^2 + P^2) */
            double energy_eV = fermi_dirac_energy(e, gp->v_i, m_eV);
            double energy = energy_eV / mult;  // energy in internal mass units

            /* Use the weighted energy instead of the mass */
            gp->mass = energy * (1.0 - gp->f_phase / gp->f_phase_i);

            /* Avoid poles */
            if (gp->mass == 0) {
              gp->mass = FLT_MIN;
            }

            // if (gp->id_or_neg_offset >= 114688-4096 &&
            //     gp->id_or_neg_offset < 114688-4096 + 5) {
            //         message("%f %f %f %f", linear_overdensity,
            //         temperature_factor, f, energy);
            // //         double p = fermi_dirac_momentum(e, gp->v_full, m_eV) /
            // //         e->cosmology->a; message("%.10e %.10e %.10e %f", p,
            // m_eV,
            // //         energy_eV, energy / gp->mass_i);
            //     }
          }
        }
      }
    }

    if (timer) TIMER_TOC(timer_weight);
  }
}
