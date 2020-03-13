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

/* Calculate first order derivatives */
__attribute__((always_inline)) INLINE static void five_point_derivative_1(
    double mesh[6][6][6], double d[3], double tx, double ty, double tz,
    double dx, double dy, double dz, double fac) {

  const int ii = 2, jj = 2, kk = 2;

  /* 5-point stencil along each axis */
  d[0] = (1. / 12.) * CIC_get(mesh, ii + 2, jj, kk, tx, ty, tz, dx, dy, dz);
  d[0] -= (2. / 3.) * CIC_get(mesh, ii + 1, jj, kk, tx, ty, tz, dx, dy, dz);
  d[0] += (2. / 3.) * CIC_get(mesh, ii - 1, jj, kk, tx, ty, tz, dx, dy, dz);
  d[0] -= (1. / 12.) * CIC_get(mesh, ii - 2, jj, kk, tx, ty, tz, dx, dy, dz);

  d[1] = (1. / 12.) * CIC_get(mesh, ii, jj + 2, kk, tx, ty, tz, dx, dy, dz);
  d[1] -= (2. / 3.) * CIC_get(mesh, ii, jj + 1, kk, tx, ty, tz, dx, dy, dz);
  d[1] += (2. / 3.) * CIC_get(mesh, ii, jj - 1, kk, tx, ty, tz, dx, dy, dz);
  d[1] -= (1. / 12.) * CIC_get(mesh, ii, jj - 2, kk, tx, ty, tz, dx, dy, dz);

  d[2] = (1. / 12.) * CIC_get(mesh, ii, jj, kk + 2, tx, ty, tz, dx, dy, dz);
  d[2] -= (2. / 3.) * CIC_get(mesh, ii, jj, kk + 1, tx, ty, tz, dx, dy, dz);
  d[2] += (2. / 3.) * CIC_get(mesh, ii, jj, kk - 1, tx, ty, tz, dx, dy, dz);
  d[2] -= (1. / 12.) * CIC_get(mesh, ii, jj, kk - 2, tx, ty, tz, dx, dy, dz);

  /* Divide by the step size */
  d[0] *= fac;
  d[1] *= fac;
  d[2] *= fac;
}

/* Calculates second order derivatives, including mixed derivatives.
 * The order is 00,11,22,01,02,12 */
__attribute__((always_inline)) INLINE static void five_point_derivative_2(
    double mesh[6][6][6], double d[6], double tx, double ty, double tz,
    double dx, double dy, double dz, double fac) {

  const int ii = 2, jj = 2, kk = 2;

  /* 5-point stencil along each axis for the second derivatives */
  d[0] = -(1. / 12.) * CIC_get(mesh, ii + 2, jj, kk, tx, ty, tz, dx, dy, dz);
  d[0] += (4. / 3.) * CIC_get(mesh, ii + 1, jj, kk, tx, ty, tz, dx, dy, dz);
  d[0] -= (5. / 2.) * CIC_get(mesh, ii, jj, kk, tx, ty, tz, dx, dy, dz);
  d[0] += (4. / 3.) * CIC_get(mesh, ii - 1, jj, kk, tx, ty, tz, dx, dy, dz);
  d[0] -= (1. / 12.) * CIC_get(mesh, ii - 2, jj, kk, tx, ty, tz, dx, dy, dz);

  d[1] = -(1. / 12.) * CIC_get(mesh, ii, jj + 2, kk, tx, ty, tz, dx, dy, dz);
  d[1] += (4. / 3.) * CIC_get(mesh, ii, jj + 1, kk, tx, ty, tz, dx, dy, dz);
  d[1] -= (5. / 2.) * CIC_get(mesh, ii, jj, kk, tx, ty, tz, dx, dy, dz);
  d[1] += (4. / 3.) * CIC_get(mesh, ii, jj - 1, kk, tx, ty, tz, dx, dy, dz);
  d[1] -= (1. / 12.) * CIC_get(mesh, ii, jj - 2, kk, tx, ty, tz, dx, dy, dz);

  d[2] = -(1. / 12.) * CIC_get(mesh, ii, jj, kk + 2, tx, ty, tz, dx, dy, dz);
  d[2] += (4. / 3.) * CIC_get(mesh, ii, jj, kk + 1, tx, ty, tz, dx, dy, dz);
  d[2] -= (5. / 2.) * CIC_get(mesh, ii, jj, kk, tx, ty, tz, dx, dy, dz);
  d[2] += (4. / 3.) * CIC_get(mesh, ii, jj, kk - 1, tx, ty, tz, dx, dy, dz);
  d[2] -= (1. / 12.) * CIC_get(mesh, ii, jj, kk - 2, tx, ty, tz, dx, dy, dz);

  /* 3-point stencil along mixed axes for the mixed second derivatives */
  d[3] = (1. / 4.) * CIC_get(mesh, ii - 1, jj - 1, kk, tx, ty, tz, dx, dy, dz);
  d[3] += (1. / 4.) * CIC_get(mesh, ii + 1, jj + 1, kk, tx, ty, tz, dx, dy, dz);
  d[3] -= (1. / 4.) * CIC_get(mesh, ii + 1, jj - 1, kk, tx, ty, tz, dx, dy, dz);
  d[3] -= (1. / 4.) * CIC_get(mesh, ii - 1, jj + 1, kk, tx, ty, tz, dx, dy, dz);

  d[4] = (1. / 4.) * CIC_get(mesh, ii - 1, jj, kk - 1, tx, ty, tz, dx, dy, dz);
  d[4] += (1. / 4.) * CIC_get(mesh, ii + 1, jj, kk + 1, tx, ty, tz, dx, dy, dz);
  d[4] -= (1. / 4.) * CIC_get(mesh, ii + 1, jj, kk - 1, tx, ty, tz, dx, dy, dz);
  d[4] -= (1. / 4.) * CIC_get(mesh, ii - 1, jj, kk + 1, tx, ty, tz, dx, dy, dz);

  d[5] = (1. / 4.) * CIC_get(mesh, ii, jj - 1, kk - 1, tx, ty, tz, dx, dy, dz);
  d[5] += (1. / 4.) * CIC_get(mesh, ii, jj + 1, kk + 1, tx, ty, tz, dx, dy, dz);
  d[5] -= (1. / 4.) * CIC_get(mesh, ii, jj + 1, kk - 1, tx, ty, tz, dx, dy, dz);
  d[5] -= (1. / 4.) * CIC_get(mesh, ii, jj - 1, kk + 1, tx, ty, tz, dx, dy, dz);

  /* Divide by the step size (second derivative!) */
  d[0] *= fac * fac;
  d[1] *= fac * fac;
  d[2] *= fac * fac;
  d[3] *= fac * fac;
  d[4] *= fac * fac;
  d[5] *= fac * fac;
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
 * @param order The Legendre polynomial order .
 * @param direction The direction in which the derivatives should be taken.
 */
double grid_get(struct gpart *gp, const double *grid, int N,
                          double fac, const double dim[3], char order,
                          const double direction[3]) {

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

  /* The result */
  double p = 0;

  /* Indices of (i,j,k) in the local copy of the mesh */
  const int ii = 2, jj = 2, kk = 2;

  /* Simple CIC for the local grid itself */
  if (order == 0) {
    p += CIC_get(subgrid, ii, jj, kk, tx, ty, tz, dx, dy, dz);
  } else if (order == 1) {
    double derivative[3];
    five_point_derivative_1(subgrid, derivative, tx, ty, tz, dx, dy, dz, fac);

    p += direction[0]*derivative[0];
    p += direction[1]*derivative[1];
    p += direction[2]*derivative[2];
  } else if (order == 2) {
    double derivative[6]; //second derivatives {00,11,22,01,02,12}
    five_point_derivative_2(subgrid, derivative, tx, ty, tz, dx, dy, dz, fac);

    /* We want the Legendre polynomial P_2(x) = (3x^2 - 1)/2 */

    /* First the zeroth order contribution */
    // p -= 0.5 * CIC_get(subgrid, ii, jj, kk, tx, ty, tz, dx, dy, dz);

    // /* Determine the normalization */
    // double norm = derivative[0] + derivative[1] + derivative[2];
    //
    // norm = 1.0;

    /* Now the second order contributions */
    p += 1.5 * direction[0] * direction[0] * derivative[0];
    p += 1.5 * direction[1] * direction[1] * derivative[1];
    p += 1.5 * direction[2] * direction[2] * derivative[2];
    p += 3.0 * direction[0] * direction[1] * derivative[3];
    p += 3.0 * direction[0] * direction[2] * derivative[4];
    p += 3.0 * direction[1] * direction[2] * derivative[5];
  } else {
    error("Higher order derivatives not implemented.");
  }


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
  /* Locate the linear theory density grid */
  const int N = e->rend->primordial_grid_N;
  const double cell_fac = e->mesh->cell_fac;
  const double *grid = e->rend->density_grid;
  const double dim[3] = {e->s->dim[0], e->s->dim[1], e->s->dim[2]};
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

  double sum1=0,sum2=0,sum3=0;
  int count=0;

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

#ifdef NEUTRINO_DELTA_F_LINEAR_THEORY
          /* Determine the normalized direction of the momentum */
          double v[3];
          v[0] = gp->v_full[0];
          v[1] = gp->v_full[1];
          v[2] = gp->v_full[2];
          double vlen = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
          v[0] /= vlen;
          v[1] /= vlen;
          v[2] /= vlen;

          /* Correction to the distribution function from linear theory */
          double Psi = 0;

          /* The perturbation theory grids */
          double *pt = e->rend->the_grids;
          int size = N * N * N; //size of the grids

          /* Calculate the momentum in units of T_nu_0 */
          double momentum = fermi_dirac_momentum(e, gp->v_full, gp->mass * mult);
          double q = momentum / T_eV;

          /* Map to the Chebyshev interval [-1,1] assuming q_max = 15 */
          const double q_max = 15.0;
          double x = 2 * q / q_max - 1.0;

          if (x < -1 || x > 1) {
              message("WARNING: momentum exceeds allowed interval: q=%f",q);
          }

          /* Calculate the Chebyshev polynomials */
          double Cheb0 = 1.0;
          double Cheb1 = x;
          double Cheb2 = 2.0 * x * x - 1.0;

          /* 0th order contributions (Chebyshev 0,1,2) */
          Psi += grid_get(gp, pt + 3 * size, N, cell_fac, dim, 0, v) * Cheb0;
          Psi += grid_get(gp, pt + 4 * size, N, cell_fac, dim, 0, v) * Cheb1;
          Psi += grid_get(gp, pt + 5 * size, N, cell_fac, dim, 0, v) * Cheb2;
          /* 1st order contributions (Chebyshev 0,1,2) */
          Psi += grid_get(gp, pt + 6 * size, N, cell_fac, dim, 1, v) * Cheb0;
          Psi += grid_get(gp, pt + 7 * size, N, cell_fac, dim, 1, v) * Cheb1;
          Psi += grid_get(gp, pt + 8 * size, N, cell_fac, dim, 1, v) * Cheb2;
          /* 2nd order contributions (Chebyshev 0,1,2) */
          // Psi += grid_get(gp, pt + 9 * size, N, cell_fac, dim, 2, v) * Cheb0;
          // Psi += grid_get(gp, pt + 10 * size, N, cell_fac, dim, 2, v) * Cheb1;
          // Psi += grid_get(gp, pt + 11 * size, N, cell_fac, dim, 2, v) * Cheb2;

          double linear_overdensity =
              grid_get(gp, grid, N, cell_fac, dim, 0, v);
          double temperature_factor = cbrt(1.0 + linear_overdensity);
          double temperature_factor2 = cbrt(1.0 + linear_overdensity);
          temperature_factor = 1.0;
#else
          double temperature_factor = 1.0;
#endif
          /* Is it the first time step? */
          if (e->step == 0) {
            /* The mass of a microscopic neutrino in eV */
            double m_eV = gp->mass * mult;
            double f,g,h;

            /* Store the initial mass & phase space density */
            f = fermi_dirac_density(e, gp->v_full, m_eV, temperature_factor);
            f *= (1.0 + Psi);
            g = fermi_dirac_density(e, gp->v_full, m_eV, temperature_factor2);
            h = fermi_dirac_density(e, gp->v_full, m_eV, temperature_factor);
            gp->mass_i = gp->mass;
            gp->f_phase = f;
            gp->f_phase_i = f;
            gp->g_phase_i = g;
            gp->h_phase_i = h;
            // gp->mass = FLT_MIN;  // dither in the first time step
          } else {
            /* The mass of a microscopic neutrino in eV */
            double m_eV = gp->mass_i * mult;
            double f,g,h;

            /* Compute the phase space density */
            f = fermi_dirac_density(e, gp->v_full, m_eV, temperature_factor);
            gp->f_phase = f * (1.0 + Psi);

            g = fermi_dirac_density(e, gp->v_full, m_eV, cbrt(1.0 + linear_overdensity));

            h = fermi_dirac_density(e, gp->v_full, m_eV, temperature_factor);

            /* We use the energy instead of the mass: M -> sqrt(M^2 + P^2) */
            double energy_eV = fermi_dirac_energy(e, gp->v_full, m_eV);
            double energy = energy_eV / mult;  // energy in internal mass units

            /* Use the weighted energy instead of the mass */
            gp->mass = energy * (gp->f_phase_i - gp->f_phase) / gp->g_phase_i;
            gp->mass = gp->mass_i;

            double w = (gp->f_phase_i - gp->f_phase) / gp->f_phase_i;
            double w2 = (gp->g_phase_i - g) / gp->g_phase_i;
            double w3 = (gp->h_phase_i - h) / gp->h_phase_i;

            sum1+=w*w;
            sum2+=w2*w2;
            sum3+=w3*w3;
            count++;

            // if (gp->id_or_neg_offset >= 110592+10000 &&
            //     gp->id_or_neg_offset < 110592+10000 + 5) {
            //         message("%f %f %f %f", Psi, w, w2, w3);
            //         // message("%f %f %f %f %f %f", linear_overdensity, (f2-f)/f, Psi, (f - gp->f_phase_i)/gp->f_phase_i, (f - gp->f_phase_i)/gp->f_phase_i/Psi, x);
            //         // double v_len = sqrt(gp->v_full[0]*gp->v_full[0] +
            //         // gp->v_full[1]*gp->v_full[1] +
            //         // gp->v_full[2]*gp->v_full[2]); message("%f %f %f %f %f %f
            //         // %f", linear_overdensity, temperature_factor, f, energy,
            //         // gp->mass_i, energy/gp->mass_i, v_len);
            // }
          }
        }
      }
    }

    if (timer) TIMER_TOC(timer_weight);
  }

  if (count>26)
  message("%f %f %f", 0.5*sum1/count, 0.5*sum2/count, 0.5*sum3/count);
}
