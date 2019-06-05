/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_EAGLE_BH_IACT_H
#define SWIFT_EAGLE_BH_IACT_H

/* Local includes */
#include "hydro.h"
#include "random.h"

/**
 * @brief Density interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param bi First particle (black hole).
 * @param pj Second particle (gas, not updated).
 * @param xpj The extended data of the second particle (not updated).
 * @param cosmo The cosmological model.
 * @param ti_current Current integer time value (for random numbers).
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_bh_density(
    const float r2, const float *dx, const float hi, const float hj,
    struct bpart *restrict bi, const struct part *restrict pj,
    const struct xpart *restrict xpj, const struct cosmology *cosmo,
    const integertime_t ti_current) {

  float wi, wi_dx;

  /* Get r and 1/r. */
  const float r_inv = 1.0f / sqrtf(r2);
  const float r = r2 * r_inv;

  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  kernel_deval(ui, &wi, &wi_dx);

  /* Compute contribution to the number of neighbours */
  bi->density.wcount += wi;
  bi->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);

  /* Contribution to the number of neighbours */
  bi->num_ngbs += 1;

  /* Neighbour gas mass */
  const float mj = hydro_get_mass(pj);

  /* Get gas particle's (physical) density */
  const float rhoj = mj * wi * cosmo->a3_inv;
  
  /* Contribution to the BH gas density */
  bi->rho_gas += mj * wi;

  /* Contribution to the total neighbour mass */
  bi->ngb_mass += mj;

  /* Neighbour sounds speed */
  const float cj = hydro_get_comoving_soundspeed(pj);

  /* Contribution to the smoothed sound speed */
  bi->sound_speed_gas += mj * cj * wi;

  /* Neighbour peculiar drifted velocity */
  const float vj[3] = {pj->v[0], pj->v[1], pj->v[2]};

  /* Contribution to the smoothed velocity */
  bi->velocity_gas[0] += mj * vj[0] * wi;
  bi->velocity_gas[1] += mj * vj[1] * wi;
  bi->velocity_gas[2] += mj * vj[2] * wi;

  const double bh_v_peculiar[3] = {bi->v[0] * cosmo->a_inv,
                                   bi->v[1] * cosmo->a_inv,
                                   bi->v[2] * cosmo->a_inv};

  const double gas_v_peculiar[3] = {vj[0] * cosmo->a_inv,
                                    vj[1] * cosmo->a_inv,
                                    vj[2] * cosmo->a_inv};

  /* Difference in peculiar velocity between the gas and the BH
   * Note that there is no need for a Hubble flow term here. We are
   * computing the gas velocity at the position of the black hole. */
  const double v_diff_peculiar[3] = {gas_v_peculiar[0] - bh_v_peculiar[0],
                                     gas_v_peculiar[1] - bh_v_peculiar[1],
                                     gas_v_peculiar[2] - bh_v_peculiar[2]};

  /* (Bulk) velocity between BH and gas particle */  
  const double v_diff_norm2 = v_diff_peculiar[0] * v_diff_peculiar[0] +
                              v_diff_peculiar[1] * v_diff_peculiar[1] +
                              v_diff_peculiar[2] * v_diff_peculiar[2];
  
  /* Contribution to the accretion rate */
  const double gas_c_phys = cj * cosmo->a_factor_sound_speed;
  const double gas_c_phys2 = gas_c_phys * gas_c_phys;
  const double denominator2 = v_diff_norm2 + gas_c_phys2;
  const double denominator_inv = 1. / sqrt(denominator2);
  printf("Bondi-denominator-sq=%g, c_phys2=%g, rhoj=%g\n",
	 denominator2, gas_c_phys2, rhoj);
   
  /* 'Accretion rate' does not yet have the constant pre-factor */
  bi->accretion_rate += rhoj * denominator_inv * denominator_inv *
                        denominator_inv;
  /*printf("Accretion rate (w/o prefix)=%g\n", bi->accretion_rate); */
  
#ifdef DEBUG_INTERACTIONS_BH
  /* Update ngb counters */
  if (si->num_ngb_density < MAX_NUM_OF_NEIGHBOURS_BH)
    bi->ids_ngbs_density[si->num_ngb_density] = pj->id;

  /* Update ngb counters */
  ++si->num_ngb_density;
#endif
}

/**
 * @brief Feedback interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param bi First particle (black hole).
 * @param pj Second particle (gas)
 * @param xpj The extended data of the second particle.
 * @param cosmo The cosmological model.
 * @param ti_current Current integer time value (for random numbers).
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_bh_feedback(const float r2, const float *dx, const float hi,
                               const float hj, struct bpart *restrict bi,
                               struct part *restrict pj,
                               struct xpart *restrict xpj,
                               const struct cosmology *cosmo,
                               const integertime_t ti_current) {

  /* Get the heating probability */
  const float prob = bi->to_distribute.AGN_heating_probability;

  /* Are we doing some feedback? */
  if (prob > 0.f) {

    /* Draw a random number (Note mixing both IDs) */
    const float rand = random_unit_interval(bi->id + pj->id, ti_current,
                                            random_number_BH_feedback);

    /* Are we lucky? */
    if (rand < prob) {

      /* Compute new energy of this particle */
      const double u_init = hydro_get_physical_internal_energy(pj, xpj, cosmo);
      const float delta_u = bi->to_distribute.AGN_delta_u;
      const double u_new = u_init + delta_u;

      hydro_set_physical_internal_energy(pj, xpj, cosmo, u_new);
      hydro_set_drifted_physical_internal_energy(pj, cosmo, u_new);

      /* Impose maximal viscosity */
      hydro_set_viscosity_alpha_max_feedback(pj);
    }
  }

#ifdef DEBUG_INTERACTIONS_BH
  /* Update ngb counters */
  if (si->num_ngb_force < MAX_NUM_OF_NEIGHBOURS_BH)
    bi->ids_ngbs_force[si->num_ngb_force] = pj->id;

  /* Update ngb counters */
  ++si->num_ngb_force;
#endif
}

#endif /* SWIFT_EAGLE_BH_IACT_H */
