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
#ifndef SWIFT_GRAVITY_PROPERTIES
#define SWIFT_GRAVITY_PROPERTIES

/* Config parameters. */
#include "../config.h"

#if defined(HAVE_HDF5)
#include <hdf5.h>
#endif

/* Local includes. */
#include "kernel_gravity.h"
#include "restart.h"

/* Forward declarations */
struct cosmology;
struct phys_const;
struct swift_params;

/**
 * @brief Contains all the constants and parameters of the self-gravity scheme
 */
struct gravity_props {

  /* -------------- Softening for the regular particles ---------------- */

  /*! Co-moving softening length for high-res. DM particles at the current
   * redshift.  */
  float epsilon_DM_cur;

  /*! Co-moving softening length for high-res. baryon particles at the current
   * redshift.  */
  float epsilon_baryon_cur;

  /* -------------- Softening for the background DM -------------------- */

  /*! Conversion factor from cbrt of particle mass to softening assuming
   * a constant fraction of the mean inter-particle separation at that mass. */
  float epsilon_background_fac;

  /* -------------- Properties of the FFM gravity ---------------------- */

  /*! Tree opening angle (Multipole acceptance criterion) */
  double theta_crit;

  /*! Square of opening angle */
  double theta_crit2;

  /*! Inverse of opening angle */
  double theta_crit_inv;

  /* ------------- Properties of the softened gravity ------------------ */

  /*! Co-moving softening length for for high-res. DM particles */
  float epsilon_DM_comoving;

  /*! Maximal softening length in physical coordinates for the high-res.
   * DM particles */
  float epsilon_DM_max_physical;

  /*! Co-moving softening length for for high-res. baryon particles */
  float epsilon_baryon_comoving;

  /*! Maximal softening length in physical coordinates for the high-res.
   * baryon particles */
  float epsilon_baryon_max_physical;

  /*! Fraction of the mean inter particle separation corresponding to the
   * co-moving softening length of the low-res. particles (DM + baryons) */
  float mean_inter_particle_fraction_high_res;

  /* ------------- Properties of the time integration  ----------------- */

  /*! Frequency of tree-rebuild in units of #gpart updates. */
  float rebuild_frequency;

  /*! Time integration dimensionless multiplier */
  float eta;

  /* ------------- Properties of the mesh-based gravity ---------------- */

  /*! Periodic long-range mesh side-length */
  int mesh_size;

  /*! Mesh smoothing scale in units of top-level cell size */
  float a_smooth;

  /*! Distance below which the truncated mesh force is Newtonian in units of
   * a_smooth */
  float r_cut_min_ratio;

  /*! Distance above which the truncated mesh force is negligible in units of
   * a_smooth */
  float r_cut_max_ratio;

  /*! Gravitational constant (in internal units, copied from the physical
   * constants) */
  float G_Newton;
};

void gravity_props_print(const struct gravity_props *p);
void gravity_props_init(struct gravity_props *p, struct swift_params *params,
                        const struct phys_const *phys_const,
                        const struct cosmology *cosmo, const int with_cosmology,
                        const int has_baryons, const int has_DM,
                        const int is_zoom_simulation, const int periodic);
void gravity_props_update(struct gravity_props *p,
                          const struct cosmology *cosmo);

#if defined(HAVE_HDF5)
void gravity_props_print_snapshot(hid_t h_grpsph,
                                  const struct gravity_props *p);
#endif

/* Dump/restore. */
void gravity_props_struct_dump(const struct gravity_props *p, FILE *stream);
void gravity_props_struct_restore(struct gravity_props *p, FILE *stream);

#endif /* SWIFT_GRAVITY_PROPERTIES */
