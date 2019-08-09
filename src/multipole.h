/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2013 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *               2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_MULTIPOLE_H
#define SWIFT_MULTIPOLE_H

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <math.h>
#include <string.h>

/* Includes. */
#include "align.h"
#include "const.h"
#include "error.h"
#include "gravity.h"
#include "gravity_derivatives.h"
#include "gravity_properties.h"
#include "gravity_softened_derivatives.h"
#include "inline.h"
#include "kernel_gravity.h"
#include "part.h"
#include "periodic.h"
#include "vector_power.h"

#define multipole_align 128

struct grav_tensor {

  /* 0th order terms */
  float F_000;

#if SELF_GRAVITY_MULTIPOLE_ORDER > 0

  /* 1st order terms */
  float F_100, F_010, F_001;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 1

  /* 2nd order terms */
  float F_200, F_020, F_002;
  float F_110, F_101, F_011;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 2

  /* 3rd order terms */
  float F_300, F_030, F_003;
  float F_210, F_201;
  float F_120, F_021;
  float F_102, F_012;
  float F_111;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 3

  /* 4th order terms */
  float F_400, F_040, F_004;
  float F_310, F_301;
  float F_130, F_031;
  float F_103, F_013;
  float F_220, F_202, F_022;
  float F_211, F_121, F_112;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 4

  /* 5th order terms */
  float F_005, F_014, F_023;
  float F_032, F_041, F_050;
  float F_104, F_113, F_122;
  float F_131, F_140, F_203;
  float F_212, F_221, F_230;
  float F_302, F_311, F_320;
  float F_401, F_410, F_500;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 5

  /* 6th order terms */
  float F_006, F_015, F_024;
  float F_033, F_042, F_051;
  float F_060, F_105, F_114;
  float F_123, F_132, F_141;
  float F_150, F_204, F_213;
  float F_222, F_231, F_240;
  float F_303, F_312, F_321;
  float F_330, F_402, F_411;
  float F_420, F_501, F_510;
#endif 
#if SELF_GRAVITY_MULTIPOLE_ORDER > 6
#error "Missing implementation for order >6"
#endif

#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_GRAVITY_FORCE_CHECKS)
  /* Total number of gpart this field tensor interacted with */
  long long num_interacted;

  /* Total number of gpart this field tensor did not interact
   * with. i.e., the distance to the cell was > r_cut_max. */
  long long num_not_interacted;
#endif

#ifdef SWIFT_DEBUG_CHECKS

  /* Last time this tensor was zeroed */
  integertime_t ti_init;

#endif

  /* Has this tensor received any contribution? */
  char interacted;
};

struct multipole {

  /*! Bulk velocity */
  float vel[3];

  /*! Maximal velocity along each axis of all #gpart */
  float max_delta_vel[3];

  /*! Minimal velocity along each axis of all #gpart */
  float min_delta_vel[3];

  /*! Maximal co-moving softening of all the #gpart in the mulipole */
  float max_softening;

#ifdef ADVANCED_OPENING_CRITERIA
  /* Array to store each order of the multipole power */
  float power[SELF_GRAVITY_MULTIPOLE_ORDER+1]; 
  
  /* Minimum (norm of the) acceleration of the gparts in the cell */
  float min_a_grav_norm;
#endif

  /* 0th order term */
  float M_000;

#if SELF_GRAVITY_MULTIPOLE_ORDER > 0

  /* 1st order terms */
  float M_100, M_010, M_001;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 1

  /* 2nd order terms */
  float M_200, M_020, M_002;
  float M_110, M_101, M_011;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 2

  /* 3rd order terms */
  float M_300, M_030, M_003;
  float M_210, M_201;
  float M_120, M_021;
  float M_102, M_012;
  float M_111;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 3

  /* 4th order terms */
  float M_400, M_040, M_004;
  float M_310, M_301;
  float M_130, M_031;
  float M_103, M_013;
  float M_220, M_202, M_022;
  float M_211, M_121, M_112;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 4

  /* 5th order terms */
  float M_005, M_014, M_023;
  float M_032, M_041, M_050;
  float M_104, M_113, M_122;
  float M_131, M_140, M_203;
  float M_212, M_221, M_230;
  float M_302, M_311, M_320;
  float M_401, M_410, M_500;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 5

  /* 6th order terms */
  float M_006, M_015, M_024;
  float M_033, M_042, M_051;
  float M_060, M_105, M_114;
  float M_123, M_132, M_141;
  float M_150, M_204, M_213;
  float M_222, M_231, M_240;
  float M_303, M_312, M_321;
  float M_330, M_402, M_411;
  float M_420, M_501, M_510;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 6
#error "Missing implementation for order >6"
#endif

#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_GRAVITY_FORCE_CHECKS)

  /* Total number of gpart in this multipole */
  long long num_gpart;

#endif
};

/**
 * @brief The multipole expansion of a mass distribution.
 */
struct gravity_tensors {

  union {

    /*! Linking pointer for "memory management". */
    struct gravity_tensors *next;

    /*! The actual content */
    struct {

      /*! Field tensor for the potential */
      struct grav_tensor pot;

      /*! Multipole mass */
      struct multipole m_pole;

      /*! Centre of mass of the matter dsitribution */
      double CoM[3];

      /*! Centre of mass of the matter dsitribution at the last rebuild */
      double CoM_rebuild[3];

      /*! Upper limit of the CoM<->gpart distance */
      double r_max;

      /*! Upper limit of the CoM<->gpart distance at the last rebuild */
      double r_max_rebuild;
    };
  };
} SWIFT_STRUCT_ALIGN;

#ifdef WITH_MPI
/* MPI datatypes for transfers */
extern MPI_Datatype multipole_mpi_type;
extern MPI_Op multipole_mpi_reduce_op;
void multipole_create_mpi_types(void);
#endif

/**
 * @brief Reset the data of a #multipole.
 *
 * @param m The #multipole.
 */
INLINE static void gravity_reset(struct gravity_tensors *m) {

  /* Just bzero the struct. */
  bzero(m, sizeof(struct gravity_tensors));
}

/**
 * @brief Drifts a #multipole forward in time.
 *
 * This uses a first-order approximation in time. We only move the CoM
 * using the bulk velocity measured at the last rebuild.
 *
 * @param m The #multipole.
 * @param dt The drift time-step.
 */
INLINE static void gravity_drift(struct gravity_tensors *m, double dt) {

  /* Motion of the centre of mass */
  const double dx = m->m_pole.vel[0] * dt;
  const double dy = m->m_pole.vel[1] * dt;
  const double dz = m->m_pole.vel[2] * dt;

  /* Move the whole thing according to bulk motion */
  m->CoM[0] += dx;
  m->CoM[1] += dy;
  m->CoM[2] += dz;

#ifdef SWIFT_DEBUG_CHECKS
  if (m->m_pole.vel[0] > m->m_pole.max_delta_vel[0] * 1.1)
    error("Invalid maximal velocity");
  if (m->m_pole.vel[0] < m->m_pole.min_delta_vel[0] * 1.1)
    error("Invalid minimal velocity");
  if (m->m_pole.vel[1] > m->m_pole.max_delta_vel[1] * 1.1)
    error("Invalid maximal velocity");
  if (m->m_pole.vel[1] < m->m_pole.min_delta_vel[1] * 1.1)
    error("Invalid minimal velocity");
  if (m->m_pole.vel[2] > m->m_pole.max_delta_vel[2] * 1.1)
    error("Invalid maximal velocity");
  if (m->m_pole.vel[2] < m->m_pole.min_delta_vel[2] * 1.1)
    error("Invalid minimal velocity");
#endif

  /* Maximal distance covered by any particle */
  float dv[3];
  dv[0] = max(m->m_pole.max_delta_vel[0] - m->m_pole.vel[0],
              m->m_pole.vel[0] - m->m_pole.min_delta_vel[0]);
  dv[1] = max(m->m_pole.max_delta_vel[1] - m->m_pole.vel[1],
              m->m_pole.vel[1] - m->m_pole.min_delta_vel[1]);
  dv[2] = max(m->m_pole.max_delta_vel[2] - m->m_pole.vel[2],
              m->m_pole.vel[2] - m->m_pole.min_delta_vel[2]);

  const float max_delta_vel =
      sqrt(dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]);
  const float x_diff = max_delta_vel * dt;

  /* Conservative change in maximal radius containing all gpart */
  m->r_max += x_diff;
}

/**
 * @brief Zeroes all the fields of a field tensor
 *
 * @param l The field tensor.
 * @param ti_current The current (integer) time (for debugging only).
 */
INLINE static void gravity_field_tensors_init(struct grav_tensor *l,
                                              integertime_t ti_current) {

  bzero(l, sizeof(struct grav_tensor));

#ifdef SWIFT_DEBUG_CHECKS
  l->ti_init = ti_current;
#endif
}

/**
 * @brief Adds a field tensor to another one (i.e. does la += lb).
 *
 * @param la The gravity tensors to add to.
 * @param lb The gravity tensors to add.
 */
INLINE static void gravity_field_tensors_add(
    struct grav_tensor *restrict la, const struct grav_tensor *restrict lb) {
#ifdef SWIFT_DEBUG_CHECKS
  if (lb->num_interacted + lb->num_not_interacted == 0)
      error("Adding tensors that did not interact");
#endif
#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_GRAVITY_FORCE_CHECKS)
  la->num_interacted += lb->num_interacted;
  la->num_not_interacted += lb->num_not_interacted;
#endif
  la->interacted = 1;

  /* Add 0th order terms */
  la->F_000 += lb->F_000;

#if SELF_GRAVITY_MULTIPOLE_ORDER > 0
  /* Add 1st order terms */
  la->F_100 += lb->F_100;
  la->F_010 += lb->F_010;
  la->F_001 += lb->F_001;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 1
  /* Add 2nd order terms */
  la->F_200 += lb->F_200;
  la->F_020 += lb->F_020;
  la->F_002 += lb->F_002;
  la->F_110 += lb->F_110;
  la->F_101 += lb->F_101;
  la->F_011 += lb->F_011;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 2
  /* Add 3rd order terms */
  la->F_300 += lb->F_300;
  la->F_030 += lb->F_030;
  la->F_003 += lb->F_003;
  la->F_210 += lb->F_210;
  la->F_201 += lb->F_201;
  la->F_120 += lb->F_120;
  la->F_021 += lb->F_021;
  la->F_102 += lb->F_102;
  la->F_012 += lb->F_012;
  la->F_111 += lb->F_111;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 3
  /* Add 4th order terms */
  la->F_400 += lb->F_400;
  la->F_040 += lb->F_040;
  la->F_004 += lb->F_004;
  la->F_310 += lb->F_310;
  la->F_301 += lb->F_301;
  la->F_130 += lb->F_130;
  la->F_031 += lb->F_031;
  la->F_103 += lb->F_103;
  la->F_013 += lb->F_013;
  la->F_220 += lb->F_220;
  la->F_202 += lb->F_202;
  la->F_022 += lb->F_022;
  la->F_211 += lb->F_211;
  la->F_121 += lb->F_121;
  la->F_112 += lb->F_112;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 4
  /* 5th order terms */
  la->F_005 += lb->F_005;
  la->F_014 += lb->F_014;
  la->F_023 += lb->F_023;
  la->F_032 += lb->F_032;
  la->F_041 += lb->F_041;
  la->F_050 += lb->F_050;
  la->F_104 += lb->F_104;
  la->F_113 += lb->F_113;
  la->F_122 += lb->F_122;
  la->F_131 += lb->F_131;
  la->F_140 += lb->F_140;
  la->F_203 += lb->F_203;
  la->F_212 += lb->F_212;
  la->F_221 += lb->F_221;
  la->F_230 += lb->F_230;
  la->F_302 += lb->F_302;
  la->F_311 += lb->F_311;
  la->F_320 += lb->F_320;
  la->F_401 += lb->F_401;
  la->F_410 += lb->F_410;
  la->F_500 += lb->F_500;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 5
  /* 6th order terms */
  la->F_006 += lb->F_006;
  la->F_015 += lb->F_015;
  la->F_024 += lb->F_024;
  la->F_033 += lb->F_033;
  la->F_042 += lb->F_042;
  la->F_051 += lb->F_051;
  la->F_060 += lb->F_060;
  la->F_105 += lb->F_105;
  la->F_114 += lb->F_114;
  la->F_123 += lb->F_123;
  la->F_132 += lb->F_132;
  la->F_141 += lb->F_141;
  la->F_150 += lb->F_150;
  la->F_204 += lb->F_204;
  la->F_213 += lb->F_213;
  la->F_222 += lb->F_222;
  la->F_231 += lb->F_231;
  la->F_240 += lb->F_240;
  la->F_303 += lb->F_303;
  la->F_312 += lb->F_312;
  la->F_321 += lb->F_321;
  la->F_330 += lb->F_330;
  la->F_402 += lb->F_402;
  la->F_411 += lb->F_411;
  la->F_420 += lb->F_420;
  la->F_501 += lb->F_501;
  la->F_510 += lb->F_510;
  la->F_600 += lb->F_600;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 6
#error "Missing implementation for order >6"
#endif
}

/**
 * @brief Prints the content of a #grav_tensor to stdout.
 *
 * Note: Uses directly printf(), not a call to message().
 *
 * @param l The #grav_tensor to print.
 */
INLINE static void gravity_field_tensors_print(const struct grav_tensor *l) {

  printf("-------------------------\n");
  printf("Interacted: %d\n", l->interacted);
  printf("F_000= %12.5e\n", l->F_000);
#if SELF_GRAVITY_MULTIPOLE_ORDER > 0
  printf("-------------------------\n");
  printf("F_100= %12.5e F_010= %12.5e F_001= %12.5e\n", l->F_100, l->F_010,
         l->F_001);
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 1
  printf("-------------------------\n");
  printf("F_200= %12.5e F_020= %12.5e F_002= %12.5e\n", l->F_200, l->F_020,
         l->F_002);
  printf("F_110= %12.5e F_101= %12.5e F_011= %12.5e\n", l->F_110, l->F_101,
         l->F_011);
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 2
  printf("-------------------------\n");
  printf("F_300= %12.5e F_030= %12.5e F_003= %12.5e\n", l->F_300, l->F_030,
         l->F_003);
  printf("F_210= %12.5e F_201= %12.5e F_120= %12.5e\n", l->F_210, l->F_201,
         l->F_120);
  printf("F_021= %12.5e F_102= %12.5e F_012= %12.5e\n", l->F_021, l->F_102,
         l->F_012);
  printf("F_111= %12.5e\n", l->F_111);
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 3
  printf("-------------------------\n");
  printf("F_400= %12.5e F_040= %12.5e F_004= %12.5e\n", l->F_400, l->F_040,
         l->F_004);
  printf("F_310= %12.5e F_301= %12.5e F_130= %12.5e\n", l->F_310, l->F_301,
         l->F_130);
  printf("F_031= %12.5e F_103= %12.5e F_013= %12.5e\n", l->F_031, l->F_103,
         l->F_013);
  printf("F_220= %12.5e F_202= %12.5e F_022= %12.5e\n", l->F_220, l->F_202,
         l->F_022);
  printf("F_211= %12.5e F_121= %12.5e F_112= %12.5e\n", l->F_211, l->F_121,
         l->F_112);
#endif
  printf("-------------------------\n");
#if SELF_GRAVITY_MULTIPOLE_ORDER > 4
#error "Missing gravity_field_tensors_print() implementation for order >4"
#endif
}

/**
 * @brief Zeroes all the fields of a multipole.
 *
 * @param m The multipole
 */
INLINE static void gravity_multipole_init(struct multipole *m) {

  bzero(m, sizeof(struct multipole));
}

/**
 * @brief Prints the content of a #multipole to stdout.
 *
 * Note: Uses directly printf(), not a call to message().
 *
 * @param m The #multipole to print.
 */
INLINE static void gravity_multipole_print(const struct multipole *m) {

  printf("eps_max = %12.5e\n", m->max_softening);
  printf("Vel= [%12.5e %12.5e %12.5e]\n", m->vel[0], m->vel[1], m->vel[2]);
  printf("-------------------------\n");
  printf("M_000= %12.5e\n", m->M_000);
#if SELF_GRAVITY_MULTIPOLE_ORDER > 0
  printf("-------------------------\n");
  printf("M_100= %12.5e M_010= %12.5e M_001= %12.5e\n", m->M_100, m->M_010,
         m->M_001);
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 1
  printf("-------------------------\n");
  printf("M_200= %12.5e M_020= %12.5e M_002= %12.5e\n", m->M_200, m->M_020,
         m->M_002);
  printf("M_110= %12.5e M_101= %12.5e M_011= %12.5e\n", m->M_110, m->M_101,
         m->M_011);
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 2
  printf("-------------------------\n");
  printf("M_300= %12.5e M_030= %12.5e M_003= %12.5e\n", m->M_300, m->M_030,
         m->M_003);
  printf("M_210= %12.5e M_201= %12.5e M_120= %12.5e\n", m->M_210, m->M_201,
         m->M_120);
  printf("M_021= %12.5e M_102= %12.5e M_012= %12.5e\n", m->M_021, m->M_102,
         m->M_012);
  printf("M_111= %12.5e\n", m->M_111);
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 3
  printf("-------------------------\n");
  printf("M_400= %12.5e M_040= %12.5e M_004= %12.5e\n", m->M_400, m->M_040,
         m->M_004);
  printf("M_310= %12.5e M_301= %12.5e M_130= %12.5e\n", m->M_310, m->M_301,
         m->M_130);
  printf("M_031= %12.5e M_103= %12.5e M_013= %12.5e\n", m->M_031, m->M_103,
         m->M_013);
  printf("M_220= %12.5e M_202= %12.5e M_022= %12.5e\n", m->M_220, m->M_202,
         m->M_022);
  printf("M_211= %12.5e M_121= %12.5e M_112= %12.5e\n", m->M_211, m->M_121,
         m->M_112);
#endif
  printf("-------------------------\n");
#if SELF_GRAVITY_MULTIPOLE_ORDER > 4
#error "Missing gravity_multipole_print() implementation for order >4"
#endif
}

/**
 * @brief Adds a #multipole to another one (i.e. does ma += mb).
 *
 * @param ma The multipole to add to.
 * @param mb The multipole to add.
 */
INLINE static void gravity_multipole_add(struct multipole *restrict ma,
                                         const struct multipole *restrict mb) {

  /* Maximum of both softenings */
  ma->max_softening = max(ma->max_softening, mb->max_softening);

  /* Add 0th order term */
  ma->M_000 += mb->M_000;

#if SELF_GRAVITY_MULTIPOLE_ORDER > 0
  /* Add 1st order terms */
  ma->M_100 += mb->M_100;
  ma->M_010 += mb->M_010;
  ma->M_001 += mb->M_001;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 1
  /* Add 2nd order terms */
  ma->M_200 += mb->M_200;
  ma->M_020 += mb->M_020;
  ma->M_002 += mb->M_002;
  ma->M_110 += mb->M_110;
  ma->M_101 += mb->M_101;
  ma->M_011 += mb->M_011;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 2
  /* Add 3rd order terms */
  ma->M_300 += mb->M_300;
  ma->M_030 += mb->M_030;
  ma->M_003 += mb->M_003;
  ma->M_210 += mb->M_210;
  ma->M_201 += mb->M_201;
  ma->M_120 += mb->M_120;
  ma->M_021 += mb->M_021;
  ma->M_102 += mb->M_102;
  ma->M_012 += mb->M_012;
  ma->M_111 += mb->M_111;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 3
  /* Add 4th order terms */
  ma->M_400 += mb->M_400;
  ma->M_040 += mb->M_040;
  ma->M_004 += mb->M_004;
  ma->M_310 += mb->M_310;
  ma->M_301 += mb->M_301;
  ma->M_130 += mb->M_130;
  ma->M_031 += mb->M_031;
  ma->M_103 += mb->M_103;
  ma->M_013 += mb->M_013;
  ma->M_220 += mb->M_220;
  ma->M_202 += mb->M_202;
  ma->M_022 += mb->M_022;
  ma->M_211 += mb->M_211;
  ma->M_121 += mb->M_121;
  ma->M_112 += mb->M_112;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 4
  /* 5th order terms */
  ma->M_005 += mb->M_005;
  ma->M_014 += mb->M_014;
  ma->M_023 += mb->M_023;
  ma->M_032 += mb->M_032;
  ma->M_041 += mb->M_041;
  ma->M_050 += mb->M_050;
  ma->M_104 += mb->M_104;
  ma->M_113 += mb->M_113;
  ma->M_122 += mb->M_122;
  ma->M_131 += mb->M_131;
  ma->M_140 += mb->M_140;
  ma->M_203 += mb->M_203;
  ma->M_212 += mb->M_212;
  ma->M_221 += mb->M_221;
  ma->M_230 += mb->M_230;
  ma->M_302 += mb->M_302;
  ma->M_311 += mb->M_311;
  ma->M_320 += mb->M_320;
  ma->M_401 += mb->M_401;
  ma->M_410 += mb->M_410;
  ma->M_500 += mb->M_500;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 5
  /* 6th order terms */
  ma->M_006 += mb->M_006;
  ma->M_015 += mb->M_015;
  ma->M_024 += mb->M_024;
  ma->M_033 += mb->M_033;
  ma->M_042 += mb->M_042;
  ma->M_051 += mb->M_051;
  ma->M_060 += mb->M_060;
  ma->M_105 += mb->M_105;
  ma->M_114 += mb->M_114;
  ma->M_123 += mb->M_123;
  ma->M_132 += mb->M_132;
  ma->M_141 += mb->M_141;
  ma->M_150 += mb->M_150;
  ma->M_204 += mb->M_204;
  ma->M_213 += mb->M_213;
  ma->M_222 += mb->M_222;
  ma->M_231 += mb->M_231;
  ma->M_240 += mb->M_240;
  ma->M_303 += mb->M_303;
  ma->M_312 += mb->M_312;
  ma->M_321 += mb->M_321;
  ma->M_330 += mb->M_330;
  ma->M_402 += mb->M_402;
  ma->M_411 += mb->M_411;
  ma->M_420 += mb->M_420;
  ma->M_501 += mb->M_501;
  ma->M_510 += mb->M_510;
  ma->M_600 += mb->M_600;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 6
#error "Missing implementation for order >6"
#endif

#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_GRAVITY_FORCE_CHECKS)
  ma->num_gpart += mb->num_gpart;
#endif
}

/**
 * @brief Verifies whether two #multipole's are equal or not.
 *
 * @param ga The first #multipole.
 * @param gb The second #multipole.
 * @param tolerance The maximal allowed relative difference for the fields.
 * @return 1 if the multipoles are equal, 0 otherwise
 */
INLINE static int gravity_multipole_equal(const struct gravity_tensors *ga,
                                          const struct gravity_tensors *gb,
                                          double tolerance) {

  /* Check CoM */
  if (fabs(ga->CoM[0] - gb->CoM[0]) / fabs(ga->CoM[0] + gb->CoM[0]) >
      tolerance) {
    message("CoM[0] different");
    return 0;
  }
  if (fabs(ga->CoM[1] - gb->CoM[1]) / fabs(ga->CoM[1] + gb->CoM[1]) >
      tolerance) {
    message("CoM[1] different");
    return 0;
  }
  if (fabs(ga->CoM[2] - gb->CoM[2]) / fabs(ga->CoM[2] + gb->CoM[2]) >
      tolerance) {
    message("CoM[2] different");
    return 0;
  }

  /* Helper pointers */
  const struct multipole *ma = &ga->m_pole;
  const struct multipole *mb = &gb->m_pole;

  const double v2 = ma->vel[0] * ma->vel[0] + ma->vel[1] * ma->vel[1] +
                    ma->vel[2] * ma->vel[2];

  /* Check maximal softening */
  if (fabsf(ma->max_softening - mb->max_softening) /
          fabsf(ma->max_softening + mb->max_softening) >
      tolerance) {
    message("max softening different!");
    return 0;
  }

  /* Check bulk velocity (if non-zero and component > 1% of norm)*/
  if (fabsf(ma->vel[0] + mb->vel[0]) > 1e-10 &&
      (ma->vel[0] * ma->vel[0]) > 0.0001 * v2 &&
      fabsf(ma->vel[0] - mb->vel[0]) / fabsf(ma->vel[0] + mb->vel[0]) >
          tolerance) {
    message("v[0] different");
    return 0;
  }
  if (fabsf(ma->vel[1] + mb->vel[1]) > 1e-10 &&
      (ma->vel[1] * ma->vel[1]) > 0.0001 * v2 &&
      fabsf(ma->vel[1] - mb->vel[1]) / fabsf(ma->vel[1] + mb->vel[1]) >
          tolerance) {
    message("v[1] different");
    return 0;
  }
  if (fabsf(ma->vel[2] + mb->vel[2]) > 1e-10 &&
      (ma->vel[2] * ma->vel[2]) > 0.0001 * v2 &&
      fabsf(ma->vel[2] - mb->vel[2]) / fabsf(ma->vel[2] + mb->vel[2]) >
          tolerance) {
    message("v[2] different");
    return 0;
  }

  /* Manhattan Norm of 0th order terms */
  const float order0_norm = fabsf(ma->M_000) + fabsf(mb->M_000);

  /* Compare 0th order terms above 1% of norm */
  if (fabsf(ma->M_000 + mb->M_000) > 0.01f * order0_norm &&
      fabsf(ma->M_000 - mb->M_000) / fabsf(ma->M_000 + mb->M_000) > tolerance) {
    message("M_000 term different");
    return 0;
  }
#if SELF_GRAVITY_MULTIPOLE_ORDER > 0
  /* Manhattan Norm of 1st order terms */
  const float order1_norm = fabsf(ma->M_001) + fabsf(mb->M_001) +
                            fabsf(ma->M_010) + fabsf(mb->M_010) +
                            fabsf(ma->M_100) + fabsf(mb->M_100);

  /* Compare 1st order terms above 1% of norm */
  if (fabsf(ma->M_001 + mb->M_001) > 0.01f * order1_norm &&
      fabsf(ma->M_001 - mb->M_001) / fabsf(ma->M_001 + mb->M_001) > tolerance) {
    message("M_001 term different");
    return 0;
  }
  if (fabsf(ma->M_010 + mb->M_010) > 0.01f * order1_norm &&
      fabsf(ma->M_010 - mb->M_010) / fabsf(ma->M_010 + mb->M_010) > tolerance) {
    message("M_010 term different");
    return 0;
  }
  if (fabsf(ma->M_100 + mb->M_100) > 0.01f * order1_norm &&
      fabsf(ma->M_100 - mb->M_100) / fabsf(ma->M_100 + mb->M_100) > tolerance) {
    message("M_100 term different");
    return 0;
  }
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 1
  /* Manhattan Norm of 2nd order terms */
  const float order2_norm =
      fabsf(ma->M_002) + fabsf(mb->M_002) + fabsf(ma->M_011) +
      fabsf(mb->M_011) + fabsf(ma->M_020) + fabsf(mb->M_020) +
      fabsf(ma->M_101) + fabsf(mb->M_101) + fabsf(ma->M_110) +
      fabsf(mb->M_110) + fabsf(ma->M_200) + fabsf(mb->M_200);

  /* Compare 2nd order terms above 1% of norm */
  if (fabsf(ma->M_002 + mb->M_002) > 0.01f * order2_norm &&
      fabsf(ma->M_002 - mb->M_002) / fabsf(ma->M_002 + mb->M_002) > tolerance) {
    message("M_002 term different");
    return 0;
  }
  if (fabsf(ma->M_011 + mb->M_011) > 0.01f * order2_norm &&
      fabsf(ma->M_011 - mb->M_011) / fabsf(ma->M_011 + mb->M_011) > tolerance) {
    message("M_011 term different");
    return 0;
  }
  if (fabsf(ma->M_020 + mb->M_020) > 0.01f * order2_norm &&
      fabsf(ma->M_020 - mb->M_020) / fabsf(ma->M_020 + mb->M_020) > tolerance) {
    message("M_020 term different");
    return 0;
  }
  if (fabsf(ma->M_101 + mb->M_101) > 0.01f * order2_norm &&
      fabsf(ma->M_101 - mb->M_101) / fabsf(ma->M_101 + mb->M_101) > tolerance) {
    message("M_101 term different");
    return 0;
  }
  if (fabsf(ma->M_110 + mb->M_110) > 0.01f * order2_norm &&
      fabsf(ma->M_110 - mb->M_110) / fabsf(ma->M_110 + mb->M_110) > tolerance) {
    message("M_110 term different");
    return 0;
  }
  if (fabsf(ma->M_200 + mb->M_200) > 0.01f * order2_norm &&
      fabsf(ma->M_200 - mb->M_200) / fabsf(ma->M_200 + mb->M_200) > tolerance) {
    message("M_200 term different");
    return 0;
  }
#endif
  tolerance *= 10.;
#if SELF_GRAVITY_MULTIPOLE_ORDER > 2
  /* Manhattan Norm of 3rd order terms */
  const float order3_norm =
      fabsf(ma->M_003) + fabsf(mb->M_003) + fabsf(ma->M_012) +
      fabsf(mb->M_012) + fabsf(ma->M_021) + fabsf(mb->M_021) +
      fabsf(ma->M_030) + fabsf(mb->M_030) + fabsf(ma->M_102) +
      fabsf(mb->M_102) + fabsf(ma->M_111) + fabsf(mb->M_111) +
      fabsf(ma->M_120) + fabsf(mb->M_120) + fabsf(ma->M_201) +
      fabsf(mb->M_201) + fabsf(ma->M_210) + fabsf(mb->M_210) +
      fabsf(ma->M_300) + fabsf(mb->M_300);

  /* Compare 3rd order terms above 1% of norm */
  if (fabsf(ma->M_003 + mb->M_003) > 0.01f * order3_norm &&
      fabsf(ma->M_003 - mb->M_003) / fabsf(ma->M_003 + mb->M_003) > tolerance) {
    message("M_003 term different");
    return 0;
  }
  if (fabsf(ma->M_012 + mb->M_012) > 0.01f * order3_norm &&
      fabsf(ma->M_012 - mb->M_012) / fabsf(ma->M_012 + mb->M_012) > tolerance) {
    message("M_012 term different");
    return 0;
  }
  if (fabsf(ma->M_021 + mb->M_021) > 0.01f * order3_norm &&
      fabsf(ma->M_021 - mb->M_021) / fabsf(ma->M_021 + mb->M_021) > tolerance) {
    message("M_021 term different");
    return 0;
  }
  if (fabsf(ma->M_030 + mb->M_030) > 0.01f * order3_norm &&
      fabsf(ma->M_030 - mb->M_030) / fabsf(ma->M_030 + mb->M_030) > tolerance) {
    message("M_030 term different");
    return 0;
  }
  if (fabsf(ma->M_102 + mb->M_102) > 0.01f * order3_norm &&
      fabsf(ma->M_102 - mb->M_102) / fabsf(ma->M_102 + mb->M_102) > tolerance) {
    message("M_102 term different");
    return 0;
  }
  if (fabsf(ma->M_111 + mb->M_111) > 0.01f * order3_norm &&
      fabsf(ma->M_111 - mb->M_111) / fabsf(ma->M_111 + mb->M_111) > tolerance) {
    message("M_111 term different");
    return 0;
  }
  if (fabsf(ma->M_120 + mb->M_120) > 0.01f * order3_norm &&
      fabsf(ma->M_120 - mb->M_120) / fabsf(ma->M_120 + mb->M_120) > tolerance) {
    message("M_120 term different");
    return 0;
  }
  if (fabsf(ma->M_201 + mb->M_201) > 0.01f * order3_norm &&
      fabsf(ma->M_201 - mb->M_201) / fabsf(ma->M_201 + mb->M_201) > tolerance) {
    message("M_201 term different");
    return 0;
  }
  if (fabsf(ma->M_210 + mb->M_210) > 0.01f * order3_norm &&
      fabsf(ma->M_210 - mb->M_210) / fabsf(ma->M_210 + mb->M_210) > tolerance) {
    message("M_210 term different");
    return 0;
  }
  if (fabsf(ma->M_300 + mb->M_300) > 0.01f * order3_norm &&
      fabsf(ma->M_300 - mb->M_300) / fabsf(ma->M_300 + mb->M_300) > tolerance) {
    message("M_300 term different");
    return 0;
  }
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 3
  /* Manhattan Norm of 4th order terms */
  const float order4_norm =
      fabsf(ma->M_004) + fabsf(mb->M_004) + fabsf(ma->M_013) +
      fabsf(mb->M_013) + fabsf(ma->M_022) + fabsf(mb->M_022) +
      fabsf(ma->M_031) + fabsf(mb->M_031) + fabsf(ma->M_040) +
      fabsf(mb->M_040) + fabsf(ma->M_103) + fabsf(mb->M_103) +
      fabsf(ma->M_112) + fabsf(mb->M_112) + fabsf(ma->M_121) +
      fabsf(mb->M_121) + fabsf(ma->M_130) + fabsf(mb->M_130) +
      fabsf(ma->M_202) + fabsf(mb->M_202) + fabsf(ma->M_211) +
      fabsf(mb->M_211) + fabsf(ma->M_220) + fabsf(mb->M_220) +
      fabsf(ma->M_301) + fabsf(mb->M_301) + fabsf(ma->M_310) +
      fabsf(mb->M_310) + fabsf(ma->M_400) + fabsf(mb->M_400);

  /* Compare 4th order terms above 1% of norm */
  if (fabsf(ma->M_004 + mb->M_004) > 0.01f * order4_norm &&
      fabsf(ma->M_004 - mb->M_004) / fabsf(ma->M_004 + mb->M_004) > tolerance) {
    message("M_004 term different");
    return 0;
  }
  if (fabsf(ma->M_013 + mb->M_013) > 0.01f * order4_norm &&
      fabsf(ma->M_013 - mb->M_013) / fabsf(ma->M_013 + mb->M_013) > tolerance) {
    message("M_013 term different");
    return 0;
  }
  if (fabsf(ma->M_022 + mb->M_022) > 0.01f * order4_norm &&
      fabsf(ma->M_022 - mb->M_022) / fabsf(ma->M_022 + mb->M_022) > tolerance) {
    message("M_022 term different");
    return 0;
  }
  if (fabsf(ma->M_031 + mb->M_031) > 0.01f * order4_norm &&
      fabsf(ma->M_031 - mb->M_031) / fabsf(ma->M_031 + mb->M_031) > tolerance) {
    message("M_031 term different");
    return 0;
  }
  if (fabsf(ma->M_040 + mb->M_040) > 0.01f * order4_norm &&
      fabsf(ma->M_040 - mb->M_040) / fabsf(ma->M_040 + mb->M_040) > tolerance) {
    message("M_040 term different");
    return 0;
  }
  if (fabsf(ma->M_103 + mb->M_103) > 0.01f * order4_norm &&
      fabsf(ma->M_103 - mb->M_103) / fabsf(ma->M_103 + mb->M_103) > tolerance) {
    message("M_103 term different");
    return 0;
  }
  if (fabsf(ma->M_112 + mb->M_112) > 0.01f * order4_norm &&
      fabsf(ma->M_112 - mb->M_112) / fabsf(ma->M_112 + mb->M_112) > tolerance) {
    message("M_112 term different");
    return 0;
  }
  if (fabsf(ma->M_121 + mb->M_121) > 0.01f * order4_norm &&
      fabsf(ma->M_121 - mb->M_121) / fabsf(ma->M_121 + mb->M_121) > tolerance) {
    message("M_121 term different");
    return 0;
  }
  if (fabsf(ma->M_130 + mb->M_130) > 0.01f * order4_norm &&
      fabsf(ma->M_130 - mb->M_130) / fabsf(ma->M_130 + mb->M_130) > tolerance) {
    message("M_130 term different");
    return 0;
  }
  if (fabsf(ma->M_202 + mb->M_202) > 0.01f * order4_norm &&
      fabsf(ma->M_202 - mb->M_202) / fabsf(ma->M_202 + mb->M_202) > tolerance) {
    message("M_202 term different");
    return 0;
  }
  if (fabsf(ma->M_211 + mb->M_211) > 0.01f * order4_norm &&
      fabsf(ma->M_211 - mb->M_211) / fabsf(ma->M_211 + mb->M_211) > tolerance) {
    message("M_211 term different");
    return 0;
  }
  if (fabsf(ma->M_220 + mb->M_220) > 0.01f * order4_norm &&
      fabsf(ma->M_220 - mb->M_220) / fabsf(ma->M_220 + mb->M_220) > tolerance) {
    message("M_220 term different");
    return 0;
  }
  if (fabsf(ma->M_301 + mb->M_301) > 0.01f * order4_norm &&
      fabsf(ma->M_301 - mb->M_301) / fabsf(ma->M_301 + mb->M_301) > tolerance) {
    message("M_301 term different");
    return 0;
  }
  if (fabsf(ma->M_310 + mb->M_310) > 0.01f * order4_norm &&
      fabsf(ma->M_310 - mb->M_310) / fabsf(ma->M_310 + mb->M_310) > tolerance) {
    message("M_310 term different");
    return 0;
  }
  if (fabsf(ma->M_400 + mb->M_400) > 0.01f * order4_norm &&
      fabsf(ma->M_400 - mb->M_400) / fabsf(ma->M_400 + mb->M_400) > tolerance) {
    message("M_400 term different");
    return 0;
  }
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 4
  /* Manhattan Norm of 5th order terms */
  const float order5_norm =
      fabsf(ma->M_005) + fabsf(mb->M_005) + fabsf(ma->M_014) +
      fabsf(mb->M_014) + fabsf(ma->M_023) + fabsf(mb->M_023) +
      fabsf(ma->M_032) + fabsf(mb->M_032) + fabsf(ma->M_041) +
      fabsf(mb->M_041) + fabsf(ma->M_050) + fabsf(mb->M_050) +
      fabsf(ma->M_104) + fabsf(mb->M_104) + fabsf(ma->M_113) +
      fabsf(mb->M_113) + fabsf(ma->M_122) + fabsf(mb->M_122) +
      fabsf(ma->M_131) + fabsf(mb->M_131) + fabsf(ma->M_140) +
      fabsf(mb->M_140) + fabsf(ma->M_203) + fabsf(mb->M_203) +
      fabsf(ma->M_212) + fabsf(mb->M_212) + fabsf(ma->M_221) +
      fabsf(mb->M_221) + fabsf(ma->M_230) + fabsf(mb->M_230) +
      fabsf(ma->M_302) + fabsf(mb->M_302) + fabsf(ma->M_311) +
      fabsf(mb->M_311) + fabsf(ma->M_320) + fabsf(mb->M_320) +
      fabsf(ma->M_401) + fabsf(mb->M_401) + fabsf(ma->M_410) +
      fabsf(mb->M_410) + fabsf(ma->M_500) + fabsf(mb->M_500);

  /* Compare 5th order terms above 1% of norm */
  if (fabsf(ma->M_005 + mb->M_005) > 0.01f * order5_norm &&
      fabsf(ma->M_005 - mb->M_005) / fabsf(ma->M_005 + mb->M_005) > tolerance) {
    message("M_005 term different");
    return 0;
  }
  if (fabsf(ma->M_014 + mb->M_014) > 0.01f * order5_norm &&
      fabsf(ma->M_014 - mb->M_014) / fabsf(ma->M_014 + mb->M_014) > tolerance) {
    message("M_014 term different");
    return 0;
  }
  if (fabsf(ma->M_023 + mb->M_023) > 0.01f * order5_norm &&
      fabsf(ma->M_023 - mb->M_023) / fabsf(ma->M_023 + mb->M_023) > tolerance) {
    message("M_023 term different");
    return 0;
  }
  if (fabsf(ma->M_032 + mb->M_032) > 0.01f * order5_norm &&
      fabsf(ma->M_032 - mb->M_032) / fabsf(ma->M_032 + mb->M_032) > tolerance) {
    message("M_032 term different");
    return 0;
  }
  if (fabsf(ma->M_041 + mb->M_041) > 0.01f * order5_norm &&
      fabsf(ma->M_041 - mb->M_041) / fabsf(ma->M_041 + mb->M_041) > tolerance) {
    message("M_041 term different");
    return 0;
  }
  if (fabsf(ma->M_050 + mb->M_050) > 0.01f * order5_norm &&
      fabsf(ma->M_050 - mb->M_050) / fabsf(ma->M_050 + mb->M_050) > tolerance) {
    message("M_050 term different");
    return 0;
  }
  if (fabsf(ma->M_104 + mb->M_104) > 0.01f * order5_norm &&
      fabsf(ma->M_104 - mb->M_104) / fabsf(ma->M_104 + mb->M_104) > tolerance) {
    message("M_104 term different");
    return 0;
  }
  if (fabsf(ma->M_113 + mb->M_113) > 0.01f * order5_norm &&
      fabsf(ma->M_113 - mb->M_113) / fabsf(ma->M_113 + mb->M_113) > tolerance) {
    message("M_113 term different");
    return 0;
  }
  if (fabsf(ma->M_122 + mb->M_122) > 0.01f * order5_norm &&
      fabsf(ma->M_122 - mb->M_122) / fabsf(ma->M_122 + mb->M_122) > tolerance) {
    message("M_122 term different");
    return 0;
  }
  if (fabsf(ma->M_131 + mb->M_131) > 0.01f * order5_norm &&
      fabsf(ma->M_131 - mb->M_131) / fabsf(ma->M_131 + mb->M_131) > tolerance) {
    message("M_131 term different");
    return 0;
  }
  if (fabsf(ma->M_140 + mb->M_140) > 0.01f * order5_norm &&
      fabsf(ma->M_140 - mb->M_140) / fabsf(ma->M_140 + mb->M_140) > tolerance) {
    message("M_140 term different");
    return 0;
  }
  if (fabsf(ma->M_203 + mb->M_203) > 0.01f * order5_norm &&
      fabsf(ma->M_203 - mb->M_203) / fabsf(ma->M_203 + mb->M_203) > tolerance) {
    message("M_203 term different");
    return 0;
  }
  if (fabsf(ma->M_212 + mb->M_212) > 0.01f * order5_norm &&
      fabsf(ma->M_212 - mb->M_212) / fabsf(ma->M_212 + mb->M_212) > tolerance) {
    message("M_212 term different");
    return 0;
  }
  if (fabsf(ma->M_221 + mb->M_221) > 0.01f * order5_norm &&
      fabsf(ma->M_221 - mb->M_221) / fabsf(ma->M_221 + mb->M_221) > tolerance) {
    message("M_221 term different");
    return 0;
  }
  if (fabsf(ma->M_230 + mb->M_230) > 0.01f * order5_norm &&
      fabsf(ma->M_230 - mb->M_230) / fabsf(ma->M_230 + mb->M_230) > tolerance) {
    message("M_230 term different");
    return 0;
  }
  if (fabsf(ma->M_302 + mb->M_302) > 0.01f * order5_norm &&
      fabsf(ma->M_302 - mb->M_302) / fabsf(ma->M_302 + mb->M_302) > tolerance) {
    message("M_302 term different");
    return 0;
  }
  if (fabsf(ma->M_311 + mb->M_311) > 0.01f * order5_norm &&
      fabsf(ma->M_311 - mb->M_311) / fabsf(ma->M_311 + mb->M_311) > tolerance) {
    message("M_311 term different");
    return 0;
  }
  if (fabsf(ma->M_320 + mb->M_320) > 0.01f * order5_norm &&
      fabsf(ma->M_320 - mb->M_320) / fabsf(ma->M_320 + mb->M_320) > tolerance) {
    message("M_320 term different");
    return 0;
  }
  if (fabsf(ma->M_401 + mb->M_401) > 0.01f * order5_norm &&
      fabsf(ma->M_401 - mb->M_401) / fabsf(ma->M_401 + mb->M_401) > tolerance) {
    message("M_401 term different");
    return 0;
  }
  if (fabsf(ma->M_410 + mb->M_410) > 0.01f * order5_norm &&
      fabsf(ma->M_410 - mb->M_410) / fabsf(ma->M_410 + mb->M_410) > tolerance) {
    message("M_410 term different");
    return 0;
  }
  if (fabsf(ma->M_500 + mb->M_500) > 0.01f * order5_norm &&
      fabsf(ma->M_500 - mb->M_500) / fabsf(ma->M_500 + mb->M_500) > tolerance) {
    message("M_500 term different");
    return 0;
  }
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 5
  /* Manhattan Norm of 6th order terms */
  const float order6_norm =
      fabsf(ma->M_006) + fabsf(mb->M_006)
    + fabsf(ma->M_015) + fabsf(mb->M_015)
    + fabsf(ma->M_024) + fabsf(mb->M_024)
    + fabsf(ma->M_033) + fabsf(mb->M_033)
    + fabsf(ma->M_042) + fabsf(mb->M_042)
    + fabsf(ma->M_051) + fabsf(mb->M_051)
    + fabsf(ma->M_060) + fabsf(mb->M_060)
    + fabsf(ma->M_105) + fabsf(mb->M_105)
    + fabsf(ma->M_114) + fabsf(mb->M_114)
    + fabsf(ma->M_123) + fabsf(mb->M_123)
    + fabsf(ma->M_132) + fabsf(mb->M_132)
    + fabsf(ma->M_141) + fabsf(mb->M_141)
    + fabsf(ma->M_150) + fabsf(mb->M_150)
    + fabsf(ma->M_204) + fabsf(mb->M_204)
    + fabsf(ma->M_213) + fabsf(mb->M_213)
    + fabsf(ma->M_222) + fabsf(mb->M_222)
    + fabsf(ma->M_231) + fabsf(mb->M_231)
    + fabsf(ma->M_240) + fabsf(mb->M_240)
    + fabsf(ma->M_303) + fabsf(mb->M_303)
    + fabsf(ma->M_312) + fabsf(mb->M_312)
    + fabsf(ma->M_321) + fabsf(mb->M_321)
    + fabsf(ma->M_330) + fabsf(mb->M_330)
    + fabsf(ma->M_402) + fabsf(mb->M_402)
    + fabsf(ma->M_411) + fabsf(mb->M_411)
    + fabsf(ma->M_420) + fabsf(mb->M_420)
    + fabsf(ma->M_501) + fabsf(mb->M_501)
    + fabsf(ma->M_510) + fabsf(mb->M_510)
    + fabsf(ma->M_600) + fabsf(mb->M_600);

  /* Compare 6th order terms above 1% of norm */
  if (fabsf(ma->M_006 + mb->M_006) > 0.01f * order6_norm &&
      fabsf(ma->M_006 - mb->M_006) / fabsf(ma->M_006 + mb->M_006) > tolerance) {
    message("M_006 term different");
    return 0;
  }
  if (fabsf(ma->M_015 + mb->M_015) > 0.01f * order6_norm &&
      fabsf(ma->M_015 - mb->M_015) / fabsf(ma->M_015 + mb->M_015) > tolerance) {
    message("M_015 term different");
    return 0;
  }
  if (fabsf(ma->M_024 + mb->M_024) > 0.01f * order6_norm &&
      fabsf(ma->M_024 - mb->M_024) / fabsf(ma->M_024 + mb->M_024) > tolerance) {
    message("M_024 term different");
    return 0;
  }
  if (fabsf(ma->M_033 + mb->M_033) > 0.01f * order6_norm &&
      fabsf(ma->M_033 - mb->M_033) / fabsf(ma->M_033 + mb->M_033) > tolerance) {
    message("M_033 term different");
    return 0;
  }
  if (fabsf(ma->M_042 + mb->M_042) > 0.01f * order6_norm &&
      fabsf(ma->M_042 - mb->M_042) / fabsf(ma->M_042 + mb->M_042) > tolerance) {
    message("M_042 term different");
    return 0;
  }
  if (fabsf(ma->M_051 + mb->M_051) > 0.01f * order6_norm &&
      fabsf(ma->M_051 - mb->M_051) / fabsf(ma->M_051 + mb->M_051) > tolerance) {
    message("M_051 term different");
    return 0;
  }
  if (fabsf(ma->M_060 + mb->M_060) > 0.01f * order6_norm &&
      fabsf(ma->M_060 - mb->M_060) / fabsf(ma->M_060 + mb->M_060) > tolerance) {
    message("M_060 term different");
    return 0;
  }
  if (fabsf(ma->M_105 + mb->M_105) > 0.01f * order6_norm &&
      fabsf(ma->M_105 - mb->M_105) / fabsf(ma->M_105 + mb->M_105) > tolerance) {
    message("M_105 term different");
    return 0;
  }
  if (fabsf(ma->M_114 + mb->M_114) > 0.01f * order6_norm &&
      fabsf(ma->M_114 - mb->M_114) / fabsf(ma->M_114 + mb->M_114) > tolerance) {
    message("M_114 term different");
    return 0;
  }
  if (fabsf(ma->M_123 + mb->M_123) > 0.01f * order6_norm &&
      fabsf(ma->M_123 - mb->M_123) / fabsf(ma->M_123 + mb->M_123) > tolerance) {
    message("M_123 term different");
    return 0;
  }
  if (fabsf(ma->M_132 + mb->M_132) > 0.01f * order6_norm &&
      fabsf(ma->M_132 - mb->M_132) / fabsf(ma->M_132 + mb->M_132) > tolerance) {
    message("M_132 term different");
    return 0;
  }
  if (fabsf(ma->M_141 + mb->M_141) > 0.01f * order6_norm &&
      fabsf(ma->M_141 - mb->M_141) / fabsf(ma->M_141 + mb->M_141) > tolerance) {
    message("M_141 term different");
    return 0;
  }
  if (fabsf(ma->M_150 + mb->M_150) > 0.01f * order6_norm &&
      fabsf(ma->M_150 - mb->M_150) / fabsf(ma->M_150 + mb->M_150) > tolerance) {
    message("M_150 term different");
    return 0;
  }
  if (fabsf(ma->M_204 + mb->M_204) > 0.01f * order6_norm &&
      fabsf(ma->M_204 - mb->M_204) / fabsf(ma->M_204 + mb->M_204) > tolerance) {
    message("M_204 term different");
    return 0;
  }
  if (fabsf(ma->M_213 + mb->M_213) > 0.01f * order6_norm &&
      fabsf(ma->M_213 - mb->M_213) / fabsf(ma->M_213 + mb->M_213) > tolerance) {
    message("M_213 term different");
    return 0;
  }
  if (fabsf(ma->M_222 + mb->M_222) > 0.01f * order6_norm &&
      fabsf(ma->M_222 - mb->M_222) / fabsf(ma->M_222 + mb->M_222) > tolerance) {
    message("M_222 term different");
    return 0;
  }
  if (fabsf(ma->M_231 + mb->M_231) > 0.01f * order6_norm &&
      fabsf(ma->M_231 - mb->M_231) / fabsf(ma->M_231 + mb->M_231) > tolerance) {
    message("M_231 term different");
    return 0;
  }
  if (fabsf(ma->M_240 + mb->M_240) > 0.01f * order6_norm &&
      fabsf(ma->M_240 - mb->M_240) / fabsf(ma->M_240 + mb->M_240) > tolerance) {
    message("M_240 term different");
    return 0;
  }
  if (fabsf(ma->M_303 + mb->M_303) > 0.01f * order6_norm &&
      fabsf(ma->M_303 - mb->M_303) / fabsf(ma->M_303 + mb->M_303) > tolerance) {
    message("M_303 term different");
    return 0;
  }
  if (fabsf(ma->M_312 + mb->M_312) > 0.01f * order6_norm &&
      fabsf(ma->M_312 - mb->M_312) / fabsf(ma->M_312 + mb->M_312) > tolerance) {
    message("M_312 term different");
    return 0;
  }
  if (fabsf(ma->M_321 + mb->M_321) > 0.01f * order6_norm &&
      fabsf(ma->M_321 - mb->M_321) / fabsf(ma->M_321 + mb->M_321) > tolerance) {
    message("M_321 term different");
    return 0;
  }
  if (fabsf(ma->M_330 + mb->M_330) > 0.01f * order6_norm &&
      fabsf(ma->M_330 - mb->M_330) / fabsf(ma->M_330 + mb->M_330) > tolerance) {
    message("M_330 term different");
    return 0;
  }
  if (fabsf(ma->M_402 + mb->M_402) > 0.01f * order6_norm &&
      fabsf(ma->M_402 - mb->M_402) / fabsf(ma->M_402 + mb->M_402) > tolerance) {
    message("M_402 term different");
    return 0;
  }
  if (fabsf(ma->M_411 + mb->M_411) > 0.01f * order6_norm &&
      fabsf(ma->M_411 - mb->M_411) / fabsf(ma->M_411 + mb->M_411) > tolerance) {
    message("M_411 term different");
    return 0;
  }
  if (fabsf(ma->M_420 + mb->M_420) > 0.01f * order6_norm &&
      fabsf(ma->M_420 - mb->M_420) / fabsf(ma->M_420 + mb->M_420) > tolerance) {
    message("M_420 term different");
    return 0;
  }
  if (fabsf(ma->M_501 + mb->M_501) > 0.01f * order6_norm &&
      fabsf(ma->M_501 - mb->M_501) / fabsf(ma->M_501 + mb->M_501) > tolerance) {
    message("M_501 term different");
    return 0;
  }
  if (fabsf(ma->M_510 + mb->M_510) > 0.01f * order6_norm &&
      fabsf(ma->M_510 - mb->M_510) / fabsf(ma->M_510 + mb->M_510) > tolerance) {
    message("M_510 term different");
    return 0;
  }
  if (fabsf(ma->M_600 + mb->M_600) > 0.01f * order6_norm &&
      fabsf(ma->M_600 - mb->M_600) / fabsf(ma->M_600 + mb->M_600) > tolerance) {
    message("M_600 term different");
    return 0;
  }
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 6
#error "Missing implementation for order >6"
#endif

  /* All is good */
  return 1;
}

/**
 * @brief XXX 
 * 
 *
 * Corresponds to equation (XXX).
 *
 * @param multi The #multipole we are computing the power for.
 */
INLINE static void compute_multipole_power(struct multipole *multi) {
#ifdef ADVANCED_OPENING_CRITERIA

  /* 0th order contribution */
  multi->power[0] = multi->M_000 * multi->M_000;

#if SELF_GRAVITY_MULTIPOLE_ORDER > 0
  /* 1st order contributions */
  multi->power[1] = 0.0;
  multi->power[1] += multi->M_100 * multi->M_100;
  multi->power[1] += multi->M_010 * multi->M_010;
  multi->power[1] += multi->M_001 * multi->M_001;
#endif

#if SELF_GRAVITY_MULTIPOLE_ORDER > 1
  /* 2nd order contributions */
  multi->power[2] = 0.0;
  multi->power[2] += multi->M_200 * multi->M_200;
  multi->power[2] += multi->M_020 * multi->M_020;
  multi->power[2] += multi->M_002 * multi->M_002;

  multi->power[2] += 0.5 * (multi->M_110 * multi->M_110);
  multi->power[2] += 0.5 * (multi->M_101 * multi->M_101);
  multi->power[2] += 0.5 * (multi->M_011 * multi->M_011);
#endif

#if SELF_GRAVITY_MULTIPOLE_ORDER > 2
  /* 3nd order contributions */
  multi->power[3] = 0.0;
  multi->power[3] += 2.0 * (multi->M_300 * multi->M_300);
  multi->power[3] += 2.0 * (multi->M_030 * multi->M_030);
  multi->power[3] += 2.0 * (multi->M_003 * multi->M_003);

  multi->power[3] += 0.6666666666666667 * (multi->M_210 * multi->M_210);
  multi->power[3] += 0.6666666666666667 * (multi->M_201 * multi->M_201);
  multi->power[3] += 0.6666666666666667 * (multi->M_120 * multi->M_120);
  multi->power[3] += 0.6666666666666667 * (multi->M_021 * multi->M_021);
  multi->power[3] += 0.6666666666666667 * (multi->M_102 * multi->M_102);
  multi->power[3] += 0.6666666666666667 * (multi->M_021 * multi->M_021);

  multi->power[3] += 0.3333333333333333 * (multi->M_111 * multi->M_111);
#endif

#if SELF_GRAVITY_MULTIPOLE_ORDER > 3
  /*4th order contributions */
  multi->power[4] = 0.0;
  multi->power[4] += 6.0 * (multi->M_400 * multi->M_400);
  multi->power[4] += 6.0 * (multi->M_040 * multi->M_040);
  multi->power[4] += 6.0 * (multi->M_004 * multi->M_004);

  multi->power[4] += 1.5 * (multi->M_310 * multi->M_310);
  multi->power[4] += 1.5 * (multi->M_301 * multi->M_301);
  multi->power[4] += 1.5 * (multi->M_130 * multi->M_130);
  multi->power[4] += 1.5 * (multi->M_031 * multi->M_031);
  multi->power[4] += 1.5 * (multi->M_103 * multi->M_103);
  multi->power[4] += 1.5 * (multi->M_013 * multi->M_013);

  multi->power[4] += multi->M_220 * multi->M_220;
  multi->power[4] += multi->M_202 * multi->M_202;
  multi->power[4] += multi->M_022 * multi->M_022;

  multi->power[4] += 0.5 * (multi->M_211 * multi->M_211);
  multi->power[4] += 0.5 * (multi->M_121 * multi->M_121);
  multi->power[4] == 0.5 * (multi->M_112 * multi->M_112);
#endif

#if SELF_GRAVITY_MULTIPOLE_ORDER > 4
  /*5th order contributions */
  multi->power[5] = 0.0;
  multi->power[5] += 24. * (multi->M_005 * multi->M_005);
  multi->power[5] += 24. * (multi->M_050 * multi->M_050);
  multi->power[5] += 24. * (multi->M_500 * multi->M_500);

  multi->power[5] += 4.8 * (multi->M_014 * multi->M_014);
  multi->power[5] += 4.8 * (multi->M_104 * multi->M_104);
  multi->power[5] += 4.8 * (multi->M_041 * multi->M_041);
  multi->power[5] += 4.8 * (multi->M_140 * multi->M_140);
  multi->power[5] += 4.8 * (multi->M_401 * multi->M_401);
  multi->power[5] += 4.8 * (multi->M_410 * multi->M_410);

  multi->power[5] += 2.4 * (multi->M_032 * multi->M_032);
  multi->power[5] += 2.4 * (multi->M_023 * multi->M_023);
  multi->power[5] += 2.4 * (multi->M_203 * multi->M_203);
  multi->power[5] += 2.4 * (multi->M_230 * multi->M_230);
  multi->power[5] += 2.4 * (multi->M_302 * multi->M_302);
  multi->power[5] += 2.4 * (multi->M_320 * multi->M_320);

  multi->power[5] += 0.8 * (multi->M_122 * multi->M_122);
  multi->power[5] += 0.8 * (multi->M_221 * multi->M_221);
  multi->power[5] += 0.8 * (multi->M_212 * multi->M_212);

  multi->power[5] += 1.2 * (multi->M_113 * multi->M_113);
  multi->power[5] += 1.2 * (multi->M_131 * multi->M_131);
  multi->power[5] += 1.2 * (multi->M_311 * multi->M_311);
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 5
#error "Missing implementation for order >5"
#endif
#endif
}

/**
 * @brief Constructs the #multipole of a bunch of particles around their
 * centre of mass.
 *
 * Corresponds to equation (28c).
 *
 * @param multi The #multipole (content will  be overwritten).
 * @param gparts The #gpart.
 * @param gcount The number of particles.
 * @param grav_props The properties of the gravity scheme.
 */
INLINE static void gravity_P2M(struct gravity_tensors *multi,
                               const struct gpart *gparts, const int gcount,
                               const struct gravity_props *const grav_props) {

  /* Temporary variables */
  float epsilon_max = 0.f;
  double mass = 0.0;
  double com[3] = {0.0, 0.0, 0.0};
  double vel[3] = {0.f, 0.f, 0.f};
#ifdef ADVANCED_OPENING_CRITERIA
  float min_a_grav_norm = FLT_MAX;
#endif

  /* Collect the particle data for CoM. */
  for (int k = 0; k < gcount; k++) {
    const double m = gparts[k].mass;
    const float epsilon = gravity_get_softening(&gparts[k], grav_props);

#ifdef SWIFT_DEBUG_CHECKS
    if (gparts[k].time_bin == time_bin_inhibited)
      error("Inhibited particle in P2M. Should have been removed earlier.");
#endif

    epsilon_max = max(epsilon_max, epsilon);
    mass += m;
    com[0] += gparts[k].x[0] * m;
    com[1] += gparts[k].x[1] * m;
    com[2] += gparts[k].x[2] * m;
    vel[0] += gparts[k].v_full[0] * m;
    vel[1] += gparts[k].v_full[1] * m;
    vel[2] += gparts[k].v_full[2] * m;

#ifdef ADVANCED_OPENING_CRITERIA
    min_a_grav_norm = min(min_a_grav_norm, gparts[k].a_grav_norm);
#endif
  }

  /* Final operation on CoM */
  const double imass = 1.0 / mass;
  com[0] *= imass;
  com[1] *= imass;
  com[2] *= imass;
  vel[0] *= imass;
  vel[1] *= imass;
  vel[2] *= imass;

  /* Prepare some local counters */
  double r_max2 = 0.;
  float max_delta_vel[3] = {0., 0., 0.};
  float min_delta_vel[3] = {0., 0., 0.};

#if SELF_GRAVITY_MULTIPOLE_ORDER > 0
  double M_100 = 0., M_010 = 0., M_001 = 0.;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 1
  double M_200 = 0., M_020 = 0., M_002 = 0.;
  double M_110 = 0., M_101 = 0., M_011 = 0.;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 2
  double M_300 = 0., M_030 = 0., M_003 = 0.;
  double M_210 = 0., M_201 = 0., M_120 = 0.;
  double M_021 = 0., M_102 = 0., M_012 = 0.;
  double M_111 = 0.;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 3
  double M_400 = 0., M_040 = 0., M_004 = 0.;
  double M_310 = 0., M_301 = 0., M_130 = 0.;
  double M_031 = 0., M_103 = 0., M_013 = 0.;
  double M_220 = 0., M_202 = 0., M_022 = 0.;
  double M_211 = 0., M_121 = 0., M_112 = 0.;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 4
  double M_005 = 0., M_014 = 0., M_023 = 0.;
  double M_032 = 0., M_041 = 0., M_050 = 0.;
  double M_104 = 0., M_113 = 0., M_122 = 0.;
  double M_131 = 0., M_140 = 0., M_203 = 0.;
  double M_212 = 0., M_221 = 0., M_230 = 0.;
  double M_302 = 0., M_311 = 0., M_320 = 0.;
  double M_401 = 0., M_410 = 0., M_500 = 0.;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 5
  double M_006 = 0., M_015 = 0., M_024 = 0.;
  double M_033 = 0., M_042 = 0., M_051 = 0.;
  double M_060 = 0., M_105 = 0., M_114 = 0.;
  double M_123 = 0., M_132 = 0., M_141 = 0.;
  double M_150 = 0., M_204 = 0., M_213 = 0.;
  double M_222 = 0., M_231 = 0., M_240 = 0.;
  double M_303 = 0., M_312 = 0., M_321 = 0.;
  double M_330 = 0., M_402 = 0., M_411 = 0.;
  double M_420 = 0., M_501 = 0., M_510 = 0.;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 6
#error "Missing implementation for order >6"
#endif

  /* Construce the higher order terms */
  for (int k = 0; k < gcount; k++) {

    const double dx[3] = {gparts[k].x[0] - com[0], gparts[k].x[1] - com[1],
                          gparts[k].x[2] - com[2]};

    /* Maximal distance CoM<->gpart */
    r_max2 = max(r_max2, dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);

    /* Store the vector of the maximal vel difference */
    max_delta_vel[0] = max(gparts[k].v_full[0], max_delta_vel[0]);
    max_delta_vel[1] = max(gparts[k].v_full[1], max_delta_vel[1]);
    max_delta_vel[2] = max(gparts[k].v_full[2], max_delta_vel[2]);

    /* Store the vector of the minimal vel difference */
    min_delta_vel[0] = min(gparts[k].v_full[0], min_delta_vel[0]);
    min_delta_vel[1] = min(gparts[k].v_full[1], min_delta_vel[1]);
    min_delta_vel[2] = min(gparts[k].v_full[2], min_delta_vel[2]);

#if SELF_GRAVITY_MULTIPOLE_ORDER > 0
    const double m = gparts[k].mass;

    /* 1st order terms */
    M_100 += -m * X_100(dx);
    M_010 += -m * X_010(dx);
    M_001 += -m * X_001(dx);
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 1

    /* 2nd order terms */
    M_200 += m * X_200(dx);
    M_020 += m * X_020(dx);
    M_002 += m * X_002(dx);
    M_110 += m * X_110(dx);
    M_101 += m * X_101(dx);
    M_011 += m * X_011(dx);
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 2

    /* 3rd order terms */
    M_300 += -m * X_300(dx);
    M_030 += -m * X_030(dx);
    M_003 += -m * X_003(dx);
    M_210 += -m * X_210(dx);
    M_201 += -m * X_201(dx);
    M_120 += -m * X_120(dx);
    M_021 += -m * X_021(dx);
    M_102 += -m * X_102(dx);
    M_012 += -m * X_012(dx);
    M_111 += -m * X_111(dx);
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 3

    /* 4th order terms */
    M_400 += m * X_400(dx);
    M_040 += m * X_040(dx);
    M_004 += m * X_004(dx);
    M_310 += m * X_310(dx);
    M_301 += m * X_301(dx);
    M_130 += m * X_130(dx);
    M_031 += m * X_031(dx);
    M_103 += m * X_103(dx);
    M_013 += m * X_013(dx);
    M_220 += m * X_220(dx);
    M_202 += m * X_202(dx);
    M_022 += m * X_022(dx);
    M_211 += m * X_211(dx);
    M_121 += m * X_121(dx);
    M_112 += m * X_112(dx);
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 4

    /* 5th order terms */
    M_005 += -m * X_005(dx);
    M_014 += -m * X_014(dx);
    M_023 += -m * X_023(dx);
    M_032 += -m * X_032(dx);
    M_041 += -m * X_041(dx);
    M_050 += -m * X_050(dx);
    M_104 += -m * X_104(dx);
    M_113 += -m * X_113(dx);
    M_122 += -m * X_122(dx);
    M_131 += -m * X_131(dx);
    M_140 += -m * X_140(dx);
    M_203 += -m * X_203(dx);
    M_212 += -m * X_212(dx);
    M_221 += -m * X_221(dx);
    M_230 += -m * X_230(dx);
    M_302 += -m * X_302(dx);
    M_311 += -m * X_311(dx);
    M_320 += -m * X_320(dx);
    M_401 += -m * X_401(dx);
    M_410 += -m * X_410(dx);
    M_500 += -m * X_500(dx);
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 5
  /* 6th order terms */
    M_006 += m * X_006(dx);
    M_015 += m * X_015(dx);
    M_024 += m * X_024(dx);
    M_033 += m * X_033(dx);
    M_042 += m * X_042(dx);
    M_051 += m * X_051(dx);
    M_060 += m * X_060(dx);
    M_105 += m * X_105(dx);
    M_114 += m * X_114(dx);
    M_123 += m * X_123(dx);
    M_132 += m * X_132(dx);
    M_141 += m * X_141(dx);
    M_150 += m * X_150(dx);
    M_204 += m * X_204(dx);
    M_213 += m * X_213(dx);
    M_222 += m * X_222(dx);
    M_231 += m * X_231(dx);
    M_240 += m * X_240(dx);
    M_303 += m * X_303(dx);
    M_312 += m * X_312(dx);
    M_321 += m * X_321(dx);
    M_330 += m * X_330(dx);
    M_402 += m * X_402(dx);
    M_411 += m * X_411(dx);
    M_420 += m * X_420(dx);
    M_501 += m * X_501(dx);
    M_510 += m * X_510(dx);
    M_600 += m * X_600(dx);
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 6
#error "Missing implementation for order >6"
#endif
  }

#if SELF_GRAVITY_MULTIPOLE_ORDER > 0

  /* We know the first-order multipole (dipole) is 0. */
  M_100 = M_010 = M_001 = 0.f;
#endif

  /* Store the data on the multipole. */
  multi->m_pole.max_softening = epsilon_max;
  multi->m_pole.M_000 = mass;
  multi->r_max = sqrt(r_max2);
  multi->CoM[0] = com[0];
  multi->CoM[1] = com[1];
  multi->CoM[2] = com[2];
  multi->m_pole.vel[0] = vel[0];
  multi->m_pole.vel[1] = vel[1];
  multi->m_pole.vel[2] = vel[2];
  multi->m_pole.max_delta_vel[0] = max_delta_vel[0];
  multi->m_pole.max_delta_vel[1] = max_delta_vel[1];
  multi->m_pole.max_delta_vel[2] = max_delta_vel[2];
  multi->m_pole.min_delta_vel[0] = min_delta_vel[0];
  multi->m_pole.min_delta_vel[1] = min_delta_vel[1];
  multi->m_pole.min_delta_vel[2] = min_delta_vel[2];
#ifdef ADVANCED_OPENING_CRITERIA
  multi->m_pole.min_a_grav_norm = sqrt(min_a_grav_norm);
#endif

#if SELF_GRAVITY_MULTIPOLE_ORDER > 0

  /* 1st order terms */
  multi->m_pole.M_100 = M_100;
  multi->m_pole.M_010 = M_010;
  multi->m_pole.M_001 = M_001;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 1

  /* 2nd order terms */
  multi->m_pole.M_200 = M_200;
  multi->m_pole.M_020 = M_020;
  multi->m_pole.M_002 = M_002;
  multi->m_pole.M_110 = M_110;
  multi->m_pole.M_101 = M_101;
  multi->m_pole.M_011 = M_011;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 2

  /* 3rd order terms */
  multi->m_pole.M_300 = M_300;
  multi->m_pole.M_030 = M_030;
  multi->m_pole.M_003 = M_003;
  multi->m_pole.M_210 = M_210;
  multi->m_pole.M_201 = M_201;
  multi->m_pole.M_120 = M_120;
  multi->m_pole.M_021 = M_021;
  multi->m_pole.M_102 = M_102;
  multi->m_pole.M_012 = M_012;
  multi->m_pole.M_111 = M_111;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 3

  /* 4th order terms */
  multi->m_pole.M_400 = M_400;
  multi->m_pole.M_040 = M_040;
  multi->m_pole.M_004 = M_004;
  multi->m_pole.M_310 = M_310;
  multi->m_pole.M_301 = M_301;
  multi->m_pole.M_130 = M_130;
  multi->m_pole.M_031 = M_031;
  multi->m_pole.M_103 = M_103;
  multi->m_pole.M_013 = M_013;
  multi->m_pole.M_220 = M_220;
  multi->m_pole.M_202 = M_202;
  multi->m_pole.M_022 = M_022;
  multi->m_pole.M_211 = M_211;
  multi->m_pole.M_121 = M_121;
  multi->m_pole.M_112 = M_112;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 4

  /* 5th order terms */
  multi->m_pole.M_005 = M_005;
  multi->m_pole.M_014 = M_014;
  multi->m_pole.M_023 = M_023;
  multi->m_pole.M_032 = M_032;
  multi->m_pole.M_041 = M_041;
  multi->m_pole.M_050 = M_050;
  multi->m_pole.M_104 = M_104;
  multi->m_pole.M_113 = M_113;
  multi->m_pole.M_122 = M_122;
  multi->m_pole.M_131 = M_131;
  multi->m_pole.M_140 = M_140;
  multi->m_pole.M_203 = M_203;
  multi->m_pole.M_212 = M_212;
  multi->m_pole.M_221 = M_221;
  multi->m_pole.M_230 = M_230;
  multi->m_pole.M_302 = M_302;
  multi->m_pole.M_311 = M_311;
  multi->m_pole.M_320 = M_320;
  multi->m_pole.M_401 = M_401;
  multi->m_pole.M_410 = M_410;
  multi->m_pole.M_500 = M_500;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 5
  /* 6th order terms */
  multi->m_pole.M_006 = M_006;
  multi->m_pole.M_015 = M_015;
  multi->m_pole.M_024 = M_024;
  multi->m_pole.M_033 = M_033;
  multi->m_pole.M_042 = M_042;
  multi->m_pole.M_051 = M_051;
  multi->m_pole.M_060 = M_060;
  multi->m_pole.M_105 = M_105;
  multi->m_pole.M_114 = M_114;
  multi->m_pole.M_123 = M_123;
  multi->m_pole.M_132 = M_132;
  multi->m_pole.M_141 = M_141;
  multi->m_pole.M_150 = M_150;
  multi->m_pole.M_204 = M_204;
  multi->m_pole.M_213 = M_213;
  multi->m_pole.M_222 = M_222;
  multi->m_pole.M_231 = M_231;
  multi->m_pole.M_240 = M_240;
  multi->m_pole.M_303 = M_303;
  multi->m_pole.M_312 = M_312;
  multi->m_pole.M_321 = M_321;
  multi->m_pole.M_330 = M_330;
  multi->m_pole.M_402 = M_402;
  multi->m_pole.M_411 = M_411;
  multi->m_pole.M_420 = M_420;
  multi->m_pole.M_501 = M_501;
  multi->m_pole.M_510 = M_510;
  multi->m_pole.M_600 = M_600;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 6
#error "Missing implementation for order >6"
#endif

#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_GRAVITY_FORCE_CHECKS)
  multi->m_pole.num_gpart = gcount;
#endif
}

/**
 * @brief Creates a copy of #multipole shifted to a new location.
 *
 * Corresponds to equation (28d).
 *
 * @param m_a The #multipole copy (content will  be overwritten).
 * @param m_b The #multipole to shift.
 * @param pos_a The position to which m_b will be shifted.
 * @param pos_b The current postion of the multipole to shift.
 */
INLINE static void gravity_M2M(struct multipole *restrict m_a,
                               const struct multipole *restrict m_b,
                               const double pos_a[3], const double pos_b[3]) {

  /* "shift" the softening */
  m_a->max_softening = m_b->max_softening;

  /* Shift 0th order term */
  m_a->M_000 = m_b->M_000;

#if SELF_GRAVITY_MULTIPOLE_ORDER > 0
  const double dx[3] = {pos_a[0] - pos_b[0], pos_a[1] - pos_b[1],
                        pos_a[2] - pos_b[2]};

  /* Shift 1st order term */
  m_a->M_100 = m_b->M_100 + X_100(dx) * m_b->M_000;
  m_a->M_010 = m_b->M_010 + X_010(dx) * m_b->M_000;
  m_a->M_001 = m_b->M_001 + X_001(dx) * m_b->M_000;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 1

  /* Shift 2nd order term */
  m_a->M_200 = m_b->M_200 + X_100(dx) * m_b->M_100 + X_200(dx) * m_b->M_000;
  m_a->M_020 = m_b->M_020 + X_010(dx) * m_b->M_010 + X_020(dx) * m_b->M_000;
  m_a->M_002 = m_b->M_002 + X_001(dx) * m_b->M_001 + X_002(dx) * m_b->M_000;
  m_a->M_110 = m_b->M_110 + X_100(dx) * m_b->M_010 + X_010(dx) * m_b->M_100 +
               X_110(dx) * m_b->M_000;
  m_a->M_101 = m_b->M_101 + X_100(dx) * m_b->M_001 + X_001(dx) * m_b->M_100 +
               X_101(dx) * m_b->M_000;
  m_a->M_011 = m_b->M_011 + X_010(dx) * m_b->M_001 + X_001(dx) * m_b->M_010 +
               X_011(dx) * m_b->M_000;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 2

  /* Shift 3rd order term */
  m_a->M_300 = m_b->M_300 + X_100(dx) * m_b->M_200 + X_200(dx) * m_b->M_100 +
               X_300(dx) * m_b->M_000;
  m_a->M_030 = m_b->M_030 + X_010(dx) * m_b->M_020 + X_020(dx) * m_b->M_010 +
               X_030(dx) * m_b->M_000;
  m_a->M_003 = m_b->M_003 + X_001(dx) * m_b->M_002 + X_002(dx) * m_b->M_001 +
               X_003(dx) * m_b->M_000;
  m_a->M_210 = m_b->M_210 + X_100(dx) * m_b->M_110 + X_010(dx) * m_b->M_200 +
               X_200(dx) * m_b->M_010 + X_110(dx) * m_b->M_100 +
               X_210(dx) * m_b->M_000;
  m_a->M_201 = m_b->M_201 + X_100(dx) * m_b->M_101 + X_001(dx) * m_b->M_200 +
               X_200(dx) * m_b->M_001 + X_101(dx) * m_b->M_100 +
               X_201(dx) * m_b->M_000;
  m_a->M_120 = m_b->M_120 + X_010(dx) * m_b->M_110 + X_100(dx) * m_b->M_020 +
               X_020(dx) * m_b->M_100 + X_110(dx) * m_b->M_010 +
               X_120(dx) * m_b->M_000;
  m_a->M_021 = m_b->M_021 + X_010(dx) * m_b->M_011 + X_001(dx) * m_b->M_020 +
               X_020(dx) * m_b->M_001 + X_011(dx) * m_b->M_010 +
               X_021(dx) * m_b->M_000;
  m_a->M_102 = m_b->M_102 + X_001(dx) * m_b->M_101 + X_100(dx) * m_b->M_002 +
               X_002(dx) * m_b->M_100 + X_101(dx) * m_b->M_001 +
               X_102(dx) * m_b->M_000;
  m_a->M_012 = m_b->M_012 + X_001(dx) * m_b->M_011 + X_010(dx) * m_b->M_002 +
               X_002(dx) * m_b->M_010 + X_011(dx) * m_b->M_001 +
               X_012(dx) * m_b->M_000;
  m_a->M_111 = m_b->M_111 + X_100(dx) * m_b->M_011 + X_010(dx) * m_b->M_101 +
               X_001(dx) * m_b->M_110 + X_110(dx) * m_b->M_001 +
               X_101(dx) * m_b->M_010 + X_011(dx) * m_b->M_100 +
               X_111(dx) * m_b->M_000;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 3
  /* Shift 4th order terms */
  m_a->M_004 = m_b->M_004 + X_001(dx) * m_b->M_003 + X_002(dx) * m_b->M_002 +
               X_003(dx) * m_b->M_001 + X_004(dx) * m_b->M_000;
  m_a->M_013 = m_b->M_013 + X_001(dx) * m_b->M_012 + X_002(dx) * m_b->M_011 +
               X_003(dx) * m_b->M_010 + X_010(dx) * m_b->M_003 +
               X_011(dx) * m_b->M_002 + X_012(dx) * m_b->M_001 +
               X_013(dx) * m_b->M_000;
  m_a->M_022 = m_b->M_022 + X_001(dx) * m_b->M_021 + X_002(dx) * m_b->M_020 +
               X_010(dx) * m_b->M_012 + X_011(dx) * m_b->M_011 +
               X_012(dx) * m_b->M_010 + X_020(dx) * m_b->M_002 +
               X_021(dx) * m_b->M_001 + X_022(dx) * m_b->M_000;
  m_a->M_031 = m_b->M_031 + X_001(dx) * m_b->M_030 + X_010(dx) * m_b->M_021 +
               X_011(dx) * m_b->M_020 + X_020(dx) * m_b->M_011 +
               X_021(dx) * m_b->M_010 + X_030(dx) * m_b->M_001 +
               X_031(dx) * m_b->M_000;
  m_a->M_040 = m_b->M_040 + X_010(dx) * m_b->M_030 + X_020(dx) * m_b->M_020 +
               X_030(dx) * m_b->M_010 + X_040(dx) * m_b->M_000;
  m_a->M_103 = m_b->M_103 + X_001(dx) * m_b->M_102 + X_002(dx) * m_b->M_101 +
               X_003(dx) * m_b->M_100 + X_100(dx) * m_b->M_003 +
               X_101(dx) * m_b->M_002 + X_102(dx) * m_b->M_001 +
               X_103(dx) * m_b->M_000;
  m_a->M_112 =
      m_b->M_112 + X_001(dx) * m_b->M_111 + X_002(dx) * m_b->M_110 +
      X_010(dx) * m_b->M_102 + X_011(dx) * m_b->M_101 + X_012(dx) * m_b->M_100 +
      X_100(dx) * m_b->M_012 + X_101(dx) * m_b->M_011 + X_102(dx) * m_b->M_010 +
      X_110(dx) * m_b->M_002 + X_111(dx) * m_b->M_001 + X_112(dx) * m_b->M_000;
  m_a->M_121 =
      m_b->M_121 + X_001(dx) * m_b->M_120 + X_010(dx) * m_b->M_111 +
      X_011(dx) * m_b->M_110 + X_020(dx) * m_b->M_101 + X_021(dx) * m_b->M_100 +
      X_100(dx) * m_b->M_021 + X_101(dx) * m_b->M_020 + X_110(dx) * m_b->M_011 +
      X_111(dx) * m_b->M_010 + X_120(dx) * m_b->M_001 + X_121(dx) * m_b->M_000;
  m_a->M_130 = m_b->M_130 + X_010(dx) * m_b->M_120 + X_020(dx) * m_b->M_110 +
               X_030(dx) * m_b->M_100 + X_100(dx) * m_b->M_030 +
               X_110(dx) * m_b->M_020 + X_120(dx) * m_b->M_010 +
               X_130(dx) * m_b->M_000;
  m_a->M_202 = m_b->M_202 + X_001(dx) * m_b->M_201 + X_002(dx) * m_b->M_200 +
               X_100(dx) * m_b->M_102 + X_101(dx) * m_b->M_101 +
               X_102(dx) * m_b->M_100 + X_200(dx) * m_b->M_002 +
               X_201(dx) * m_b->M_001 + X_202(dx) * m_b->M_000;
  m_a->M_211 =
      m_b->M_211 + X_001(dx) * m_b->M_210 + X_010(dx) * m_b->M_201 +
      X_011(dx) * m_b->M_200 + X_100(dx) * m_b->M_111 + X_101(dx) * m_b->M_110 +
      X_110(dx) * m_b->M_101 + X_111(dx) * m_b->M_100 + X_200(dx) * m_b->M_011 +
      X_201(dx) * m_b->M_010 + X_210(dx) * m_b->M_001 + X_211(dx) * m_b->M_000;
  m_a->M_220 = m_b->M_220 + X_010(dx) * m_b->M_210 + X_020(dx) * m_b->M_200 +
               X_100(dx) * m_b->M_120 + X_110(dx) * m_b->M_110 +
               X_120(dx) * m_b->M_100 + X_200(dx) * m_b->M_020 +
               X_210(dx) * m_b->M_010 + X_220(dx) * m_b->M_000;
  m_a->M_301 = m_b->M_301 + X_001(dx) * m_b->M_300 + X_100(dx) * m_b->M_201 +
               X_101(dx) * m_b->M_200 + X_200(dx) * m_b->M_101 +
               X_201(dx) * m_b->M_100 + X_300(dx) * m_b->M_001 +
               X_301(dx) * m_b->M_000;
  m_a->M_310 = m_b->M_310 + X_010(dx) * m_b->M_300 + X_100(dx) * m_b->M_210 +
               X_110(dx) * m_b->M_200 + X_200(dx) * m_b->M_110 +
               X_210(dx) * m_b->M_100 + X_300(dx) * m_b->M_010 +
               X_310(dx) * m_b->M_000;
  m_a->M_400 = m_b->M_400 + X_100(dx) * m_b->M_300 + X_200(dx) * m_b->M_200 +
               X_300(dx) * m_b->M_100 + X_400(dx) * m_b->M_000;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 4
  /* Shift 5th order terms */
  m_a->M_005 = m_b->M_005 + X_001(dx) * m_b->M_004 + X_002(dx) * m_b->M_003 +
               X_003(dx) * m_b->M_002 + X_004(dx) * m_b->M_001 +
               X_005(dx) * m_b->M_000;
  m_a->M_014 = m_b->M_014 + X_001(dx) * m_b->M_013 + X_002(dx) * m_b->M_012 +
               X_003(dx) * m_b->M_011 + X_004(dx) * m_b->M_010 +
               X_010(dx) * m_b->M_004 + X_011(dx) * m_b->M_003 +
               X_012(dx) * m_b->M_002 + X_013(dx) * m_b->M_001 +
               X_014(dx) * m_b->M_000;
  m_a->M_023 =
      m_b->M_023 + X_001(dx) * m_b->M_022 + X_002(dx) * m_b->M_021 +
      X_003(dx) * m_b->M_020 + X_010(dx) * m_b->M_013 + X_011(dx) * m_b->M_012 +
      X_012(dx) * m_b->M_011 + X_013(dx) * m_b->M_010 + X_020(dx) * m_b->M_003 +
      X_021(dx) * m_b->M_002 + X_022(dx) * m_b->M_001 + X_023(dx) * m_b->M_000;
  m_a->M_032 =
      m_b->M_032 + X_001(dx) * m_b->M_031 + X_002(dx) * m_b->M_030 +
      X_010(dx) * m_b->M_022 + X_011(dx) * m_b->M_021 + X_012(dx) * m_b->M_020 +
      X_020(dx) * m_b->M_012 + X_021(dx) * m_b->M_011 + X_022(dx) * m_b->M_010 +
      X_030(dx) * m_b->M_002 + X_031(dx) * m_b->M_001 + X_032(dx) * m_b->M_000;
  m_a->M_041 = m_b->M_041 + X_001(dx) * m_b->M_040 + X_010(dx) * m_b->M_031 +
               X_011(dx) * m_b->M_030 + X_020(dx) * m_b->M_021 +
               X_021(dx) * m_b->M_020 + X_030(dx) * m_b->M_011 +
               X_031(dx) * m_b->M_010 + X_040(dx) * m_b->M_001 +
               X_041(dx) * m_b->M_000;
  m_a->M_050 = m_b->M_050 + X_010(dx) * m_b->M_040 + X_020(dx) * m_b->M_030 +
               X_030(dx) * m_b->M_020 + X_040(dx) * m_b->M_010 +
               X_050(dx) * m_b->M_000;
  m_a->M_104 = m_b->M_104 + X_001(dx) * m_b->M_103 + X_002(dx) * m_b->M_102 +
               X_003(dx) * m_b->M_101 + X_004(dx) * m_b->M_100 +
               X_100(dx) * m_b->M_004 + X_101(dx) * m_b->M_003 +
               X_102(dx) * m_b->M_002 + X_103(dx) * m_b->M_001 +
               X_104(dx) * m_b->M_000;
  m_a->M_113 =
      m_b->M_113 + X_001(dx) * m_b->M_112 + X_002(dx) * m_b->M_111 +
      X_003(dx) * m_b->M_110 + X_010(dx) * m_b->M_103 + X_011(dx) * m_b->M_102 +
      X_012(dx) * m_b->M_101 + X_013(dx) * m_b->M_100 + X_100(dx) * m_b->M_013 +
      X_101(dx) * m_b->M_012 + X_102(dx) * m_b->M_011 + X_103(dx) * m_b->M_010 +
      X_110(dx) * m_b->M_003 + X_111(dx) * m_b->M_002 + X_112(dx) * m_b->M_001 +
      X_113(dx) * m_b->M_000;
  m_a->M_122 =
      m_b->M_122 + X_001(dx) * m_b->M_121 + X_002(dx) * m_b->M_120 +
      X_010(dx) * m_b->M_112 + X_011(dx) * m_b->M_111 + X_012(dx) * m_b->M_110 +
      X_020(dx) * m_b->M_102 + X_021(dx) * m_b->M_101 + X_022(dx) * m_b->M_100 +
      X_100(dx) * m_b->M_022 + X_101(dx) * m_b->M_021 + X_102(dx) * m_b->M_020 +
      X_110(dx) * m_b->M_012 + X_111(dx) * m_b->M_011 + X_112(dx) * m_b->M_010 +
      X_120(dx) * m_b->M_002 + X_121(dx) * m_b->M_001 + X_122(dx) * m_b->M_000;
  m_a->M_131 =
      m_b->M_131 + X_001(dx) * m_b->M_130 + X_010(dx) * m_b->M_121 +
      X_011(dx) * m_b->M_120 + X_020(dx) * m_b->M_111 + X_021(dx) * m_b->M_110 +
      X_030(dx) * m_b->M_101 + X_031(dx) * m_b->M_100 + X_100(dx) * m_b->M_031 +
      X_101(dx) * m_b->M_030 + X_110(dx) * m_b->M_021 + X_111(dx) * m_b->M_020 +
      X_120(dx) * m_b->M_011 + X_121(dx) * m_b->M_010 + X_130(dx) * m_b->M_001 +
      X_131(dx) * m_b->M_000;
  m_a->M_140 = m_b->M_140 + X_010(dx) * m_b->M_130 + X_020(dx) * m_b->M_120 +
               X_030(dx) * m_b->M_110 + X_040(dx) * m_b->M_100 +
               X_100(dx) * m_b->M_040 + X_110(dx) * m_b->M_030 +
               X_120(dx) * m_b->M_020 + X_130(dx) * m_b->M_010 +
               X_140(dx) * m_b->M_000;
  m_a->M_203 =
      m_b->M_203 + X_001(dx) * m_b->M_202 + X_002(dx) * m_b->M_201 +
      X_003(dx) * m_b->M_200 + X_100(dx) * m_b->M_103 + X_101(dx) * m_b->M_102 +
      X_102(dx) * m_b->M_101 + X_103(dx) * m_b->M_100 + X_200(dx) * m_b->M_003 +
      X_201(dx) * m_b->M_002 + X_202(dx) * m_b->M_001 + X_203(dx) * m_b->M_000;
  m_a->M_212 =
      m_b->M_212 + X_001(dx) * m_b->M_211 + X_002(dx) * m_b->M_210 +
      X_010(dx) * m_b->M_202 + X_011(dx) * m_b->M_201 + X_012(dx) * m_b->M_200 +
      X_100(dx) * m_b->M_112 + X_101(dx) * m_b->M_111 + X_102(dx) * m_b->M_110 +
      X_110(dx) * m_b->M_102 + X_111(dx) * m_b->M_101 + X_112(dx) * m_b->M_100 +
      X_200(dx) * m_b->M_012 + X_201(dx) * m_b->M_011 + X_202(dx) * m_b->M_010 +
      X_210(dx) * m_b->M_002 + X_211(dx) * m_b->M_001 + X_212(dx) * m_b->M_000;
  m_a->M_221 =
      m_b->M_221 + X_001(dx) * m_b->M_220 + X_010(dx) * m_b->M_211 +
      X_011(dx) * m_b->M_210 + X_020(dx) * m_b->M_201 + X_021(dx) * m_b->M_200 +
      X_100(dx) * m_b->M_121 + X_101(dx) * m_b->M_120 + X_110(dx) * m_b->M_111 +
      X_111(dx) * m_b->M_110 + X_120(dx) * m_b->M_101 + X_121(dx) * m_b->M_100 +
      X_200(dx) * m_b->M_021 + X_201(dx) * m_b->M_020 + X_210(dx) * m_b->M_011 +
      X_211(dx) * m_b->M_010 + X_220(dx) * m_b->M_001 + X_221(dx) * m_b->M_000;
  m_a->M_230 =
      m_b->M_230 + X_010(dx) * m_b->M_220 + X_020(dx) * m_b->M_210 +
      X_030(dx) * m_b->M_200 + X_100(dx) * m_b->M_130 + X_110(dx) * m_b->M_120 +
      X_120(dx) * m_b->M_110 + X_130(dx) * m_b->M_100 + X_200(dx) * m_b->M_030 +
      X_210(dx) * m_b->M_020 + X_220(dx) * m_b->M_010 + X_230(dx) * m_b->M_000;
  m_a->M_302 =
      m_b->M_302 + X_001(dx) * m_b->M_301 + X_002(dx) * m_b->M_300 +
      X_100(dx) * m_b->M_202 + X_101(dx) * m_b->M_201 + X_102(dx) * m_b->M_200 +
      X_200(dx) * m_b->M_102 + X_201(dx) * m_b->M_101 + X_202(dx) * m_b->M_100 +
      X_300(dx) * m_b->M_002 + X_301(dx) * m_b->M_001 + X_302(dx) * m_b->M_000;
  m_a->M_311 =
      m_b->M_311 + X_001(dx) * m_b->M_310 + X_010(dx) * m_b->M_301 +
      X_011(dx) * m_b->M_300 + X_100(dx) * m_b->M_211 + X_101(dx) * m_b->M_210 +
      X_110(dx) * m_b->M_201 + X_111(dx) * m_b->M_200 + X_200(dx) * m_b->M_111 +
      X_201(dx) * m_b->M_110 + X_210(dx) * m_b->M_101 + X_211(dx) * m_b->M_100 +
      X_300(dx) * m_b->M_011 + X_301(dx) * m_b->M_010 + X_310(dx) * m_b->M_001 +
      X_311(dx) * m_b->M_000;
  m_a->M_320 =
      m_b->M_320 + X_010(dx) * m_b->M_310 + X_020(dx) * m_b->M_300 +
      X_100(dx) * m_b->M_220 + X_110(dx) * m_b->M_210 + X_120(dx) * m_b->M_200 +
      X_200(dx) * m_b->M_120 + X_210(dx) * m_b->M_110 + X_220(dx) * m_b->M_100 +
      X_300(dx) * m_b->M_020 + X_310(dx) * m_b->M_010 + X_320(dx) * m_b->M_000;
  m_a->M_401 = m_b->M_401 + X_001(dx) * m_b->M_400 + X_100(dx) * m_b->M_301 +
               X_101(dx) * m_b->M_300 + X_200(dx) * m_b->M_201 +
               X_201(dx) * m_b->M_200 + X_300(dx) * m_b->M_101 +
               X_301(dx) * m_b->M_100 + X_400(dx) * m_b->M_001 +
               X_401(dx) * m_b->M_000;
  m_a->M_410 = m_b->M_410 + X_010(dx) * m_b->M_400 + X_100(dx) * m_b->M_310 +
               X_110(dx) * m_b->M_300 + X_200(dx) * m_b->M_210 +
               X_210(dx) * m_b->M_200 + X_300(dx) * m_b->M_110 +
               X_310(dx) * m_b->M_100 + X_400(dx) * m_b->M_010 +
               X_410(dx) * m_b->M_000;
  m_a->M_500 = m_b->M_500 + X_100(dx) * m_b->M_400 + X_200(dx) * m_b->M_300 +
               X_300(dx) * m_b->M_200 + X_400(dx) * m_b->M_100 +
               X_500(dx) * m_b->M_000;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 5
  /* Shift 6th order terms */
  m_a->M_006 =
    m_b->M_006 + X_001(dx) * m_b->M_005 + X_002(dx) * m_b->M_004
    + X_003(dx) * m_b->M_003 + X_004(dx) * m_b->M_002
    + X_005(dx) * m_b->M_001 + X_006(dx) * m_b->M_000;
  m_a->M_015 =
    m_b->M_015 + X_001(dx) * m_b->M_014 + X_002(dx) * m_b->M_013
    + X_003(dx) * m_b->M_012 + X_004(dx) * m_b->M_011
    + X_005(dx) * m_b->M_010 + X_010(dx) * m_b->M_005
    + X_011(dx) * m_b->M_004 + X_012(dx) * m_b->M_003
    + X_013(dx) * m_b->M_002 + X_014(dx) * m_b->M_001
    + X_015(dx) * m_b->M_000;
  m_a->M_024 =
    m_b->M_024 + X_001(dx) * m_b->M_023 + X_002(dx) * m_b->M_022
    + X_003(dx) * m_b->M_021 + X_004(dx) * m_b->M_020
    + X_010(dx) * m_b->M_014 + X_011(dx) * m_b->M_013
    + X_012(dx) * m_b->M_012 + X_013(dx) * m_b->M_011
    + X_014(dx) * m_b->M_010 + X_020(dx) * m_b->M_004
    + X_021(dx) * m_b->M_003 + X_022(dx) * m_b->M_002
    + X_023(dx) * m_b->M_001 + X_024(dx) * m_b->M_000;
  m_a->M_033 =
    m_b->M_033 + X_001(dx) * m_b->M_032 + X_002(dx) * m_b->M_031
    + X_003(dx) * m_b->M_030 + X_010(dx) * m_b->M_023
    + X_011(dx) * m_b->M_022 + X_012(dx) * m_b->M_021
    + X_013(dx) * m_b->M_020 + X_020(dx) * m_b->M_013
    + X_021(dx) * m_b->M_012 + X_022(dx) * m_b->M_011
    + X_023(dx) * m_b->M_010 + X_030(dx) * m_b->M_003
    + X_031(dx) * m_b->M_002 + X_032(dx) * m_b->M_001
    + X_033(dx) * m_b->M_000;
  m_a->M_042 =
    m_b->M_042 + X_001(dx) * m_b->M_041 + X_002(dx) * m_b->M_040
    + X_010(dx) * m_b->M_032 + X_011(dx) * m_b->M_031
    + X_012(dx) * m_b->M_030 + X_020(dx) * m_b->M_022
    + X_021(dx) * m_b->M_021 + X_022(dx) * m_b->M_020
    + X_030(dx) * m_b->M_012 + X_031(dx) * m_b->M_011
    + X_032(dx) * m_b->M_010 + X_040(dx) * m_b->M_002
    + X_041(dx) * m_b->M_001 + X_042(dx) * m_b->M_000;
  m_a->M_051 =
    m_b->M_051 + X_001(dx) * m_b->M_050 + X_010(dx) * m_b->M_041
    + X_011(dx) * m_b->M_040 + X_020(dx) * m_b->M_031
    + X_021(dx) * m_b->M_030 + X_030(dx) * m_b->M_021
    + X_031(dx) * m_b->M_020 + X_040(dx) * m_b->M_011
    + X_041(dx) * m_b->M_010 + X_050(dx) * m_b->M_001
    + X_051(dx) * m_b->M_000;
  m_a->M_060 =
    m_b->M_060 + X_010(dx) * m_b->M_050 + X_020(dx) * m_b->M_040
    + X_030(dx) * m_b->M_030 + X_040(dx) * m_b->M_020
    + X_050(dx) * m_b->M_010 + X_060(dx) * m_b->M_000;
  m_a->M_105 =
    m_b->M_105 + X_001(dx) * m_b->M_104 + X_002(dx) * m_b->M_103
    + X_003(dx) * m_b->M_102 + X_004(dx) * m_b->M_101
    + X_005(dx) * m_b->M_100 + X_100(dx) * m_b->M_005
    + X_101(dx) * m_b->M_004 + X_102(dx) * m_b->M_003
    + X_103(dx) * m_b->M_002 + X_104(dx) * m_b->M_001
    + X_105(dx) * m_b->M_000;
  m_a->M_114 =
    m_b->M_114 + X_001(dx) * m_b->M_113 + X_002(dx) * m_b->M_112
    + X_003(dx) * m_b->M_111 + X_004(dx) * m_b->M_110
    + X_010(dx) * m_b->M_104 + X_011(dx) * m_b->M_103
    + X_012(dx) * m_b->M_102 + X_013(dx) * m_b->M_101
    + X_014(dx) * m_b->M_100 + X_100(dx) * m_b->M_014
    + X_101(dx) * m_b->M_013 + X_102(dx) * m_b->M_012
    + X_103(dx) * m_b->M_011 + X_104(dx) * m_b->M_010
    + X_110(dx) * m_b->M_004 + X_111(dx) * m_b->M_003
    + X_112(dx) * m_b->M_002 + X_113(dx) * m_b->M_001
    + X_114(dx) * m_b->M_000;
  m_a->M_123 =
    m_b->M_123 + X_001(dx) * m_b->M_122 + X_002(dx) * m_b->M_121
    + X_003(dx) * m_b->M_120 + X_010(dx) * m_b->M_113
    + X_011(dx) * m_b->M_112 + X_012(dx) * m_b->M_111
    + X_013(dx) * m_b->M_110 + X_020(dx) * m_b->M_103
    + X_021(dx) * m_b->M_102 + X_022(dx) * m_b->M_101
    + X_023(dx) * m_b->M_100 + X_100(dx) * m_b->M_023
    + X_101(dx) * m_b->M_022 + X_102(dx) * m_b->M_021
    + X_103(dx) * m_b->M_020 + X_110(dx) * m_b->M_013
    + X_111(dx) * m_b->M_012 + X_112(dx) * m_b->M_011
    + X_113(dx) * m_b->M_010 + X_120(dx) * m_b->M_003
    + X_121(dx) * m_b->M_002 + X_122(dx) * m_b->M_001
    + X_123(dx) * m_b->M_000;
  m_a->M_132 =
    m_b->M_132 + X_001(dx) * m_b->M_131 + X_002(dx) * m_b->M_130
    + X_010(dx) * m_b->M_122 + X_011(dx) * m_b->M_121
    + X_012(dx) * m_b->M_120 + X_020(dx) * m_b->M_112
    + X_021(dx) * m_b->M_111 + X_022(dx) * m_b->M_110
    + X_030(dx) * m_b->M_102 + X_031(dx) * m_b->M_101
    + X_032(dx) * m_b->M_100 + X_100(dx) * m_b->M_032
    + X_101(dx) * m_b->M_031 + X_102(dx) * m_b->M_030
    + X_110(dx) * m_b->M_022 + X_111(dx) * m_b->M_021
    + X_112(dx) * m_b->M_020 + X_120(dx) * m_b->M_012
    + X_121(dx) * m_b->M_011 + X_122(dx) * m_b->M_010
    + X_130(dx) * m_b->M_002 + X_131(dx) * m_b->M_001
    + X_132(dx) * m_b->M_000;
  m_a->M_141 =
    m_b->M_141 + X_001(dx) * m_b->M_140 + X_010(dx) * m_b->M_131
    + X_011(dx) * m_b->M_130 + X_020(dx) * m_b->M_121
    + X_021(dx) * m_b->M_120 + X_030(dx) * m_b->M_111
    + X_031(dx) * m_b->M_110 + X_040(dx) * m_b->M_101
    + X_041(dx) * m_b->M_100 + X_100(dx) * m_b->M_041
    + X_101(dx) * m_b->M_040 + X_110(dx) * m_b->M_031
    + X_111(dx) * m_b->M_030 + X_120(dx) * m_b->M_021
    + X_121(dx) * m_b->M_020 + X_130(dx) * m_b->M_011
    + X_131(dx) * m_b->M_010 + X_140(dx) * m_b->M_001
    + X_141(dx) * m_b->M_000;
  m_a->M_150 =
    m_b->M_150 + X_010(dx) * m_b->M_140 + X_020(dx) * m_b->M_130
    + X_030(dx) * m_b->M_120 + X_040(dx) * m_b->M_110
    + X_050(dx) * m_b->M_100 + X_100(dx) * m_b->M_050
    + X_110(dx) * m_b->M_040 + X_120(dx) * m_b->M_030
    + X_130(dx) * m_b->M_020 + X_140(dx) * m_b->M_010
    + X_150(dx) * m_b->M_000;
  m_a->M_204 =
    m_b->M_204 + X_001(dx) * m_b->M_203 + X_002(dx) * m_b->M_202
    + X_003(dx) * m_b->M_201 + X_004(dx) * m_b->M_200
    + X_100(dx) * m_b->M_104 + X_101(dx) * m_b->M_103
    + X_102(dx) * m_b->M_102 + X_103(dx) * m_b->M_101
    + X_104(dx) * m_b->M_100 + X_200(dx) * m_b->M_004
    + X_201(dx) * m_b->M_003 + X_202(dx) * m_b->M_002
    + X_203(dx) * m_b->M_001 + X_204(dx) * m_b->M_000;
  m_a->M_213 =
    m_b->M_213 + X_001(dx) * m_b->M_212 + X_002(dx) * m_b->M_211
    + X_003(dx) * m_b->M_210 + X_010(dx) * m_b->M_203
    + X_011(dx) * m_b->M_202 + X_012(dx) * m_b->M_201
    + X_013(dx) * m_b->M_200 + X_100(dx) * m_b->M_113
    + X_101(dx) * m_b->M_112 + X_102(dx) * m_b->M_111
    + X_103(dx) * m_b->M_110 + X_110(dx) * m_b->M_103
    + X_111(dx) * m_b->M_102 + X_112(dx) * m_b->M_101
    + X_113(dx) * m_b->M_100 + X_200(dx) * m_b->M_013
    + X_201(dx) * m_b->M_012 + X_202(dx) * m_b->M_011
    + X_203(dx) * m_b->M_010 + X_210(dx) * m_b->M_003
    + X_211(dx) * m_b->M_002 + X_212(dx) * m_b->M_001
    + X_213(dx) * m_b->M_000;
  m_a->M_222 =
    m_b->M_222 + X_001(dx) * m_b->M_221 + X_002(dx) * m_b->M_220
    + X_010(dx) * m_b->M_212 + X_011(dx) * m_b->M_211
    + X_012(dx) * m_b->M_210 + X_020(dx) * m_b->M_202
    + X_021(dx) * m_b->M_201 + X_022(dx) * m_b->M_200
    + X_100(dx) * m_b->M_122 + X_101(dx) * m_b->M_121
    + X_102(dx) * m_b->M_120 + X_110(dx) * m_b->M_112
    + X_111(dx) * m_b->M_111 + X_112(dx) * m_b->M_110
    + X_120(dx) * m_b->M_102 + X_121(dx) * m_b->M_101
    + X_122(dx) * m_b->M_100 + X_200(dx) * m_b->M_022
    + X_201(dx) * m_b->M_021 + X_202(dx) * m_b->M_020
    + X_210(dx) * m_b->M_012 + X_211(dx) * m_b->M_011
    + X_212(dx) * m_b->M_010 + X_220(dx) * m_b->M_002
    + X_221(dx) * m_b->M_001 + X_222(dx) * m_b->M_000;
  m_a->M_231 =
    m_b->M_231 + X_001(dx) * m_b->M_230 + X_010(dx) * m_b->M_221
    + X_011(dx) * m_b->M_220 + X_020(dx) * m_b->M_211
    + X_021(dx) * m_b->M_210 + X_030(dx) * m_b->M_201
    + X_031(dx) * m_b->M_200 + X_100(dx) * m_b->M_131
    + X_101(dx) * m_b->M_130 + X_110(dx) * m_b->M_121
    + X_111(dx) * m_b->M_120 + X_120(dx) * m_b->M_111
    + X_121(dx) * m_b->M_110 + X_130(dx) * m_b->M_101
    + X_131(dx) * m_b->M_100 + X_200(dx) * m_b->M_031
    + X_201(dx) * m_b->M_030 + X_210(dx) * m_b->M_021
    + X_211(dx) * m_b->M_020 + X_220(dx) * m_b->M_011
    + X_221(dx) * m_b->M_010 + X_230(dx) * m_b->M_001
    + X_231(dx) * m_b->M_000;
  m_a->M_240 =
    m_b->M_240 + X_010(dx) * m_b->M_230 + X_020(dx) * m_b->M_220
    + X_030(dx) * m_b->M_210 + X_040(dx) * m_b->M_200
    + X_100(dx) * m_b->M_140 + X_110(dx) * m_b->M_130
    + X_120(dx) * m_b->M_120 + X_130(dx) * m_b->M_110
    + X_140(dx) * m_b->M_100 + X_200(dx) * m_b->M_040
    + X_210(dx) * m_b->M_030 + X_220(dx) * m_b->M_020
    + X_230(dx) * m_b->M_010 + X_240(dx) * m_b->M_000;
  m_a->M_303 =
    m_b->M_303 + X_001(dx) * m_b->M_302 + X_002(dx) * m_b->M_301
    + X_003(dx) * m_b->M_300 + X_100(dx) * m_b->M_203
    + X_101(dx) * m_b->M_202 + X_102(dx) * m_b->M_201
    + X_103(dx) * m_b->M_200 + X_200(dx) * m_b->M_103
    + X_201(dx) * m_b->M_102 + X_202(dx) * m_b->M_101
    + X_203(dx) * m_b->M_100 + X_300(dx) * m_b->M_003
    + X_301(dx) * m_b->M_002 + X_302(dx) * m_b->M_001
    + X_303(dx) * m_b->M_000;
  m_a->M_312 =
    m_b->M_312 + X_001(dx) * m_b->M_311 + X_002(dx) * m_b->M_310
    + X_010(dx) * m_b->M_302 + X_011(dx) * m_b->M_301
    + X_012(dx) * m_b->M_300 + X_100(dx) * m_b->M_212
    + X_101(dx) * m_b->M_211 + X_102(dx) * m_b->M_210
    + X_110(dx) * m_b->M_202 + X_111(dx) * m_b->M_201
    + X_112(dx) * m_b->M_200 + X_200(dx) * m_b->M_112
    + X_201(dx) * m_b->M_111 + X_202(dx) * m_b->M_110
    + X_210(dx) * m_b->M_102 + X_211(dx) * m_b->M_101
    + X_212(dx) * m_b->M_100 + X_300(dx) * m_b->M_012
    + X_301(dx) * m_b->M_011 + X_302(dx) * m_b->M_010
    + X_310(dx) * m_b->M_002 + X_311(dx) * m_b->M_001
    + X_312(dx) * m_b->M_000;
  m_a->M_321 =
    m_b->M_321 + X_001(dx) * m_b->M_320 + X_010(dx) * m_b->M_311
    + X_011(dx) * m_b->M_310 + X_020(dx) * m_b->M_301
    + X_021(dx) * m_b->M_300 + X_100(dx) * m_b->M_221
    + X_101(dx) * m_b->M_220 + X_110(dx) * m_b->M_211
    + X_111(dx) * m_b->M_210 + X_120(dx) * m_b->M_201
    + X_121(dx) * m_b->M_200 + X_200(dx) * m_b->M_121
    + X_201(dx) * m_b->M_120 + X_210(dx) * m_b->M_111
    + X_211(dx) * m_b->M_110 + X_220(dx) * m_b->M_101
    + X_221(dx) * m_b->M_100 + X_300(dx) * m_b->M_021
    + X_301(dx) * m_b->M_020 + X_310(dx) * m_b->M_011
    + X_311(dx) * m_b->M_010 + X_320(dx) * m_b->M_001
    + X_321(dx) * m_b->M_000;
  m_a->M_330 =
    m_b->M_330 + X_010(dx) * m_b->M_320 + X_020(dx) * m_b->M_310
    + X_030(dx) * m_b->M_300 + X_100(dx) * m_b->M_230
    + X_110(dx) * m_b->M_220 + X_120(dx) * m_b->M_210
    + X_130(dx) * m_b->M_200 + X_200(dx) * m_b->M_130
    + X_210(dx) * m_b->M_120 + X_220(dx) * m_b->M_110
    + X_230(dx) * m_b->M_100 + X_300(dx) * m_b->M_030
    + X_310(dx) * m_b->M_020 + X_320(dx) * m_b->M_010
    + X_330(dx) * m_b->M_000;
  m_a->M_402 =
    m_b->M_402 + X_001(dx) * m_b->M_401 + X_002(dx) * m_b->M_400
    + X_100(dx) * m_b->M_302 + X_101(dx) * m_b->M_301
    + X_102(dx) * m_b->M_300 + X_200(dx) * m_b->M_202
    + X_201(dx) * m_b->M_201 + X_202(dx) * m_b->M_200
    + X_300(dx) * m_b->M_102 + X_301(dx) * m_b->M_101
    + X_302(dx) * m_b->M_100 + X_400(dx) * m_b->M_002
    + X_401(dx) * m_b->M_001 + X_402(dx) * m_b->M_000;
  m_a->M_411 =
    m_b->M_411 + X_001(dx) * m_b->M_410 + X_010(dx) * m_b->M_401
    + X_011(dx) * m_b->M_400 + X_100(dx) * m_b->M_311
    + X_101(dx) * m_b->M_310 + X_110(dx) * m_b->M_301
    + X_111(dx) * m_b->M_300 + X_200(dx) * m_b->M_211
    + X_201(dx) * m_b->M_210 + X_210(dx) * m_b->M_201
    + X_211(dx) * m_b->M_200 + X_300(dx) * m_b->M_111
    + X_301(dx) * m_b->M_110 + X_310(dx) * m_b->M_101
    + X_311(dx) * m_b->M_100 + X_400(dx) * m_b->M_011
    + X_401(dx) * m_b->M_010 + X_410(dx) * m_b->M_001
    + X_411(dx) * m_b->M_000;
  m_a->M_420 =
    m_b->M_420 + X_010(dx) * m_b->M_410 + X_020(dx) * m_b->M_400
    + X_100(dx) * m_b->M_320 + X_110(dx) * m_b->M_310
    + X_120(dx) * m_b->M_300 + X_200(dx) * m_b->M_220
    + X_210(dx) * m_b->M_210 + X_220(dx) * m_b->M_200
    + X_300(dx) * m_b->M_120 + X_310(dx) * m_b->M_110
    + X_320(dx) * m_b->M_100 + X_400(dx) * m_b->M_020
    + X_410(dx) * m_b->M_010 + X_420(dx) * m_b->M_000;
  m_a->M_501 =
    m_b->M_501 + X_001(dx) * m_b->M_500 + X_100(dx) * m_b->M_401
    + X_101(dx) * m_b->M_400 + X_200(dx) * m_b->M_301
    + X_201(dx) * m_b->M_300 + X_300(dx) * m_b->M_201
    + X_301(dx) * m_b->M_200 + X_400(dx) * m_b->M_101
    + X_401(dx) * m_b->M_100 + X_500(dx) * m_b->M_001
    + X_501(dx) * m_b->M_000;
  m_a->M_510 =
    m_b->M_510 + X_010(dx) * m_b->M_500 + X_100(dx) * m_b->M_410
    + X_110(dx) * m_b->M_400 + X_200(dx) * m_b->M_310
    + X_210(dx) * m_b->M_300 + X_300(dx) * m_b->M_210
    + X_310(dx) * m_b->M_200 + X_400(dx) * m_b->M_110
    + X_410(dx) * m_b->M_100 + X_500(dx) * m_b->M_010
    + X_510(dx) * m_b->M_000;
  m_a->M_600 =
    m_b->M_600 + X_100(dx) * m_b->M_500 + X_200(dx) * m_b->M_400
    + X_300(dx) * m_b->M_300 + X_400(dx) * m_b->M_200
    + X_500(dx) * m_b->M_100 + X_600(dx) * m_b->M_000;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 6
#error "Missing implementation for order >6"
#endif

#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_GRAVITY_FORCE_CHECKS)
  m_a->num_gpart = m_b->num_gpart;
#endif
}

/**
 * @brief Compute the field tensors due to a multipole.
 *
 * Corresponds to equation (28b).
 *
 * @param l_b The field tensor to compute.
 * @param m_a The multipole creating the field.
 * @param pot The derivatives of the potential.
 */
INLINE static void gravity_M2L_apply(
    struct grav_tensor *restrict l_b, const struct multipole *restrict m_a,
    const struct potential_derivatives_M2L *pot) {

#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_GRAVITY_FORCE_CHECKS)
  /* Count interactions */
  l_b->num_interacted += m_a->num_gpart;
#endif

  /* Record that this tensor has received contributions */
  l_b->interacted = 1;

  const float M_000 = m_a->M_000;
  const float D_000 = pot->D_000;

  /*  0th order term */
  l_b->F_000 += M_000 * D_000;

#if SELF_GRAVITY_MULTIPOLE_ORDER > 0

  /* The dipole term is zero when using the CoM */
  /* The compiler will optimize out the terms in the equations */
  /* below. We keep them written to maintain the logical structure. */
  const float M_100 = 0.f;
  const float M_010 = 0.f;
  const float M_001 = 0.f;

  const float D_100 = pot->D_100;
  const float D_010 = pot->D_010;
  const float D_001 = pot->D_001;

  /*  1st order multipole term (addition to rank 0)*/
  l_b->F_000 += M_100 * D_100 + M_010 * D_010 + M_001 * D_001;

  /*  1st order multipole term (addition to rank 1)*/
  l_b->F_100 += M_000 * D_100;
  l_b->F_010 += M_000 * D_010;
  l_b->F_001 += M_000 * D_001;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 1

  const float M_200 = m_a->M_200;
  const float M_020 = m_a->M_020;
  const float M_002 = m_a->M_002;
  const float M_110 = m_a->M_110;
  const float M_101 = m_a->M_101;
  const float M_011 = m_a->M_011;

  const float D_200 = pot->D_200;
  const float D_020 = pot->D_020;
  const float D_002 = pot->D_002;
  const float D_110 = pot->D_110;
  const float D_101 = pot->D_101;
  const float D_011 = pot->D_011;

  /*  2nd order multipole term (addition to rank 0)*/
  l_b->F_000 += M_200 * D_200 + M_020 * D_020 + M_002 * D_002;
  l_b->F_000 += M_110 * D_110 + M_101 * D_101 + M_011 * D_011;

  /*  2nd order multipole term (addition to rank 1)*/
  l_b->F_100 += M_100 * D_200 + M_010 * D_110 + M_001 * D_101;
  l_b->F_010 += M_100 * D_110 + M_010 * D_020 + M_001 * D_011;
  l_b->F_001 += M_100 * D_101 + M_010 * D_011 + M_001 * D_002;

  /*  2nd order multipole term (addition to rank 2)*/
  l_b->F_200 += M_000 * D_200;
  l_b->F_020 += M_000 * D_020;
  l_b->F_002 += M_000 * D_002;
  l_b->F_110 += M_000 * D_110;
  l_b->F_101 += M_000 * D_101;
  l_b->F_011 += M_000 * D_011;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 2

  const float M_300 = m_a->M_300;
  const float M_030 = m_a->M_030;
  const float M_003 = m_a->M_003;
  const float M_210 = m_a->M_210;
  const float M_201 = m_a->M_201;
  const float M_021 = m_a->M_021;
  const float M_120 = m_a->M_120;
  const float M_012 = m_a->M_012;
  const float M_102 = m_a->M_102;
  const float M_111 = m_a->M_111;

  const float D_300 = pot->D_300;
  const float D_030 = pot->D_030;
  const float D_003 = pot->D_003;
  const float D_210 = pot->D_210;
  const float D_201 = pot->D_201;
  const float D_021 = pot->D_021;
  const float D_120 = pot->D_120;
  const float D_012 = pot->D_012;
  const float D_102 = pot->D_102;
  const float D_111 = pot->D_111;

  /*  3rd order multipole term (addition to rank 0)*/
  l_b->F_000 += M_300 * D_300 + M_030 * D_030 + M_003 * D_003;
  l_b->F_000 += M_210 * D_210 + M_201 * D_201 + M_120 * D_120;
  l_b->F_000 += M_021 * D_021 + M_102 * D_102 + M_012 * D_012;
  l_b->F_000 += M_111 * D_111;

  /*  3rd order multipole term (addition to rank 1)*/
  l_b->F_100 += M_200 * D_300 + M_020 * D_120 + M_002 * D_102;
  l_b->F_100 += M_110 * D_210 + M_101 * D_201 + M_011 * D_111;
  l_b->F_010 += M_200 * D_210 + M_020 * D_030 + M_002 * D_012;
  l_b->F_010 += M_110 * D_120 + M_101 * D_111 + M_011 * D_021;
  l_b->F_001 += M_200 * D_201 + M_020 * D_021 + M_002 * D_003;
  l_b->F_001 += M_110 * D_111 + M_101 * D_102 + M_011 * D_012;

  /*  3rd order multipole term (addition to rank 2)*/
  l_b->F_200 += M_100 * D_300 + M_010 * D_210 + M_001 * D_201;
  l_b->F_020 += M_100 * D_120 + M_010 * D_030 + M_001 * D_021;
  l_b->F_002 += M_100 * D_102 + M_010 * D_012 + M_001 * D_003;
  l_b->F_110 += M_100 * D_210 + M_010 * D_120 + M_001 * D_111;
  l_b->F_101 += M_100 * D_201 + M_010 * D_111 + M_001 * D_102;
  l_b->F_011 += M_100 * D_111 + M_010 * D_021 + M_001 * D_012;

  /*  3rd order multipole term (addition to rank 3)*/
  l_b->F_300 += M_000 * D_300;
  l_b->F_030 += M_000 * D_030;
  l_b->F_003 += M_000 * D_003;
  l_b->F_210 += M_000 * D_210;
  l_b->F_201 += M_000 * D_201;
  l_b->F_120 += M_000 * D_120;
  l_b->F_021 += M_000 * D_021;
  l_b->F_102 += M_000 * D_102;
  l_b->F_012 += M_000 * D_012;
  l_b->F_111 += M_000 * D_111;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 3

  const float M_400 = m_a->M_400;
  const float M_040 = m_a->M_040;
  const float M_004 = m_a->M_004;
  const float M_310 = m_a->M_310;
  const float M_301 = m_a->M_301;
  const float M_031 = m_a->M_031;
  const float M_130 = m_a->M_130;
  const float M_013 = m_a->M_013;
  const float M_103 = m_a->M_103;
  const float M_220 = m_a->M_220;
  const float M_202 = m_a->M_202;
  const float M_022 = m_a->M_022;
  const float M_211 = m_a->M_211;
  const float M_121 = m_a->M_121;
  const float M_112 = m_a->M_112;

  const float D_400 = pot->D_400;
  const float D_040 = pot->D_040;
  const float D_004 = pot->D_004;
  const float D_310 = pot->D_310;
  const float D_301 = pot->D_301;
  const float D_031 = pot->D_031;
  const float D_130 = pot->D_130;
  const float D_013 = pot->D_013;
  const float D_103 = pot->D_103;
  const float D_220 = pot->D_220;
  const float D_202 = pot->D_202;
  const float D_022 = pot->D_022;
  const float D_211 = pot->D_211;
  const float D_121 = pot->D_121;
  const float D_112 = pot->D_112;

  /* Compute 4th order field tensor terms (addition to rank 0) */
  l_b->F_000 += M_004 * D_004 + M_013 * D_013 + M_022 * D_022 + M_031 * D_031 +
                M_040 * D_040 + M_103 * D_103 + M_112 * D_112 + M_121 * D_121 +
                M_130 * D_130 + M_202 * D_202 + M_211 * D_211 + M_220 * D_220 +
                M_301 * D_301 + M_310 * D_310 + M_400 * D_400;

  /* Compute 4th order field tensor terms (addition to rank 1) */
  l_b->F_001 += M_003 * D_004 + M_012 * D_013 + M_021 * D_022 + M_030 * D_031 +
                M_102 * D_103 + M_111 * D_112 + M_120 * D_121 + M_201 * D_202 +
                M_210 * D_211 + M_300 * D_301;
  l_b->F_010 += M_003 * D_013 + M_012 * D_022 + M_021 * D_031 + M_030 * D_040 +
                M_102 * D_112 + M_111 * D_121 + M_120 * D_130 + M_201 * D_211 +
                M_210 * D_220 + M_300 * D_310;
  l_b->F_100 += M_003 * D_103 + M_012 * D_112 + M_021 * D_121 + M_030 * D_130 +
                M_102 * D_202 + M_111 * D_211 + M_120 * D_220 + M_201 * D_301 +
                M_210 * D_310 + M_300 * D_400;

  /* Compute 4th order field tensor terms (addition to rank 2) */
  l_b->F_002 += M_002 * D_004 + M_011 * D_013 + M_020 * D_022 + M_101 * D_103 +
                M_110 * D_112 + M_200 * D_202;
  l_b->F_011 += M_002 * D_013 + M_011 * D_022 + M_020 * D_031 + M_101 * D_112 +
                M_110 * D_121 + M_200 * D_211;
  l_b->F_020 += M_002 * D_022 + M_011 * D_031 + M_020 * D_040 + M_101 * D_121 +
                M_110 * D_130 + M_200 * D_220;
  l_b->F_101 += M_002 * D_103 + M_011 * D_112 + M_020 * D_121 + M_101 * D_202 +
                M_110 * D_211 + M_200 * D_301;
  l_b->F_110 += M_002 * D_112 + M_011 * D_121 + M_020 * D_130 + M_101 * D_211 +
                M_110 * D_220 + M_200 * D_310;
  l_b->F_200 += M_002 * D_202 + M_011 * D_211 + M_020 * D_220 + M_101 * D_301 +
                M_110 * D_310 + M_200 * D_400;

  /* Compute 4th order field tensor terms (addition to rank 3) */
  l_b->F_003 += M_001 * D_004 + M_010 * D_013 + M_100 * D_103;
  l_b->F_012 += M_001 * D_013 + M_010 * D_022 + M_100 * D_112;
  l_b->F_021 += M_001 * D_022 + M_010 * D_031 + M_100 * D_121;
  l_b->F_030 += M_001 * D_031 + M_010 * D_040 + M_100 * D_130;
  l_b->F_102 += M_001 * D_103 + M_010 * D_112 + M_100 * D_202;
  l_b->F_111 += M_001 * D_112 + M_010 * D_121 + M_100 * D_211;
  l_b->F_120 += M_001 * D_121 + M_010 * D_130 + M_100 * D_220;
  l_b->F_201 += M_001 * D_202 + M_010 * D_211 + M_100 * D_301;
  l_b->F_210 += M_001 * D_211 + M_010 * D_220 + M_100 * D_310;
  l_b->F_300 += M_001 * D_301 + M_010 * D_310 + M_100 * D_400;

  /* Compute 4th order field tensor terms (addition to rank 4) */
  l_b->F_004 += M_000 * D_004;
  l_b->F_013 += M_000 * D_013;
  l_b->F_022 += M_000 * D_022;
  l_b->F_031 += M_000 * D_031;
  l_b->F_040 += M_000 * D_040;
  l_b->F_103 += M_000 * D_103;
  l_b->F_112 += M_000 * D_112;
  l_b->F_121 += M_000 * D_121;
  l_b->F_130 += M_000 * D_130;
  l_b->F_202 += M_000 * D_202;
  l_b->F_211 += M_000 * D_211;
  l_b->F_220 += M_000 * D_220;
  l_b->F_301 += M_000 * D_301;
  l_b->F_310 += M_000 * D_310;
  l_b->F_400 += M_000 * D_400;

#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 4

  const float M_500 = m_a->M_500;
  const float M_050 = m_a->M_050;
  const float M_005 = m_a->M_005;
  const float M_410 = m_a->M_410;
  const float M_401 = m_a->M_401;
  const float M_041 = m_a->M_041;
  const float M_140 = m_a->M_140;
  const float M_014 = m_a->M_014;
  const float M_104 = m_a->M_104;
  const float M_320 = m_a->M_320;
  const float M_302 = m_a->M_302;
  const float M_230 = m_a->M_230;
  const float M_032 = m_a->M_032;
  const float M_203 = m_a->M_203;
  const float M_023 = m_a->M_023;
  const float M_122 = m_a->M_122;
  const float M_212 = m_a->M_212;
  const float M_221 = m_a->M_221;
  const float M_311 = m_a->M_311;
  const float M_131 = m_a->M_131;
  const float M_113 = m_a->M_113;

  const float D_500 = pot->D_500;
  const float D_050 = pot->D_050;
  const float D_005 = pot->D_005;
  const float D_410 = pot->D_410;
  const float D_401 = pot->D_401;
  const float D_041 = pot->D_041;
  const float D_140 = pot->D_140;
  const float D_014 = pot->D_014;
  const float D_104 = pot->D_104;
  const float D_320 = pot->D_320;
  const float D_302 = pot->D_302;
  const float D_230 = pot->D_230;
  const float D_032 = pot->D_032;
  const float D_203 = pot->D_203;
  const float D_023 = pot->D_023;
  const float D_122 = pot->D_122;
  const float D_212 = pot->D_212;
  const float D_221 = pot->D_221;
  const float D_311 = pot->D_311;
  const float D_131 = pot->D_131;
  const float D_113 = pot->D_113;

  /* Compute 5th order field tensor terms (addition to rank 0) */
  l_b->F_000 += M_005 * D_005 + M_014 * D_014 + M_023 * D_023 + M_032 * D_032 +
                M_041 * D_041 + M_050 * D_050 + M_104 * D_104 + M_113 * D_113 +
                M_122 * D_122 + M_131 * D_131 + M_140 * D_140 + M_203 * D_203 +
                M_212 * D_212 + M_221 * D_221 + M_230 * D_230 + M_302 * D_302 +
                M_311 * D_311 + M_320 * D_320 + M_401 * D_401 + M_410 * D_410 +
                M_500 * D_500;

  /* Compute 5th order field tensor terms (addition to rank 1) */
  l_b->F_001 += M_004 * D_005 + M_013 * D_014 + M_022 * D_023 + M_031 * D_032 +
                M_040 * D_041 + M_103 * D_104 + M_112 * D_113 + M_121 * D_122 +
                M_130 * D_131 + M_202 * D_203 + M_211 * D_212 + M_220 * D_221 +
                M_301 * D_302 + M_310 * D_311 + M_400 * D_401;
  l_b->F_010 += M_004 * D_014 + M_013 * D_023 + M_022 * D_032 + M_031 * D_041 +
                M_040 * D_050 + M_103 * D_113 + M_112 * D_122 + M_121 * D_131 +
                M_130 * D_140 + M_202 * D_212 + M_211 * D_221 + M_220 * D_230 +
                M_301 * D_311 + M_310 * D_320 + M_400 * D_410;
  l_b->F_100 += M_004 * D_104 + M_013 * D_113 + M_022 * D_122 + M_031 * D_131 +
                M_040 * D_140 + M_103 * D_203 + M_112 * D_212 + M_121 * D_221 +
                M_130 * D_230 + M_202 * D_302 + M_211 * D_311 + M_220 * D_320 +
                M_301 * D_401 + M_310 * D_410 + M_400 * D_500;

  /* Compute 5th order field tensor terms (addition to rank 2) */
  l_b->F_002 += M_003 * D_005 + M_012 * D_014 + M_021 * D_023 + M_030 * D_032 +
                M_102 * D_104 + M_111 * D_113 + M_120 * D_122 + M_201 * D_203 +
                M_210 * D_212 + M_300 * D_302;
  l_b->F_011 += M_003 * D_014 + M_012 * D_023 + M_021 * D_032 + M_030 * D_041 +
                M_102 * D_113 + M_111 * D_122 + M_120 * D_131 + M_201 * D_212 +
                M_210 * D_221 + M_300 * D_311;
  l_b->F_020 += M_003 * D_023 + M_012 * D_032 + M_021 * D_041 + M_030 * D_050 +
                M_102 * D_122 + M_111 * D_131 + M_120 * D_140 + M_201 * D_221 +
                M_210 * D_230 + M_300 * D_320;
  l_b->F_101 += M_003 * D_104 + M_012 * D_113 + M_021 * D_122 + M_030 * D_131 +
                M_102 * D_203 + M_111 * D_212 + M_120 * D_221 + M_201 * D_302 +
                M_210 * D_311 + M_300 * D_401;
  l_b->F_110 += M_003 * D_113 + M_012 * D_122 + M_021 * D_131 + M_030 * D_140 +
                M_102 * D_212 + M_111 * D_221 + M_120 * D_230 + M_201 * D_311 +
                M_210 * D_320 + M_300 * D_410;
  l_b->F_200 += M_003 * D_203 + M_012 * D_212 + M_021 * D_221 + M_030 * D_230 +
                M_102 * D_302 + M_111 * D_311 + M_120 * D_320 + M_201 * D_401 +
                M_210 * D_410 + M_300 * D_500;

  /* Compute 5th order field tensor terms (addition to rank 3) */
  l_b->F_003 += M_002 * D_005 + M_011 * D_014 + M_020 * D_023 + M_101 * D_104 +
                M_110 * D_113 + M_200 * D_203;
  l_b->F_012 += M_002 * D_014 + M_011 * D_023 + M_020 * D_032 + M_101 * D_113 +
                M_110 * D_122 + M_200 * D_212;
  l_b->F_021 += M_002 * D_023 + M_011 * D_032 + M_020 * D_041 + M_101 * D_122 +
                M_110 * D_131 + M_200 * D_221;
  l_b->F_030 += M_002 * D_032 + M_011 * D_041 + M_020 * D_050 + M_101 * D_131 +
                M_110 * D_140 + M_200 * D_230;
  l_b->F_102 += M_002 * D_104 + M_011 * D_113 + M_020 * D_122 + M_101 * D_203 +
                M_110 * D_212 + M_200 * D_302;
  l_b->F_111 += M_002 * D_113 + M_011 * D_122 + M_020 * D_131 + M_101 * D_212 +
                M_110 * D_221 + M_200 * D_311;
  l_b->F_120 += M_002 * D_122 + M_011 * D_131 + M_020 * D_140 + M_101 * D_221 +
                M_110 * D_230 + M_200 * D_320;
  l_b->F_201 += M_002 * D_203 + M_011 * D_212 + M_020 * D_221 + M_101 * D_302 +
                M_110 * D_311 + M_200 * D_401;
  l_b->F_210 += M_002 * D_212 + M_011 * D_221 + M_020 * D_230 + M_101 * D_311 +
                M_110 * D_320 + M_200 * D_410;
  l_b->F_300 += M_002 * D_302 + M_011 * D_311 + M_020 * D_320 + M_101 * D_401 +
                M_110 * D_410 + M_200 * D_500;

  /* Compute 5th order field tensor terms (addition to rank 4) */
  l_b->F_004 += M_001 * D_005 + M_010 * D_014 + M_100 * D_104;
  l_b->F_013 += M_001 * D_014 + M_010 * D_023 + M_100 * D_113;
  l_b->F_022 += M_001 * D_023 + M_010 * D_032 + M_100 * D_122;
  l_b->F_031 += M_001 * D_032 + M_010 * D_041 + M_100 * D_131;
  l_b->F_040 += M_001 * D_041 + M_010 * D_050 + M_100 * D_140;
  l_b->F_103 += M_001 * D_104 + M_010 * D_113 + M_100 * D_203;
  l_b->F_112 += M_001 * D_113 + M_010 * D_122 + M_100 * D_212;
  l_b->F_121 += M_001 * D_122 + M_010 * D_131 + M_100 * D_221;
  l_b->F_130 += M_001 * D_131 + M_010 * D_140 + M_100 * D_230;
  l_b->F_202 += M_001 * D_203 + M_010 * D_212 + M_100 * D_302;
  l_b->F_211 += M_001 * D_212 + M_010 * D_221 + M_100 * D_311;
  l_b->F_220 += M_001 * D_221 + M_010 * D_230 + M_100 * D_320;
  l_b->F_301 += M_001 * D_302 + M_010 * D_311 + M_100 * D_401;
  l_b->F_310 += M_001 * D_311 + M_010 * D_320 + M_100 * D_410;
  l_b->F_400 += M_001 * D_401 + M_010 * D_410 + M_100 * D_500;

  /* Compute 5th order field tensor terms (addition to rank 5) */
  l_b->F_005 += M_000 * D_005;
  l_b->F_014 += M_000 * D_014;
  l_b->F_023 += M_000 * D_023;
  l_b->F_032 += M_000 * D_032;
  l_b->F_041 += M_000 * D_041;
  l_b->F_050 += M_000 * D_050;
  l_b->F_104 += M_000 * D_104;
  l_b->F_113 += M_000 * D_113;
  l_b->F_122 += M_000 * D_122;
  l_b->F_131 += M_000 * D_131;
  l_b->F_140 += M_000 * D_140;
  l_b->F_203 += M_000 * D_203;
  l_b->F_212 += M_000 * D_212;
  l_b->F_221 += M_000 * D_221;
  l_b->F_230 += M_000 * D_230;
  l_b->F_302 += M_000 * D_302;
  l_b->F_311 += M_000 * D_311;
  l_b->F_320 += M_000 * D_320;
  l_b->F_401 += M_000 * D_401;
  l_b->F_410 += M_000 * D_410;
  l_b->F_500 += M_000 * D_500;

#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 5

  const float M_006 = m_a->M_006;
  const float M_015 = m_a->M_015;
  const float M_024 = m_a->M_024;
  const float M_033 = m_a->M_033;
  const float M_042 = m_a->M_042;
  const float M_051 = m_a->M_051;
  const float M_060 = m_a->M_060;
  const float M_105 = m_a->M_105;
  const float M_114 = m_a->M_114;
  const float M_123 = m_a->M_123;
  const float M_132 = m_a->M_132;
  const float M_141 = m_a->M_141;
  const float M_150 = m_a->M_150;
  const float M_204 = m_a->M_204;
  const float M_213 = m_a->M_213;
  const float M_222 = m_a->M_222;
  const float M_231 = m_a->M_231;
  const float M_240 = m_a->M_240;
  const float M_303 = m_a->M_303;
  const float M_312 = m_a->M_312;
  const float M_321 = m_a->M_321;
  const float M_330 = m_a->M_330;
  const float M_402 = m_a->M_402;
  const float M_411 = m_a->M_411;
  const float M_420 = m_a->M_420;
  const float M_501 = m_a->M_501;
  const float M_510 = m_a->M_510;
  const float M_600 = m_a->M_600;

  const float D_006 = pot->D_006;
  const float D_015 = pot->D_015;
  const float D_024 = pot->D_024;
  const float D_033 = pot->D_033;
  const float D_042 = pot->D_042;
  const float D_051 = pot->D_051;
  const float D_060 = pot->D_060;
  const float D_105 = pot->D_105;
  const float D_114 = pot->D_114;
  const float D_123 = pot->D_123;
  const float D_132 = pot->D_132;
  const float D_141 = pot->D_141;
  const float D_150 = pot->D_150;
  const float D_204 = pot->D_204;
  const float D_213 = pot->D_213;
  const float D_222 = pot->D_222;
  const float D_231 = pot->D_231;
  const float D_240 = pot->D_240;
  const float D_303 = pot->D_303;
  const float D_312 = pot->D_312;
  const float D_321 = pot->D_321;
  const float D_330 = pot->D_330;
  const float D_402 = pot->D_402;
  const float D_411 = pot->D_411;
  const float D_420 = pot->D_420;
  const float D_501 = pot->D_501;
  const float D_510 = pot->D_510;
  const float D_600 = pot->D_600;

  /* Compute 6th order field tensor terms (addition to rank 0) */
  l_b->F_000 += M_006 * D_006 + M_015 * D_015 + M_024 * D_024
    + M_033 * D_033 + M_042 * D_042 + M_051 * D_051
    + M_060 * D_060 + M_105 * D_105 + M_114 * D_114
    + M_123 * D_123 + M_132 * D_132 + M_141 * D_141
    + M_150 * D_150 + M_204 * D_204 + M_213 * D_213
    + M_222 * D_222 + M_231 * D_231 + M_240 * D_240
    + M_303 * D_303 + M_312 * D_312 + M_321 * D_321
    + M_330 * D_330 + M_402 * D_402 + M_411 * D_411
    + M_420 * D_420 + M_501 * D_501 + M_510 * D_510
    + M_600 * D_600;

  /* Compute 6th order field tensor terms (addition to rank 1) */
  l_b->F_001 += M_005 * D_006 + M_014 * D_015 + M_023 * D_024
    + M_032 * D_033 + M_041 * D_042 + M_050 * D_051
    + M_104 * D_105 + M_113 * D_114 + M_122 * D_123
    + M_131 * D_132 + M_140 * D_141 + M_203 * D_204
    + M_212 * D_213 + M_221 * D_222 + M_230 * D_231
    + M_302 * D_303 + M_311 * D_312 + M_320 * D_321
    + M_401 * D_402 + M_410 * D_411 + M_500 * D_501;
  l_b->F_010 += M_005 * D_015 + M_014 * D_024 + M_023 * D_033
    + M_032 * D_042 + M_041 * D_051 + M_050 * D_060
    + M_104 * D_114 + M_113 * D_123 + M_122 * D_132
    + M_131 * D_141 + M_140 * D_150 + M_203 * D_213
    + M_212 * D_222 + M_221 * D_231 + M_230 * D_240
    + M_302 * D_312 + M_311 * D_321 + M_320 * D_330
    + M_401 * D_411 + M_410 * D_420 + M_500 * D_510;
  l_b->F_100 += M_005 * D_105 + M_014 * D_114 + M_023 * D_123
    + M_032 * D_132 + M_041 * D_141 + M_050 * D_150
    + M_104 * D_204 + M_113 * D_213 + M_122 * D_222
    + M_131 * D_231 + M_140 * D_240 + M_203 * D_303
    + M_212 * D_312 + M_221 * D_321 + M_230 * D_330
    + M_302 * D_402 + M_311 * D_411 + M_320 * D_420
    + M_401 * D_501 + M_410 * D_510 + M_500 * D_600;

  /* Compute 6th order field tensor terms (addition to rank 2) */
  l_b->F_002 += M_004 * D_006 + M_013 * D_015 + M_022 * D_024
    + M_031 * D_033 + M_040 * D_042 + M_103 * D_105
    + M_112 * D_114 + M_121 * D_123 + M_130 * D_132
    + M_202 * D_204 + M_211 * D_213 + M_220 * D_222
    + M_301 * D_303 + M_310 * D_312 + M_400 * D_402;
  l_b->F_011 += M_004 * D_015 + M_013 * D_024 + M_022 * D_033
    + M_031 * D_042 + M_040 * D_051 + M_103 * D_114
    + M_112 * D_123 + M_121 * D_132 + M_130 * D_141
    + M_202 * D_213 + M_211 * D_222 + M_220 * D_231
    + M_301 * D_312 + M_310 * D_321 + M_400 * D_411;
  l_b->F_020 += M_004 * D_024 + M_013 * D_033 + M_022 * D_042
    + M_031 * D_051 + M_040 * D_060 + M_103 * D_123
    + M_112 * D_132 + M_121 * D_141 + M_130 * D_150
    + M_202 * D_222 + M_211 * D_231 + M_220 * D_240
    + M_301 * D_321 + M_310 * D_330 + M_400 * D_420;
  l_b->F_101 += M_004 * D_105 + M_013 * D_114 + M_022 * D_123
    + M_031 * D_132 + M_040 * D_141 + M_103 * D_204
    + M_112 * D_213 + M_121 * D_222 + M_130 * D_231
    + M_202 * D_303 + M_211 * D_312 + M_220 * D_321
    + M_301 * D_402 + M_310 * D_411 + M_400 * D_501;
  l_b->F_110 += M_004 * D_114 + M_013 * D_123 + M_022 * D_132
    + M_031 * D_141 + M_040 * D_150 + M_103 * D_213
    + M_112 * D_222 + M_121 * D_231 + M_130 * D_240
    + M_202 * D_312 + M_211 * D_321 + M_220 * D_330
    + M_301 * D_411 + M_310 * D_420 + M_400 * D_510;
  l_b->F_200 += M_004 * D_204 + M_013 * D_213 + M_022 * D_222
    + M_031 * D_231 + M_040 * D_240 + M_103 * D_303
    + M_112 * D_312 + M_121 * D_321 + M_130 * D_330
    + M_202 * D_402 + M_211 * D_411 + M_220 * D_420
    + M_301 * D_501 + M_310 * D_510 + M_400 * D_600;

  /* Compute 6th order field tensor terms (addition to rank 3) */
  l_b->F_003 += M_003 * D_006 + M_012 * D_015 + M_021 * D_024
    + M_030 * D_033 + M_102 * D_105 + M_111 * D_114
    + M_120 * D_123 + M_201 * D_204 + M_210 * D_213
    + M_300 * D_303;
  l_b->F_012 += M_003 * D_015 + M_012 * D_024 + M_021 * D_033
    + M_030 * D_042 + M_102 * D_114 + M_111 * D_123
    + M_120 * D_132 + M_201 * D_213 + M_210 * D_222
    + M_300 * D_312;
  l_b->F_021 += M_003 * D_024 + M_012 * D_033 + M_021 * D_042
    + M_030 * D_051 + M_102 * D_123 + M_111 * D_132
    + M_120 * D_141 + M_201 * D_222 + M_210 * D_231
    + M_300 * D_321;
  l_b->F_030 += M_003 * D_033 + M_012 * D_042 + M_021 * D_051
    + M_030 * D_060 + M_102 * D_132 + M_111 * D_141
    + M_120 * D_150 + M_201 * D_231 + M_210 * D_240
    + M_300 * D_330;
  l_b->F_102 += M_003 * D_105 + M_012 * D_114 + M_021 * D_123
    + M_030 * D_132 + M_102 * D_204 + M_111 * D_213
    + M_120 * D_222 + M_201 * D_303 + M_210 * D_312
    + M_300 * D_402;
  l_b->F_111 += M_003 * D_114 + M_012 * D_123 + M_021 * D_132
    + M_030 * D_141 + M_102 * D_213 + M_111 * D_222
    + M_120 * D_231 + M_201 * D_312 + M_210 * D_321
    + M_300 * D_411;
  l_b->F_120 += M_003 * D_123 + M_012 * D_132 + M_021 * D_141
    + M_030 * D_150 + M_102 * D_222 + M_111 * D_231
    + M_120 * D_240 + M_201 * D_321 + M_210 * D_330
    + M_300 * D_420;
  l_b->F_201 += M_003 * D_204 + M_012 * D_213 + M_021 * D_222
    + M_030 * D_231 + M_102 * D_303 + M_111 * D_312
    + M_120 * D_321 + M_201 * D_402 + M_210 * D_411
    + M_300 * D_501;
  l_b->F_210 += M_003 * D_213 + M_012 * D_222 + M_021 * D_231
    + M_030 * D_240 + M_102 * D_312 + M_111 * D_321
    + M_120 * D_330 + M_201 * D_411 + M_210 * D_420
    + M_300 * D_510;
  l_b->F_300 += M_003 * D_303 + M_012 * D_312 + M_021 * D_321
    + M_030 * D_330 + M_102 * D_402 + M_111 * D_411
    + M_120 * D_420 + M_201 * D_501 + M_210 * D_510
    + M_300 * D_600;

  /* Compute 6th order field tensor terms (addition to rank 4) */
  l_b->F_004 += M_002 * D_006 + M_011 * D_015 + M_020 * D_024
    + M_101 * D_105 + M_110 * D_114 + M_200 * D_204;
  l_b->F_013 += M_002 * D_015 + M_011 * D_024 + M_020 * D_033
    + M_101 * D_114 + M_110 * D_123 + M_200 * D_213;
  l_b->F_022 += M_002 * D_024 + M_011 * D_033 + M_020 * D_042
    + M_101 * D_123 + M_110 * D_132 + M_200 * D_222;
  l_b->F_031 += M_002 * D_033 + M_011 * D_042 + M_020 * D_051
    + M_101 * D_132 + M_110 * D_141 + M_200 * D_231;
  l_b->F_040 += M_002 * D_042 + M_011 * D_051 + M_020 * D_060
    + M_101 * D_141 + M_110 * D_150 + M_200 * D_240;
  l_b->F_103 += M_002 * D_105 + M_011 * D_114 + M_020 * D_123
    + M_101 * D_204 + M_110 * D_213 + M_200 * D_303;
  l_b->F_112 += M_002 * D_114 + M_011 * D_123 + M_020 * D_132
    + M_101 * D_213 + M_110 * D_222 + M_200 * D_312;
  l_b->F_121 += M_002 * D_123 + M_011 * D_132 + M_020 * D_141
    + M_101 * D_222 + M_110 * D_231 + M_200 * D_321;
  l_b->F_130 += M_002 * D_132 + M_011 * D_141 + M_020 * D_150
    + M_101 * D_231 + M_110 * D_240 + M_200 * D_330;
  l_b->F_202 += M_002 * D_204 + M_011 * D_213 + M_020 * D_222
    + M_101 * D_303 + M_110 * D_312 + M_200 * D_402;
  l_b->F_211 += M_002 * D_213 + M_011 * D_222 + M_020 * D_231
    + M_101 * D_312 + M_110 * D_321 + M_200 * D_411;
  l_b->F_220 += M_002 * D_222 + M_011 * D_231 + M_020 * D_240
    + M_101 * D_321 + M_110 * D_330 + M_200 * D_420;
  l_b->F_301 += M_002 * D_303 + M_011 * D_312 + M_020 * D_321
    + M_101 * D_402 + M_110 * D_411 + M_200 * D_501;
  l_b->F_310 += M_002 * D_312 + M_011 * D_321 + M_020 * D_330
    + M_101 * D_411 + M_110 * D_420 + M_200 * D_510;
  l_b->F_400 += M_002 * D_402 + M_011 * D_411 + M_020 * D_420
    + M_101 * D_501 + M_110 * D_510 + M_200 * D_600;

  /* Compute 6th order field tensor terms (addition to rank 5) */
  l_b->F_005 += M_001 * D_006 + M_010 * D_015 + M_100 * D_105;
  l_b->F_014 += M_001 * D_015 + M_010 * D_024 + M_100 * D_114;
  l_b->F_023 += M_001 * D_024 + M_010 * D_033 + M_100 * D_123;
  l_b->F_032 += M_001 * D_033 + M_010 * D_042 + M_100 * D_132;
  l_b->F_041 += M_001 * D_042 + M_010 * D_051 + M_100 * D_141;
  l_b->F_050 += M_001 * D_051 + M_010 * D_060 + M_100 * D_150;
  l_b->F_104 += M_001 * D_105 + M_010 * D_114 + M_100 * D_204;
  l_b->F_113 += M_001 * D_114 + M_010 * D_123 + M_100 * D_213;
  l_b->F_122 += M_001 * D_123 + M_010 * D_132 + M_100 * D_222;
  l_b->F_131 += M_001 * D_132 + M_010 * D_141 + M_100 * D_231;
  l_b->F_140 += M_001 * D_141 + M_010 * D_150 + M_100 * D_240;
  l_b->F_203 += M_001 * D_204 + M_010 * D_213 + M_100 * D_303;
  l_b->F_212 += M_001 * D_213 + M_010 * D_222 + M_100 * D_312;
  l_b->F_221 += M_001 * D_222 + M_010 * D_231 + M_100 * D_321;
  l_b->F_230 += M_001 * D_231 + M_010 * D_240 + M_100 * D_330;
  l_b->F_302 += M_001 * D_303 + M_010 * D_312 + M_100 * D_402;
  l_b->F_311 += M_001 * D_312 + M_010 * D_321 + M_100 * D_411;
  l_b->F_320 += M_001 * D_321 + M_010 * D_330 + M_100 * D_420;
  l_b->F_401 += M_001 * D_402 + M_010 * D_411 + M_100 * D_501;
  l_b->F_410 += M_001 * D_411 + M_010 * D_420 + M_100 * D_510;
  l_b->F_500 += M_001 * D_501 + M_010 * D_510 + M_100 * D_600;

  /* Compute 6th order field tensor terms (addition to rank 6) */
  l_b->F_006 += M_000 * D_006;
  l_b->F_015 += M_000 * D_015;
  l_b->F_024 += M_000 * D_024;
  l_b->F_033 += M_000 * D_033;
  l_b->F_042 += M_000 * D_042;
  l_b->F_051 += M_000 * D_051;
  l_b->F_060 += M_000 * D_060;
  l_b->F_105 += M_000 * D_105;
  l_b->F_114 += M_000 * D_114;
  l_b->F_123 += M_000 * D_123;
  l_b->F_132 += M_000 * D_132;
  l_b->F_141 += M_000 * D_141;
  l_b->F_150 += M_000 * D_150;
  l_b->F_204 += M_000 * D_204;
  l_b->F_213 += M_000 * D_213;
  l_b->F_222 += M_000 * D_222;
  l_b->F_231 += M_000 * D_231;
  l_b->F_240 += M_000 * D_240;
  l_b->F_303 += M_000 * D_303;
  l_b->F_312 += M_000 * D_312;
  l_b->F_321 += M_000 * D_321;
  l_b->F_330 += M_000 * D_330;
  l_b->F_402 += M_000 * D_402;
  l_b->F_411 += M_000 * D_411;
  l_b->F_420 += M_000 * D_420;
  l_b->F_501 += M_000 * D_501;
  l_b->F_510 += M_000 * D_510;
  l_b->F_600 += M_000 * D_600;

#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 6
#error "Missing implementation for order >6"
#endif
}

/**
 * @brief Compute the field tensor due to a multipole.
 *
 * @param l_b The field tensor to compute.
 * @param m_a The multipole.
 * @param pos_b The position of the field tensor.
 * @param pos_a The position of the multipole.
 * @param props The #gravity_props of this calculation.
 * @param periodic Is the calculation periodic ?
 * @param dim The size of the simulation box.
 * @param rs_inv The inverse of the gravity mesh-smoothing scale.
 */
INLINE static void gravity_M2L_nonsym(
    struct grav_tensor *l_b, const struct multipole *m_a, const double pos_b[3],
    const double pos_a[3], const struct gravity_props *props,
    const int periodic, const double dim[3], const float rs_inv) {

  /* Recover some constants */
  const float eps = m_a->max_softening;
  const float eps_inv = 1.f / eps;

  /* Compute distance vector */
  float dx = (float)(pos_b[0] - pos_a[0]);
  float dy = (float)(pos_b[1] - pos_a[1]);
  float dz = (float)(pos_b[2] - pos_a[2]);

  /* Apply BC */
  if (periodic) {
    dx = nearest(dx, dim[0]);
    dy = nearest(dy, dim[1]);
    dz = nearest(dz, dim[2]);
  }

  /* Compute distance */
  const float r2 = dx * dx + dy * dy + dz * dz;
  const float r_inv = 1. / sqrtf(r2);

  /* Compute all derivatives */
  struct potential_derivatives_M2L pot;
  potential_derivatives_compute_M2L(dx, dy, dz, r2, r_inv, eps, eps_inv,
                                    periodic, rs_inv, &pot);

  /* Do the M2L tensor multiplication */
  gravity_M2L_apply(l_b, m_a, &pot);
}

/**
 * @brief Compute the field tensor due to a multipole and the symmetric
 * equivalent.
 *
 * @param l_a The first field tensor to compute.
 * @param l_b The second field tensor to compute.
 * @param m_a The first multipole.
 * @param m_b The second multipole.
 * @param pos_a The position of the first m-pole and field tensor.
 * @param pos_b The position of the second m-pole and field tensor.
 * @param props The #gravity_props of this calculation.
 * @param periodic Is the calculation periodic ?
 * @param dim The size of the simulation box.
 * @param rs_inv The inverse of the gravity mesh-smoothing scale.
 */
INLINE static void gravity_M2L_symmetric(
    struct grav_tensor *restrict l_a, struct grav_tensor *restrict l_b,
    const struct multipole *restrict m_a, const struct multipole *restrict m_b,
    const double pos_a[3], const double pos_b[3],
    const struct gravity_props *props, const int periodic, const double dim[3],
    const float rs_inv) {

  /* Recover some constants */
  const float eps = m_a->max_softening;
  const float eps_inv = 1.f / eps;

  /* Compute distance vector */
  float dx = (float)(pos_b[0] - pos_a[0]);
  float dy = (float)(pos_b[1] - pos_a[1]);
  float dz = (float)(pos_b[2] - pos_a[2]);

  /* Apply BC */
  if (periodic) {
    dx = nearest(dx, dim[0]);
    dy = nearest(dy, dim[1]);
    dz = nearest(dz, dim[2]);
  }

  /* Compute distance */
  const float r2 = dx * dx + dy * dy + dz * dz;
  const float r_inv = 1. / sqrtf(r2);

  /* Compute all derivatives */
  struct potential_derivatives_M2L pot;
  potential_derivatives_compute_M2L(dx, dy, dz, r2, r_inv, eps, eps_inv,
                                    periodic, rs_inv, &pot);

  /* Do the first M2L tensor multiplication */
  gravity_M2L_apply(l_b, m_a, &pot);

  /* Flip the signs of odd derivatives */
  potential_derivatives_flip_signs(&pot);

  /* Do the second M2L tensor multiplication */
  gravity_M2L_apply(l_a, m_b, &pot);
}

/**
 * @brief Creates a copy of #grav_tensor shifted to a new location.
 *
 * Corresponds to equation (28e).
 *
 * @param la The #grav_tensor copy (content will  be overwritten).
 * @param lb The #grav_tensor to shift.
 * @param pos_a The position to which m_b will be shifted.
 * @param pos_b The current postion of the multipole to shift.
 */
INLINE static void gravity_L2L(struct grav_tensor *restrict la,
                               const struct grav_tensor *restrict lb,
                               const double pos_a[3], const double pos_b[3]) {

  /* Initialise everything to zero */
  gravity_field_tensors_init(la, 0);

#ifdef SWIFT_DEBUG_CHECKS
  if (lb->num_interacted + lb->num_not_interacted == 0)
      error("Shifting tensors that did not interact");
#endif
#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_GRAVITY_FORCE_CHECKS)
  la->num_interacted = lb->num_interacted;
  la->num_not_interacted = lb->num_not_interacted;
#endif

  /* Distance to shift by */
  const double dx[3] = {pos_a[0] - pos_b[0], pos_a[1] - pos_b[1],
                        pos_a[2] - pos_b[2]};

  /* Shift 0th order term */
  la->F_000 += X_000(dx) * lb->F_000;

#if SELF_GRAVITY_MULTIPOLE_ORDER > 0
  /* Shift 1st order multipole term (addition to rank 0)*/
  la->F_000 +=
      X_100(dx) * lb->F_100 + X_010(dx) * lb->F_010 + X_001(dx) * lb->F_001;

  /* Shift 1st order multipole term (addition to rank 1)*/
  la->F_100 += X_000(dx) * lb->F_100;
  la->F_010 += X_000(dx) * lb->F_010;
  la->F_001 += X_000(dx) * lb->F_001;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 1

  /* Shift 2nd order multipole term (addition to rank 0)*/
  la->F_000 +=
      X_200(dx) * lb->F_200 + X_020(dx) * lb->F_020 + X_002(dx) * lb->F_002;
  la->F_000 +=
      X_110(dx) * lb->F_110 + X_101(dx) * lb->F_101 + X_011(dx) * lb->F_011;

  /* Shift 2nd order multipole term (addition to rank 1)*/
  la->F_100 +=
      X_100(dx) * lb->F_200 + X_010(dx) * lb->F_110 + X_001(dx) * lb->F_101;
  la->F_010 +=
      X_100(dx) * lb->F_110 + X_010(dx) * lb->F_020 + X_001(dx) * lb->F_011;
  la->F_001 +=
      X_100(dx) * lb->F_101 + X_010(dx) * lb->F_011 + X_001(dx) * lb->F_002;

  /* Shift 2nd order multipole term (addition to rank 2)*/
  la->F_200 += X_000(dx) * lb->F_200;
  la->F_020 += X_000(dx) * lb->F_020;
  la->F_002 += X_000(dx) * lb->F_002;
  la->F_110 += X_000(dx) * lb->F_110;
  la->F_101 += X_000(dx) * lb->F_101;
  la->F_011 += X_000(dx) * lb->F_011;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 2

  /* Shift 3rd order multipole term (addition to rank 0)*/
  la->F_000 +=
      X_300(dx) * lb->F_300 + X_030(dx) * lb->F_030 + X_003(dx) * lb->F_003;
  la->F_000 +=
      X_210(dx) * lb->F_210 + X_201(dx) * lb->F_201 + X_120(dx) * lb->F_120;
  la->F_000 +=
      X_021(dx) * lb->F_021 + X_102(dx) * lb->F_102 + X_012(dx) * lb->F_012;
  la->F_000 += X_111(dx) * lb->F_111;

  /* Shift 3rd order multipole term (addition to rank 1)*/
  la->F_100 +=
      X_200(dx) * lb->F_300 + X_020(dx) * lb->F_120 + X_002(dx) * lb->F_102;
  la->F_100 +=
      X_110(dx) * lb->F_210 + X_101(dx) * lb->F_201 + X_011(dx) * lb->F_111;
  la->F_010 +=
      X_200(dx) * lb->F_210 + X_020(dx) * lb->F_030 + X_002(dx) * lb->F_012;
  la->F_010 +=
      X_110(dx) * lb->F_120 + X_101(dx) * lb->F_111 + X_011(dx) * lb->F_021;
  la->F_001 +=
      X_200(dx) * lb->F_201 + X_020(dx) * lb->F_021 + X_002(dx) * lb->F_003;
  la->F_001 +=
      X_110(dx) * lb->F_111 + X_101(dx) * lb->F_102 + X_011(dx) * lb->F_012;

  /* Shift 3rd order multipole term (addition to rank 2)*/
  la->F_200 +=
      X_100(dx) * lb->F_300 + X_010(dx) * lb->F_210 + X_001(dx) * lb->F_201;
  la->F_020 +=
      X_100(dx) * lb->F_120 + X_010(dx) * lb->F_030 + X_001(dx) * lb->F_021;
  la->F_002 +=
      X_100(dx) * lb->F_102 + X_010(dx) * lb->F_012 + X_001(dx) * lb->F_003;
  la->F_110 +=
      X_100(dx) * lb->F_210 + X_010(dx) * lb->F_120 + X_001(dx) * lb->F_111;
  la->F_101 +=
      X_100(dx) * lb->F_201 + X_010(dx) * lb->F_111 + X_001(dx) * lb->F_102;
  la->F_011 +=
      X_100(dx) * lb->F_111 + X_010(dx) * lb->F_021 + X_001(dx) * lb->F_012;

  /* Shift 3rd order multipole term (addition to rank 2)*/
  la->F_300 += X_000(dx) * lb->F_300;
  la->F_030 += X_000(dx) * lb->F_030;
  la->F_003 += X_000(dx) * lb->F_003;
  la->F_210 += X_000(dx) * lb->F_210;
  la->F_201 += X_000(dx) * lb->F_201;
  la->F_120 += X_000(dx) * lb->F_120;
  la->F_021 += X_000(dx) * lb->F_021;
  la->F_102 += X_000(dx) * lb->F_102;
  la->F_012 += X_000(dx) * lb->F_012;
  la->F_111 += X_000(dx) * lb->F_111;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 3

  /* Shift 4th order field tensor terms (addition to rank 0) */
  la->F_000 +=
      X_004(dx) * lb->F_004 + X_013(dx) * lb->F_013 + X_022(dx) * lb->F_022 +
      X_031(dx) * lb->F_031 + X_040(dx) * lb->F_040 + X_103(dx) * lb->F_103 +
      X_112(dx) * lb->F_112 + X_121(dx) * lb->F_121 + X_130(dx) * lb->F_130 +
      X_202(dx) * lb->F_202 + X_211(dx) * lb->F_211 + X_220(dx) * lb->F_220 +
      X_301(dx) * lb->F_301 + X_310(dx) * lb->F_310 + X_400(dx) * lb->F_400;

  /* Shift 4th order field tensor terms (addition to rank 1) */
  la->F_001 += X_003(dx) * lb->F_004 + X_012(dx) * lb->F_013 +
               X_021(dx) * lb->F_022 + X_030(dx) * lb->F_031 +
               X_102(dx) * lb->F_103 + X_111(dx) * lb->F_112 +
               X_120(dx) * lb->F_121 + X_201(dx) * lb->F_202 +
               X_210(dx) * lb->F_211 + X_300(dx) * lb->F_301;
  la->F_010 += X_003(dx) * lb->F_013 + X_012(dx) * lb->F_022 +
               X_021(dx) * lb->F_031 + X_030(dx) * lb->F_040 +
               X_102(dx) * lb->F_112 + X_111(dx) * lb->F_121 +
               X_120(dx) * lb->F_130 + X_201(dx) * lb->F_211 +
               X_210(dx) * lb->F_220 + X_300(dx) * lb->F_310;
  la->F_100 += X_003(dx) * lb->F_103 + X_012(dx) * lb->F_112 +
               X_021(dx) * lb->F_121 + X_030(dx) * lb->F_130 +
               X_102(dx) * lb->F_202 + X_111(dx) * lb->F_211 +
               X_120(dx) * lb->F_220 + X_201(dx) * lb->F_301 +
               X_210(dx) * lb->F_310 + X_300(dx) * lb->F_400;

  /* Shift 4th order field tensor terms (addition to rank 2) */
  la->F_002 += X_002(dx) * lb->F_004 + X_011(dx) * lb->F_013 +
               X_020(dx) * lb->F_022 + X_101(dx) * lb->F_103 +
               X_110(dx) * lb->F_112 + X_200(dx) * lb->F_202;
  la->F_011 += X_002(dx) * lb->F_013 + X_011(dx) * lb->F_022 +
               X_020(dx) * lb->F_031 + X_101(dx) * lb->F_112 +
               X_110(dx) * lb->F_121 + X_200(dx) * lb->F_211;
  la->F_020 += X_002(dx) * lb->F_022 + X_011(dx) * lb->F_031 +
               X_020(dx) * lb->F_040 + X_101(dx) * lb->F_121 +
               X_110(dx) * lb->F_130 + X_200(dx) * lb->F_220;
  la->F_101 += X_002(dx) * lb->F_103 + X_011(dx) * lb->F_112 +
               X_020(dx) * lb->F_121 + X_101(dx) * lb->F_202 +
               X_110(dx) * lb->F_211 + X_200(dx) * lb->F_301;
  la->F_110 += X_002(dx) * lb->F_112 + X_011(dx) * lb->F_121 +
               X_020(dx) * lb->F_130 + X_101(dx) * lb->F_211 +
               X_110(dx) * lb->F_220 + X_200(dx) * lb->F_310;
  la->F_200 += X_002(dx) * lb->F_202 + X_011(dx) * lb->F_211 +
               X_020(dx) * lb->F_220 + X_101(dx) * lb->F_301 +
               X_110(dx) * lb->F_310 + X_200(dx) * lb->F_400;

  /* Shift 4th order field tensor terms (addition to rank 3) */
  la->F_003 +=
      X_001(dx) * lb->F_004 + X_010(dx) * lb->F_013 + X_100(dx) * lb->F_103;
  la->F_012 +=
      X_001(dx) * lb->F_013 + X_010(dx) * lb->F_022 + X_100(dx) * lb->F_112;
  la->F_021 +=
      X_001(dx) * lb->F_022 + X_010(dx) * lb->F_031 + X_100(dx) * lb->F_121;
  la->F_030 +=
      X_001(dx) * lb->F_031 + X_010(dx) * lb->F_040 + X_100(dx) * lb->F_130;
  la->F_102 +=
      X_001(dx) * lb->F_103 + X_010(dx) * lb->F_112 + X_100(dx) * lb->F_202;
  la->F_111 +=
      X_001(dx) * lb->F_112 + X_010(dx) * lb->F_121 + X_100(dx) * lb->F_211;
  la->F_120 +=
      X_001(dx) * lb->F_121 + X_010(dx) * lb->F_130 + X_100(dx) * lb->F_220;
  la->F_201 +=
      X_001(dx) * lb->F_202 + X_010(dx) * lb->F_211 + X_100(dx) * lb->F_301;
  la->F_210 +=
      X_001(dx) * lb->F_211 + X_010(dx) * lb->F_220 + X_100(dx) * lb->F_310;
  la->F_300 +=
      X_001(dx) * lb->F_301 + X_010(dx) * lb->F_310 + X_100(dx) * lb->F_400;

  /* Shift 4th order field tensor terms (addition to rank 4) */
  la->F_004 += X_000(dx) * lb->F_004;
  la->F_013 += X_000(dx) * lb->F_013;
  la->F_022 += X_000(dx) * lb->F_022;
  la->F_031 += X_000(dx) * lb->F_031;
  la->F_040 += X_000(dx) * lb->F_040;
  la->F_103 += X_000(dx) * lb->F_103;
  la->F_112 += X_000(dx) * lb->F_112;
  la->F_121 += X_000(dx) * lb->F_121;
  la->F_130 += X_000(dx) * lb->F_130;
  la->F_202 += X_000(dx) * lb->F_202;
  la->F_211 += X_000(dx) * lb->F_211;
  la->F_220 += X_000(dx) * lb->F_220;
  la->F_301 += X_000(dx) * lb->F_301;
  la->F_310 += X_000(dx) * lb->F_310;
  la->F_400 += X_000(dx) * lb->F_400;

#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 4

  /* Shift 5th order field tensor terms (addition to rank 0) */
  la->F_000 +=
      X_005(dx) * lb->F_005 + X_014(dx) * lb->F_014 + X_023(dx) * lb->F_023 +
      X_032(dx) * lb->F_032 + X_041(dx) * lb->F_041 + X_050(dx) * lb->F_050 +
      X_104(dx) * lb->F_104 + X_113(dx) * lb->F_113 + X_122(dx) * lb->F_122 +
      X_131(dx) * lb->F_131 + X_140(dx) * lb->F_140 + X_203(dx) * lb->F_203 +
      X_212(dx) * lb->F_212 + X_221(dx) * lb->F_221 + X_230(dx) * lb->F_230 +
      X_302(dx) * lb->F_302 + X_311(dx) * lb->F_311 + X_320(dx) * lb->F_320 +
      X_401(dx) * lb->F_401 + X_410(dx) * lb->F_410 + X_500(dx) * lb->F_500;

  /* Shift 5th order field tensor terms (addition to rank 1) */
  la->F_001 +=
      X_004(dx) * lb->F_005 + X_013(dx) * lb->F_014 + X_022(dx) * lb->F_023 +
      X_031(dx) * lb->F_032 + X_040(dx) * lb->F_041 + X_103(dx) * lb->F_104 +
      X_112(dx) * lb->F_113 + X_121(dx) * lb->F_122 + X_130(dx) * lb->F_131 +
      X_202(dx) * lb->F_203 + X_211(dx) * lb->F_212 + X_220(dx) * lb->F_221 +
      X_301(dx) * lb->F_302 + X_310(dx) * lb->F_311 + X_400(dx) * lb->F_401;
  la->F_010 +=
      X_004(dx) * lb->F_014 + X_013(dx) * lb->F_023 + X_022(dx) * lb->F_032 +
      X_031(dx) * lb->F_041 + X_040(dx) * lb->F_050 + X_103(dx) * lb->F_113 +
      X_112(dx) * lb->F_122 + X_121(dx) * lb->F_131 + X_130(dx) * lb->F_140 +
      X_202(dx) * lb->F_212 + X_211(dx) * lb->F_221 + X_220(dx) * lb->F_230 +
      X_301(dx) * lb->F_311 + X_310(dx) * lb->F_320 + X_400(dx) * lb->F_410;
  la->F_100 +=
      X_004(dx) * lb->F_104 + X_013(dx) * lb->F_113 + X_022(dx) * lb->F_122 +
      X_031(dx) * lb->F_131 + X_040(dx) * lb->F_140 + X_103(dx) * lb->F_203 +
      X_112(dx) * lb->F_212 + X_121(dx) * lb->F_221 + X_130(dx) * lb->F_230 +
      X_202(dx) * lb->F_302 + X_211(dx) * lb->F_311 + X_220(dx) * lb->F_320 +
      X_301(dx) * lb->F_401 + X_310(dx) * lb->F_410 + X_400(dx) * lb->F_500;

  /* Shift 5th order field tensor terms (addition to rank 2) */
  la->F_002 += X_003(dx) * lb->F_005 + X_012(dx) * lb->F_014 +
               X_021(dx) * lb->F_023 + X_030(dx) * lb->F_032 +
               X_102(dx) * lb->F_104 + X_111(dx) * lb->F_113 +
               X_120(dx) * lb->F_122 + X_201(dx) * lb->F_203 +
               X_210(dx) * lb->F_212 + X_300(dx) * lb->F_302;
  la->F_011 += X_003(dx) * lb->F_014 + X_012(dx) * lb->F_023 +
               X_021(dx) * lb->F_032 + X_030(dx) * lb->F_041 +
               X_102(dx) * lb->F_113 + X_111(dx) * lb->F_122 +
               X_120(dx) * lb->F_131 + X_201(dx) * lb->F_212 +
               X_210(dx) * lb->F_221 + X_300(dx) * lb->F_311;
  la->F_020 += X_003(dx) * lb->F_023 + X_012(dx) * lb->F_032 +
               X_021(dx) * lb->F_041 + X_030(dx) * lb->F_050 +
               X_102(dx) * lb->F_122 + X_111(dx) * lb->F_131 +
               X_120(dx) * lb->F_140 + X_201(dx) * lb->F_221 +
               X_210(dx) * lb->F_230 + X_300(dx) * lb->F_320;
  la->F_101 += X_003(dx) * lb->F_104 + X_012(dx) * lb->F_113 +
               X_021(dx) * lb->F_122 + X_030(dx) * lb->F_131 +
               X_102(dx) * lb->F_203 + X_111(dx) * lb->F_212 +
               X_120(dx) * lb->F_221 + X_201(dx) * lb->F_302 +
               X_210(dx) * lb->F_311 + X_300(dx) * lb->F_401;
  la->F_110 += X_003(dx) * lb->F_113 + X_012(dx) * lb->F_122 +
               X_021(dx) * lb->F_131 + X_030(dx) * lb->F_140 +
               X_102(dx) * lb->F_212 + X_111(dx) * lb->F_221 +
               X_120(dx) * lb->F_230 + X_201(dx) * lb->F_311 +
               X_210(dx) * lb->F_320 + X_300(dx) * lb->F_410;
  la->F_200 += X_003(dx) * lb->F_203 + X_012(dx) * lb->F_212 +
               X_021(dx) * lb->F_221 + X_030(dx) * lb->F_230 +
               X_102(dx) * lb->F_302 + X_111(dx) * lb->F_311 +
               X_120(dx) * lb->F_320 + X_201(dx) * lb->F_401 +
               X_210(dx) * lb->F_410 + X_300(dx) * lb->F_500;

  /* Shift 5th order field tensor terms (addition to rank 3) */
  la->F_003 += X_002(dx) * lb->F_005 + X_011(dx) * lb->F_014 +
               X_020(dx) * lb->F_023 + X_101(dx) * lb->F_104 +
               X_110(dx) * lb->F_113 + X_200(dx) * lb->F_203;
  la->F_012 += X_002(dx) * lb->F_014 + X_011(dx) * lb->F_023 +
               X_020(dx) * lb->F_032 + X_101(dx) * lb->F_113 +
               X_110(dx) * lb->F_122 + X_200(dx) * lb->F_212;
  la->F_021 += X_002(dx) * lb->F_023 + X_011(dx) * lb->F_032 +
               X_020(dx) * lb->F_041 + X_101(dx) * lb->F_122 +
               X_110(dx) * lb->F_131 + X_200(dx) * lb->F_221;
  la->F_030 += X_002(dx) * lb->F_032 + X_011(dx) * lb->F_041 +
               X_020(dx) * lb->F_050 + X_101(dx) * lb->F_131 +
               X_110(dx) * lb->F_140 + X_200(dx) * lb->F_230;
  la->F_102 += X_002(dx) * lb->F_104 + X_011(dx) * lb->F_113 +
               X_020(dx) * lb->F_122 + X_101(dx) * lb->F_203 +
               X_110(dx) * lb->F_212 + X_200(dx) * lb->F_302;
  la->F_111 += X_002(dx) * lb->F_113 + X_011(dx) * lb->F_122 +
               X_020(dx) * lb->F_131 + X_101(dx) * lb->F_212 +
               X_110(dx) * lb->F_221 + X_200(dx) * lb->F_311;
  la->F_120 += X_002(dx) * lb->F_122 + X_011(dx) * lb->F_131 +
               X_020(dx) * lb->F_140 + X_101(dx) * lb->F_221 +
               X_110(dx) * lb->F_230 + X_200(dx) * lb->F_320;
  la->F_201 += X_002(dx) * lb->F_203 + X_011(dx) * lb->F_212 +
               X_020(dx) * lb->F_221 + X_101(dx) * lb->F_302 +
               X_110(dx) * lb->F_311 + X_200(dx) * lb->F_401;
  la->F_210 += X_002(dx) * lb->F_212 + X_011(dx) * lb->F_221 +
               X_020(dx) * lb->F_230 + X_101(dx) * lb->F_311 +
               X_110(dx) * lb->F_320 + X_200(dx) * lb->F_410;
  la->F_300 += X_002(dx) * lb->F_302 + X_011(dx) * lb->F_311 +
               X_020(dx) * lb->F_320 + X_101(dx) * lb->F_401 +
               X_110(dx) * lb->F_410 + X_200(dx) * lb->F_500;

  /* Shift 5th order field tensor terms (addition to rank 4) */
  la->F_004 +=
      X_001(dx) * lb->F_005 + X_010(dx) * lb->F_014 + X_100(dx) * lb->F_104;
  la->F_013 +=
      X_001(dx) * lb->F_014 + X_010(dx) * lb->F_023 + X_100(dx) * lb->F_113;
  la->F_022 +=
      X_001(dx) * lb->F_023 + X_010(dx) * lb->F_032 + X_100(dx) * lb->F_122;
  la->F_031 +=
      X_001(dx) * lb->F_032 + X_010(dx) * lb->F_041 + X_100(dx) * lb->F_131;
  la->F_040 +=
      X_001(dx) * lb->F_041 + X_010(dx) * lb->F_050 + X_100(dx) * lb->F_140;
  la->F_103 +=
      X_001(dx) * lb->F_104 + X_010(dx) * lb->F_113 + X_100(dx) * lb->F_203;
  la->F_112 +=
      X_001(dx) * lb->F_113 + X_010(dx) * lb->F_122 + X_100(dx) * lb->F_212;
  la->F_121 +=
      X_001(dx) * lb->F_122 + X_010(dx) * lb->F_131 + X_100(dx) * lb->F_221;
  la->F_130 +=
      X_001(dx) * lb->F_131 + X_010(dx) * lb->F_140 + X_100(dx) * lb->F_230;
  la->F_202 +=
      X_001(dx) * lb->F_203 + X_010(dx) * lb->F_212 + X_100(dx) * lb->F_302;
  la->F_211 +=
      X_001(dx) * lb->F_212 + X_010(dx) * lb->F_221 + X_100(dx) * lb->F_311;
  la->F_220 +=
      X_001(dx) * lb->F_221 + X_010(dx) * lb->F_230 + X_100(dx) * lb->F_320;
  la->F_301 +=
      X_001(dx) * lb->F_302 + X_010(dx) * lb->F_311 + X_100(dx) * lb->F_401;
  la->F_310 +=
      X_001(dx) * lb->F_311 + X_010(dx) * lb->F_320 + X_100(dx) * lb->F_410;
  la->F_400 +=
      X_001(dx) * lb->F_401 + X_010(dx) * lb->F_410 + X_100(dx) * lb->F_500;

  /* Shift 5th order field tensor terms (addition to rank 5) */
  la->F_005 += X_000(dx) * lb->F_005;
  la->F_014 += X_000(dx) * lb->F_014;
  la->F_023 += X_000(dx) * lb->F_023;
  la->F_032 += X_000(dx) * lb->F_032;
  la->F_041 += X_000(dx) * lb->F_041;
  la->F_050 += X_000(dx) * lb->F_050;
  la->F_104 += X_000(dx) * lb->F_104;
  la->F_113 += X_000(dx) * lb->F_113;
  la->F_122 += X_000(dx) * lb->F_122;
  la->F_131 += X_000(dx) * lb->F_131;
  la->F_140 += X_000(dx) * lb->F_140;
  la->F_203 += X_000(dx) * lb->F_203;
  la->F_212 += X_000(dx) * lb->F_212;
  la->F_221 += X_000(dx) * lb->F_221;
  la->F_230 += X_000(dx) * lb->F_230;
  la->F_302 += X_000(dx) * lb->F_302;
  la->F_311 += X_000(dx) * lb->F_311;
  la->F_320 += X_000(dx) * lb->F_320;
  la->F_401 += X_000(dx) * lb->F_401;
  la->F_410 += X_000(dx) * lb->F_410;
  la->F_500 += X_000(dx) * lb->F_500;

#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 5

  /* Shift 6th order field tensor terms (addition to rank 0) */
  la->F_000 +=  X_006(dx) * lb->F_006 X_015(dx) * lb->F_015
    + X_024(dx) * lb->F_024 X_033(dx) * lb->F_033
    + X_042(dx) * lb->F_042 X_051(dx) * lb->F_051
    + X_060(dx) * lb->F_060 X_105(dx) * lb->F_105
    + X_114(dx) * lb->F_114 X_123(dx) * lb->F_123
    + X_132(dx) * lb->F_132 X_141(dx) * lb->F_141
    + X_150(dx) * lb->F_150 X_204(dx) * lb->F_204
    + X_213(dx) * lb->F_213 X_222(dx) * lb->F_222
    + X_231(dx) * lb->F_231 X_240(dx) * lb->F_240
    + X_303(dx) * lb->F_303 X_312(dx) * lb->F_312
    + X_321(dx) * lb->F_321 X_330(dx) * lb->F_330
    + X_402(dx) * lb->F_402 X_411(dx) * lb->F_411
    + X_420(dx) * lb->F_420 X_501(dx) * lb->F_501
    + X_510(dx) * lb->F_510 X_600(dx) * lb->F_600;

  /* Shift 6th order field tensor terms (addition to rank 1) */
  la->F_001 +=  X_005(dx) * lb->F_006 X_014(dx) * lb->F_015
    + X_023(dx) * lb->F_024 X_032(dx) * lb->F_033
    + X_041(dx) * lb->F_042 X_050(dx) * lb->F_051
    + X_104(dx) * lb->F_105 X_113(dx) * lb->F_114
    + X_122(dx) * lb->F_123 X_131(dx) * lb->F_132
    + X_140(dx) * lb->F_141 X_203(dx) * lb->F_204
    + X_212(dx) * lb->F_213 X_221(dx) * lb->F_222
    + X_230(dx) * lb->F_231 X_302(dx) * lb->F_303
    + X_311(dx) * lb->F_312 X_320(dx) * lb->F_321
    + X_401(dx) * lb->F_402 X_410(dx) * lb->F_411
    + X_500(dx) * lb->F_501;
  la->F_010 +=  X_005(dx) * lb->F_015 X_014(dx) * lb->F_024
    + X_023(dx) * lb->F_033 X_032(dx) * lb->F_042
    + X_041(dx) * lb->F_051 X_050(dx) * lb->F_060
    + X_104(dx) * lb->F_114 X_113(dx) * lb->F_123
    + X_122(dx) * lb->F_132 X_131(dx) * lb->F_141
    + X_140(dx) * lb->F_150 X_203(dx) * lb->F_213
    + X_212(dx) * lb->F_222 X_221(dx) * lb->F_231
    + X_230(dx) * lb->F_240 X_302(dx) * lb->F_312
    + X_311(dx) * lb->F_321 X_320(dx) * lb->F_330
    + X_401(dx) * lb->F_411 X_410(dx) * lb->F_420
    + X_500(dx) * lb->F_510;
  la->F_100 +=  X_005(dx) * lb->F_105 X_014(dx) * lb->F_114
    + X_023(dx) * lb->F_123 X_032(dx) * lb->F_132
    + X_041(dx) * lb->F_141 X_050(dx) * lb->F_150
    + X_104(dx) * lb->F_204 X_113(dx) * lb->F_213
    + X_122(dx) * lb->F_222 X_131(dx) * lb->F_231
    + X_140(dx) * lb->F_240 X_203(dx) * lb->F_303
    + X_212(dx) * lb->F_312 X_221(dx) * lb->F_321
    + X_230(dx) * lb->F_330 X_302(dx) * lb->F_402
    + X_311(dx) * lb->F_411 X_320(dx) * lb->F_420
    + X_401(dx) * lb->F_501 X_410(dx) * lb->F_510
    + X_500(dx) * lb->F_600;

  /* Shift 6th order field tensor terms (addition to rank 2) */
  la->F_002 +=  X_004(dx) * lb->F_006 X_013(dx) * lb->F_015
    + X_022(dx) * lb->F_024 X_031(dx) * lb->F_033
    + X_040(dx) * lb->F_042 X_103(dx) * lb->F_105
    + X_112(dx) * lb->F_114 X_121(dx) * lb->F_123
    + X_130(dx) * lb->F_132 X_202(dx) * lb->F_204
    + X_211(dx) * lb->F_213 X_220(dx) * lb->F_222
    + X_301(dx) * lb->F_303 X_310(dx) * lb->F_312
    + X_400(dx) * lb->F_402;
  la->F_011 +=  X_004(dx) * lb->F_015 X_013(dx) * lb->F_024
    + X_022(dx) * lb->F_033 X_031(dx) * lb->F_042
    + X_040(dx) * lb->F_051 X_103(dx) * lb->F_114
    + X_112(dx) * lb->F_123 X_121(dx) * lb->F_132
    + X_130(dx) * lb->F_141 X_202(dx) * lb->F_213
    + X_211(dx) * lb->F_222 X_220(dx) * lb->F_231
    + X_301(dx) * lb->F_312 X_310(dx) * lb->F_321
    + X_400(dx) * lb->F_411;
  la->F_020 +=  X_004(dx) * lb->F_024 X_013(dx) * lb->F_033
    + X_022(dx) * lb->F_042 X_031(dx) * lb->F_051
    + X_040(dx) * lb->F_060 X_103(dx) * lb->F_123
    + X_112(dx) * lb->F_132 X_121(dx) * lb->F_141
    + X_130(dx) * lb->F_150 X_202(dx) * lb->F_222
    + X_211(dx) * lb->F_231 X_220(dx) * lb->F_240
    + X_301(dx) * lb->F_321 X_310(dx) * lb->F_330
    + X_400(dx) * lb->F_420;
  la->F_101 +=  X_004(dx) * lb->F_105 X_013(dx) * lb->F_114
    + X_022(dx) * lb->F_123 X_031(dx) * lb->F_132
    + X_040(dx) * lb->F_141 X_103(dx) * lb->F_204
    + X_112(dx) * lb->F_213 X_121(dx) * lb->F_222
    + X_130(dx) * lb->F_231 X_202(dx) * lb->F_303
    + X_211(dx) * lb->F_312 X_220(dx) * lb->F_321
    + X_301(dx) * lb->F_402 X_310(dx) * lb->F_411
    + X_400(dx) * lb->F_501;
  la->F_110 +=  X_004(dx) * lb->F_114 X_013(dx) * lb->F_123
    + X_022(dx) * lb->F_132 X_031(dx) * lb->F_141
    + X_040(dx) * lb->F_150 X_103(dx) * lb->F_213
    + X_112(dx) * lb->F_222 X_121(dx) * lb->F_231
    + X_130(dx) * lb->F_240 X_202(dx) * lb->F_312
    + X_211(dx) * lb->F_321 X_220(dx) * lb->F_330
    + X_301(dx) * lb->F_411 X_310(dx) * lb->F_420
    + X_400(dx) * lb->F_510;
  la->F_200 +=  X_004(dx) * lb->F_204 X_013(dx) * lb->F_213
    + X_022(dx) * lb->F_222 X_031(dx) * lb->F_231
    + X_040(dx) * lb->F_240 X_103(dx) * lb->F_303
    + X_112(dx) * lb->F_312 X_121(dx) * lb->F_321
    + X_130(dx) * lb->F_330 X_202(dx) * lb->F_402
    + X_211(dx) * lb->F_411 X_220(dx) * lb->F_420
    + X_301(dx) * lb->F_501 X_310(dx) * lb->F_510
    + X_400(dx) * lb->F_600;

  /* Shift 6th order field tensor terms (addition to rank 3) */
  la->F_003 +=  X_003(dx) * lb->F_006 X_012(dx) * lb->F_015
    + X_021(dx) * lb->F_024 X_030(dx) * lb->F_033
    + X_102(dx) * lb->F_105 X_111(dx) * lb->F_114
    + X_120(dx) * lb->F_123 X_201(dx) * lb->F_204
    + X_210(dx) * lb->F_213 X_300(dx) * lb->F_303;
  la->F_012 +=  X_003(dx) * lb->F_015 X_012(dx) * lb->F_024
    + X_021(dx) * lb->F_033 X_030(dx) * lb->F_042
    + X_102(dx) * lb->F_114 X_111(dx) * lb->F_123
    + X_120(dx) * lb->F_132 X_201(dx) * lb->F_213
    + X_210(dx) * lb->F_222 X_300(dx) * lb->F_312;
  la->F_021 +=  X_003(dx) * lb->F_024 X_012(dx) * lb->F_033
    + X_021(dx) * lb->F_042 X_030(dx) * lb->F_051
    + X_102(dx) * lb->F_123 X_111(dx) * lb->F_132
    + X_120(dx) * lb->F_141 X_201(dx) * lb->F_222
    + X_210(dx) * lb->F_231 X_300(dx) * lb->F_321;
  la->F_030 +=  X_003(dx) * lb->F_033 X_012(dx) * lb->F_042
    + X_021(dx) * lb->F_051 X_030(dx) * lb->F_060
    + X_102(dx) * lb->F_132 X_111(dx) * lb->F_141
    + X_120(dx) * lb->F_150 X_201(dx) * lb->F_231
    + X_210(dx) * lb->F_240 X_300(dx) * lb->F_330;
  la->F_102 +=  X_003(dx) * lb->F_105 X_012(dx) * lb->F_114
    + X_021(dx) * lb->F_123 X_030(dx) * lb->F_132
    + X_102(dx) * lb->F_204 X_111(dx) * lb->F_213
    + X_120(dx) * lb->F_222 X_201(dx) * lb->F_303
    + X_210(dx) * lb->F_312 X_300(dx) * lb->F_402;
  la->F_111 +=  X_003(dx) * lb->F_114 X_012(dx) * lb->F_123
    + X_021(dx) * lb->F_132 X_030(dx) * lb->F_141
    + X_102(dx) * lb->F_213 X_111(dx) * lb->F_222
    + X_120(dx) * lb->F_231 X_201(dx) * lb->F_312
    + X_210(dx) * lb->F_321 X_300(dx) * lb->F_411;
  la->F_120 +=  X_003(dx) * lb->F_123 X_012(dx) * lb->F_132
    + X_021(dx) * lb->F_141 X_030(dx) * lb->F_150
    + X_102(dx) * lb->F_222 X_111(dx) * lb->F_231
    + X_120(dx) * lb->F_240 X_201(dx) * lb->F_321
    + X_210(dx) * lb->F_330 X_300(dx) * lb->F_420;
  la->F_201 +=  X_003(dx) * lb->F_204 X_012(dx) * lb->F_213
    + X_021(dx) * lb->F_222 X_030(dx) * lb->F_231
    + X_102(dx) * lb->F_303 X_111(dx) * lb->F_312
    + X_120(dx) * lb->F_321 X_201(dx) * lb->F_402
    + X_210(dx) * lb->F_411 X_300(dx) * lb->F_501;
  la->F_210 +=  X_003(dx) * lb->F_213 X_012(dx) * lb->F_222
    + X_021(dx) * lb->F_231 X_030(dx) * lb->F_240
    + X_102(dx) * lb->F_312 X_111(dx) * lb->F_321
    + X_120(dx) * lb->F_330 X_201(dx) * lb->F_411
    + X_210(dx) * lb->F_420 X_300(dx) * lb->F_510;
  la->F_300 +=  X_003(dx) * lb->F_303 X_012(dx) * lb->F_312
    + X_021(dx) * lb->F_321 X_030(dx) * lb->F_330
    + X_102(dx) * lb->F_402 X_111(dx) * lb->F_411
    + X_120(dx) * lb->F_420 X_201(dx) * lb->F_501
    + X_210(dx) * lb->F_510 X_300(dx) * lb->F_600;

  /* Shift 6th order field tensor terms (addition to rank 4) */
  la->F_004 +=  X_002(dx) * lb->F_006 X_011(dx) * lb->F_015
    + X_020(dx) * lb->F_024 X_101(dx) * lb->F_105
    + X_110(dx) * lb->F_114 X_200(dx) * lb->F_204;
  la->F_013 +=  X_002(dx) * lb->F_015 X_011(dx) * lb->F_024
    + X_020(dx) * lb->F_033 X_101(dx) * lb->F_114
    + X_110(dx) * lb->F_123 X_200(dx) * lb->F_213;
  la->F_022 +=  X_002(dx) * lb->F_024 X_011(dx) * lb->F_033
    + X_020(dx) * lb->F_042 X_101(dx) * lb->F_123
    + X_110(dx) * lb->F_132 X_200(dx) * lb->F_222;
  la->F_031 +=  X_002(dx) * lb->F_033 X_011(dx) * lb->F_042
    + X_020(dx) * lb->F_051 X_101(dx) * lb->F_132
    + X_110(dx) * lb->F_141 X_200(dx) * lb->F_231;
  la->F_040 +=  X_002(dx) * lb->F_042 X_011(dx) * lb->F_051
    + X_020(dx) * lb->F_060 X_101(dx) * lb->F_141
    + X_110(dx) * lb->F_150 X_200(dx) * lb->F_240;
  la->F_103 +=  X_002(dx) * lb->F_105 X_011(dx) * lb->F_114
    + X_020(dx) * lb->F_123 X_101(dx) * lb->F_204
    + X_110(dx) * lb->F_213 X_200(dx) * lb->F_303;
  la->F_112 +=  X_002(dx) * lb->F_114 X_011(dx) * lb->F_123
    + X_020(dx) * lb->F_132 X_101(dx) * lb->F_213
    + X_110(dx) * lb->F_222 X_200(dx) * lb->F_312;
  la->F_121 +=  X_002(dx) * lb->F_123 X_011(dx) * lb->F_132
    + X_020(dx) * lb->F_141 X_101(dx) * lb->F_222
    + X_110(dx) * lb->F_231 X_200(dx) * lb->F_321;
  la->F_130 +=  X_002(dx) * lb->F_132 X_011(dx) * lb->F_141
    + X_020(dx) * lb->F_150 X_101(dx) * lb->F_231
    + X_110(dx) * lb->F_240 X_200(dx) * lb->F_330;
  la->F_202 +=  X_002(dx) * lb->F_204 X_011(dx) * lb->F_213
    + X_020(dx) * lb->F_222 X_101(dx) * lb->F_303
    + X_110(dx) * lb->F_312 X_200(dx) * lb->F_402;
  la->F_211 +=  X_002(dx) * lb->F_213 X_011(dx) * lb->F_222
    + X_020(dx) * lb->F_231 X_101(dx) * lb->F_312
    + X_110(dx) * lb->F_321 X_200(dx) * lb->F_411;
  la->F_220 +=  X_002(dx) * lb->F_222 X_011(dx) * lb->F_231
    + X_020(dx) * lb->F_240 X_101(dx) * lb->F_321
    + X_110(dx) * lb->F_330 X_200(dx) * lb->F_420;
  la->F_301 +=  X_002(dx) * lb->F_303 X_011(dx) * lb->F_312
    + X_020(dx) * lb->F_321 X_101(dx) * lb->F_402
    + X_110(dx) * lb->F_411 X_200(dx) * lb->F_501;
  la->F_310 +=  X_002(dx) * lb->F_312 X_011(dx) * lb->F_321
    + X_020(dx) * lb->F_330 X_101(dx) * lb->F_411
    + X_110(dx) * lb->F_420 X_200(dx) * lb->F_510;
  la->F_400 +=  X_002(dx) * lb->F_402 X_011(dx) * lb->F_411
    + X_020(dx) * lb->F_420 X_101(dx) * lb->F_501
    + X_110(dx) * lb->F_510 X_200(dx) * lb->F_600;

  /* Shift 6th order field tensor terms (addition to rank 5) */
  la->F_005 +=  X_001(dx) * lb->F_006 X_010(dx) * lb->F_015
    + X_100(dx) * lb->F_105;
  la->F_014 +=  X_001(dx) * lb->F_015 X_010(dx) * lb->F_024
    + X_100(dx) * lb->F_114;
  la->F_023 +=  X_001(dx) * lb->F_024 X_010(dx) * lb->F_033
    + X_100(dx) * lb->F_123;
  la->F_032 +=  X_001(dx) * lb->F_033 X_010(dx) * lb->F_042
    + X_100(dx) * lb->F_132;
  la->F_041 +=  X_001(dx) * lb->F_042 X_010(dx) * lb->F_051
    + X_100(dx) * lb->F_141;
  la->F_050 +=  X_001(dx) * lb->F_051 X_010(dx) * lb->F_060
    + X_100(dx) * lb->F_150;
  la->F_104 +=  X_001(dx) * lb->F_105 X_010(dx) * lb->F_114
    + X_100(dx) * lb->F_204;
  la->F_113 +=  X_001(dx) * lb->F_114 X_010(dx) * lb->F_123
    + X_100(dx) * lb->F_213;
  la->F_122 +=  X_001(dx) * lb->F_123 X_010(dx) * lb->F_132
    + X_100(dx) * lb->F_222;
  la->F_131 +=  X_001(dx) * lb->F_132 X_010(dx) * lb->F_141
    + X_100(dx) * lb->F_231;
  la->F_140 +=  X_001(dx) * lb->F_141 X_010(dx) * lb->F_150
    + X_100(dx) * lb->F_240;
  la->F_203 +=  X_001(dx) * lb->F_204 X_010(dx) * lb->F_213
    + X_100(dx) * lb->F_303;
  la->F_212 +=  X_001(dx) * lb->F_213 X_010(dx) * lb->F_222
    + X_100(dx) * lb->F_312;
  la->F_221 +=  X_001(dx) * lb->F_222 X_010(dx) * lb->F_231
    + X_100(dx) * lb->F_321;
  la->F_230 +=  X_001(dx) * lb->F_231 X_010(dx) * lb->F_240
    + X_100(dx) * lb->F_330;
  la->F_302 +=  X_001(dx) * lb->F_303 X_010(dx) * lb->F_312
    + X_100(dx) * lb->F_402;
  la->F_311 +=  X_001(dx) * lb->F_312 X_010(dx) * lb->F_321
    + X_100(dx) * lb->F_411;
  la->F_320 +=  X_001(dx) * lb->F_321 X_010(dx) * lb->F_330
    + X_100(dx) * lb->F_420;
  la->F_401 +=  X_001(dx) * lb->F_402 X_010(dx) * lb->F_411
    + X_100(dx) * lb->F_501;
  la->F_410 +=  X_001(dx) * lb->F_411 X_010(dx) * lb->F_420
    + X_100(dx) * lb->F_510;
  la->F_500 +=  X_001(dx) * lb->F_501 X_010(dx) * lb->F_510
    + X_100(dx) * lb->F_600;

  /* Shift 6th order field tensor terms (addition to rank 6) */
  la->F_006 +=  X_000(dx) * lb->F_006;
  la->F_015 +=  X_000(dx) * lb->F_015;
  la->F_024 +=  X_000(dx) * lb->F_024;
  la->F_033 +=  X_000(dx) * lb->F_033;
  la->F_042 +=  X_000(dx) * lb->F_042;
  la->F_051 +=  X_000(dx) * lb->F_051;
  la->F_060 +=  X_000(dx) * lb->F_060;
  la->F_105 +=  X_000(dx) * lb->F_105;
  la->F_114 +=  X_000(dx) * lb->F_114;
  la->F_123 +=  X_000(dx) * lb->F_123;
  la->F_132 +=  X_000(dx) * lb->F_132;
  la->F_141 +=  X_000(dx) * lb->F_141;
  la->F_150 +=  X_000(dx) * lb->F_150;
  la->F_204 +=  X_000(dx) * lb->F_204;
  la->F_213 +=  X_000(dx) * lb->F_213;
  la->F_222 +=  X_000(dx) * lb->F_222;
  la->F_231 +=  X_000(dx) * lb->F_231;
  la->F_240 +=  X_000(dx) * lb->F_240;
  la->F_303 +=  X_000(dx) * lb->F_303;
  la->F_312 +=  X_000(dx) * lb->F_312;
  la->F_321 +=  X_000(dx) * lb->F_321;
  la->F_330 +=  X_000(dx) * lb->F_330;
  la->F_402 +=  X_000(dx) * lb->F_402;
  la->F_411 +=  X_000(dx) * lb->F_411;
  la->F_420 +=  X_000(dx) * lb->F_420;
  la->F_501 +=  X_000(dx) * lb->F_501;
  la->F_510 +=  X_000(dx) * lb->F_510;
  la->F_600 +=  X_000(dx) * lb->F_600;

#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 6
#error "Missing implementation for order >6"
#endif
}

/**
 * @brief Applies the  #grav_tensor to a  #gpart.
 *
 * Corresponds to equation (28a).
 *
 * @param lb The gravity field tensor to apply.
 * @param loc The position of the gravity field tensor.
 * @param gp The #gpart to update.
 */
INLINE static void gravity_L2P(const struct grav_tensor *lb,
                               const double loc[3], struct gpart *gp) {

#ifdef SWIFT_DEBUG_CHECKS
  if (lb->num_interacted + lb->num_not_interacted == 0)
      error("Interacting with empty field tensor");
#endif
#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_GRAVITY_FORCE_CHECKS)
  gp->num_interacted += lb->num_interacted + lb->num_not_interacted;
#endif
#ifdef SWIFT_GRAVITY_FORCE_CHECKS
  gp->num_interacted_m2m += lb->num_interacted;
  gp->num_not_interacted += lb->num_not_interacted;
#endif

  /* Local accumulator */
  double a_grav[3] = {0., 0., 0.};
  double pot = 0.;

  /* Distance to the multipole */
  const double dx[3] = {gp->x[0] - loc[0], gp->x[1] - loc[1],
                        gp->x[2] - loc[2]};
  /* 0th order contributions */
  pot -= X_000(dx) * lb->F_000;

#if SELF_GRAVITY_MULTIPOLE_ORDER > 0
  /* 1st order contributions */
  a_grav[0] += X_000(dx) * lb->F_100;
  a_grav[1] += X_000(dx) * lb->F_010;
  a_grav[2] += X_000(dx) * lb->F_001;

  pot -= X_001(dx) * lb->F_001 + X_010(dx) * lb->F_010 + X_100(dx) * lb->F_100;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 1

  /* 2nd order contributions */
  a_grav[0] +=
      X_100(dx) * lb->F_200 + X_010(dx) * lb->F_110 + X_001(dx) * lb->F_101;
  a_grav[1] +=
      X_100(dx) * lb->F_110 + X_010(dx) * lb->F_020 + X_001(dx) * lb->F_011;
  a_grav[2] +=
      X_100(dx) * lb->F_101 + X_010(dx) * lb->F_011 + X_001(dx) * lb->F_002;

  pot -= X_002(dx) * lb->F_002 + X_011(dx) * lb->F_011 + X_020(dx) * lb->F_020 +
         X_101(dx) * lb->F_101 + X_110(dx) * lb->F_110 + X_200(dx) * lb->F_200;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 2

  /* 3rd order contributions */
  a_grav[0] +=
      X_200(dx) * lb->F_300 + X_020(dx) * lb->F_120 + X_002(dx) * lb->F_102;
  a_grav[0] +=
      X_110(dx) * lb->F_210 + X_101(dx) * lb->F_201 + X_011(dx) * lb->F_111;
  a_grav[1] +=
      X_200(dx) * lb->F_210 + X_020(dx) * lb->F_030 + X_002(dx) * lb->F_012;
  a_grav[1] +=
      X_110(dx) * lb->F_120 + X_101(dx) * lb->F_111 + X_011(dx) * lb->F_021;
  a_grav[2] +=
      X_200(dx) * lb->F_201 + X_020(dx) * lb->F_021 + X_002(dx) * lb->F_003;
  a_grav[2] +=
      X_110(dx) * lb->F_111 + X_101(dx) * lb->F_102 + X_011(dx) * lb->F_012;

  pot -= X_003(dx) * lb->F_003 + X_012(dx) * lb->F_012 + X_021(dx) * lb->F_021 +
         X_030(dx) * lb->F_030 + X_102(dx) * lb->F_102 + X_111(dx) * lb->F_111 +
         X_120(dx) * lb->F_120 + X_201(dx) * lb->F_201 + X_210(dx) * lb->F_210 +
         X_300(dx) * lb->F_300;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 3

  /* 4th order contributions */
  a_grav[0] += X_003(dx) * lb->F_103 + X_012(dx) * lb->F_112 +
               X_021(dx) * lb->F_121 + X_030(dx) * lb->F_130 +
               X_102(dx) * lb->F_202 + X_111(dx) * lb->F_211 +
               X_120(dx) * lb->F_220 + X_201(dx) * lb->F_301 +
               X_210(dx) * lb->F_310 + X_300(dx) * lb->F_400;
  a_grav[1] += X_003(dx) * lb->F_013 + X_012(dx) * lb->F_022 +
               X_021(dx) * lb->F_031 + X_030(dx) * lb->F_040 +
               X_102(dx) * lb->F_112 + X_111(dx) * lb->F_121 +
               X_120(dx) * lb->F_130 + X_201(dx) * lb->F_211 +
               X_210(dx) * lb->F_220 + X_300(dx) * lb->F_310;
  a_grav[2] += X_003(dx) * lb->F_004 + X_012(dx) * lb->F_013 +
               X_021(dx) * lb->F_022 + X_030(dx) * lb->F_031 +
               X_102(dx) * lb->F_103 + X_111(dx) * lb->F_112 +
               X_120(dx) * lb->F_121 + X_201(dx) * lb->F_202 +
               X_210(dx) * lb->F_211 + X_300(dx) * lb->F_301;

  pot -= X_004(dx) * lb->F_004 + X_013(dx) * lb->F_013 + X_022(dx) * lb->F_022 +
         X_031(dx) * lb->F_031 + X_040(dx) * lb->F_040 + X_103(dx) * lb->F_103 +
         X_112(dx) * lb->F_112 + X_121(dx) * lb->F_121 + X_130(dx) * lb->F_130 +
         X_202(dx) * lb->F_202 + X_211(dx) * lb->F_211 + X_220(dx) * lb->F_220 +
         X_301(dx) * lb->F_301 + X_310(dx) * lb->F_310 + X_400(dx) * lb->F_400;

#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 4

  /* 5th order contributions */
  a_grav[0] +=
      X_004(dx) * lb->F_104 + X_013(dx) * lb->F_113 + X_022(dx) * lb->F_122 +
      X_031(dx) * lb->F_131 + X_040(dx) * lb->F_140 + X_103(dx) * lb->F_203 +
      X_112(dx) * lb->F_212 + X_121(dx) * lb->F_221 + X_130(dx) * lb->F_230 +
      X_202(dx) * lb->F_302 + X_211(dx) * lb->F_311 + X_220(dx) * lb->F_320 +
      X_301(dx) * lb->F_401 + X_310(dx) * lb->F_410 + X_400(dx) * lb->F_500;
  a_grav[1] +=
      X_004(dx) * lb->F_014 + X_013(dx) * lb->F_023 + X_022(dx) * lb->F_032 +
      X_031(dx) * lb->F_041 + X_040(dx) * lb->F_050 + X_103(dx) * lb->F_113 +
      X_112(dx) * lb->F_122 + X_121(dx) * lb->F_131 + X_130(dx) * lb->F_140 +
      X_202(dx) * lb->F_212 + X_211(dx) * lb->F_221 + X_220(dx) * lb->F_230 +
      X_301(dx) * lb->F_311 + X_310(dx) * lb->F_320 + X_400(dx) * lb->F_410;
  a_grav[2] +=
      X_004(dx) * lb->F_005 + X_013(dx) * lb->F_014 + X_022(dx) * lb->F_023 +
      X_031(dx) * lb->F_032 + X_040(dx) * lb->F_041 + X_103(dx) * lb->F_104 +
      X_112(dx) * lb->F_113 + X_121(dx) * lb->F_122 + X_130(dx) * lb->F_131 +
      X_202(dx) * lb->F_203 + X_211(dx) * lb->F_212 + X_220(dx) * lb->F_221 +
      X_301(dx) * lb->F_302 + X_310(dx) * lb->F_311 + X_400(dx) * lb->F_401;

  pot -= X_005(dx) * lb->F_005 + X_014(dx) * lb->F_014 + X_023(dx) * lb->F_023 +
         X_032(dx) * lb->F_032 + X_041(dx) * lb->F_041 + X_050(dx) * lb->F_050 +
         X_104(dx) * lb->F_104 + X_113(dx) * lb->F_113 + X_122(dx) * lb->F_122 +
         X_131(dx) * lb->F_131 + X_140(dx) * lb->F_140 + X_203(dx) * lb->F_203 +
         X_212(dx) * lb->F_212 + X_221(dx) * lb->F_221 + X_230(dx) * lb->F_230 +
         X_302(dx) * lb->F_302 + X_311(dx) * lb->F_311 + X_320(dx) * lb->F_320 +
         X_401(dx) * lb->F_401 + X_410(dx) * lb->F_410 + X_500(dx) * lb->F_500;

#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 5

  /* 6th order contributions */
  a_grav[0] += X_005(dx) * lb->F_105 + X_014(dx) * lb->F_114
    + X_023(dx) * lb->F_123 + X_032(dx) * lb->F_132
    + X_041(dx) * lb->F_141 + X_050(dx) * lb->F_150
    + X_104(dx) * lb->F_204 + X_113(dx) * lb->F_213
    + X_122(dx) * lb->F_222 + X_131(dx) * lb->F_231
    + X_140(dx) * lb->F_240 + X_203(dx) * lb->F_303
    + X_212(dx) * lb->F_312 + X_221(dx) * lb->F_321
    + X_230(dx) * lb->F_330 + X_302(dx) * lb->F_402
    + X_311(dx) * lb->F_411 + X_320(dx) * lb->F_420
    + X_401(dx) * lb->F_501 + X_410(dx) * lb->F_510
    + X_500(dx) * lb->F_600;
  a_grav[1] += X_005(dx) * lb->F_015 + X_014(dx) * lb->F_024
    + X_023(dx) * lb->F_033 + X_032(dx) * lb->F_042
    + X_041(dx) * lb->F_051 + X_050(dx) * lb->F_060
    + X_104(dx) * lb->F_114 + X_113(dx) * lb->F_123
    + X_122(dx) * lb->F_132 + X_131(dx) * lb->F_141
    + X_140(dx) * lb->F_150 + X_203(dx) * lb->F_213
    + X_212(dx) * lb->F_222 + X_221(dx) * lb->F_231
    + X_230(dx) * lb->F_240 + X_302(dx) * lb->F_312
    + X_311(dx) * lb->F_321 + X_320(dx) * lb->F_330
    + X_401(dx) * lb->F_411 + X_410(dx) * lb->F_420
    + X_500(dx) * lb->F_510;
  a_grav[2] += X_005(dx) * lb->F_006 + X_014(dx) * lb->F_015
    + X_023(dx) * lb->F_024 + X_032(dx) * lb->F_033
    + X_041(dx) * lb->F_042 + X_050(dx) * lb->F_051
    + X_104(dx) * lb->F_105 + X_113(dx) * lb->F_114
    + X_122(dx) * lb->F_123 + X_131(dx) * lb->F_132
    + X_140(dx) * lb->F_141 + X_203(dx) * lb->F_204
    + X_212(dx) * lb->F_213 + X_221(dx) * lb->F_222
    + X_230(dx) * lb->F_231 + X_302(dx) * lb->F_303
    + X_311(dx) * lb->F_312 + X_320(dx) * lb->F_321
    + `X_401(dx) * lb->F_402 + X_410(dx) * lb->F_411
    + X_500(dx) * lb->F_501;

  pot -= X_006(dx) * lb->F_006 + X_015(dx) * lb->F_015
    + X_024(dx) * lb->F_024 + X_033(dx) * lb->F_033
    + X_042(dx) * lb->F_042 + X_051(dx) * lb->F_051
    + X_060(dx) * lb->F_060 + X_105(dx) * lb->F_105
    + X_114(dx) * lb->F_114 + X_123(dx) * lb->F_123
    + X_132(dx) * lb->F_132 + X_141(dx) * lb->F_141
    + X_150(dx) * lb->F_150 + X_204(dx) * lb->F_204
    + X_213(dx) * lb->F_213 + X_222(dx) * lb->F_222
    + X_231(dx) * lb->F_231 + X_240(dx) * lb->F_240
    + X_303(dx) * lb->F_303 + X_312(dx) * lb->F_312
    + X_321(dx) * lb->F_321 + X_330(dx) * lb->F_330
    + X_402(dx) * lb->F_402 + X_411(dx) * lb->F_411
    + X_420(dx) * lb->F_420 + X_501(dx) * lb->F_501
    + X_510(dx) * lb->F_510 + X_600(dx) * lb->F_600;

#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 6
#error "Missing implementation for order >6"
#endif

  /* Update the particle */
  gp->a_grav[0] += a_grav[0];
  gp->a_grav[1] += a_grav[1];
  gp->a_grav[2] += a_grav[2];
  gravity_add_comoving_potential(gp, pot);

#ifdef SWIFT_GRAVITY_FORCE_CHECKS
  gp->a_grav_m2m[0] += a_grav[0];
  gp->a_grav_m2m[1] += a_grav[1];
  gp->a_grav_m2m[2] += a_grav[2];
#endif
}

/**
 * @brief Checks whether a cell-cell interaction can be appromixated by a M-M
 * interaction using the distance and cell radius.
 *
 * We use the multipole acceptance criterion of Dehnen, 2002, JCoPh, Volume 179,
 * Issue 1, pp.27-42, equation 10.
 *
 * @param r_crit_a The size of the multipole A.
 * @param r_crit_b The size of the multipole B.
 * @param theta_crit2 The square of the critical opening angle.
 * @param r2 Square of the distance (periodically wrapped) between the
 * multipoles.
 */
__attribute__((always_inline, const)) INLINE static int gravity_M2L_accept(
    const double r_crit_a, const double r_crit_b, const double theta_crit2,
    const double r2) {

  const double size = r_crit_a + r_crit_b;
  const double size2 = size * size;

  // MATTHIEU: Make this mass-dependent ?

  /* Multipole acceptance criterion (Dehnen 2002, eq.10) */
  return (r2 * theta_crit2 > size2);
}

/**
 * @brief Checks whether a particle-cell interaction can be appromixated by a
 * M2P interaction using the distance and cell radius.
 *
 * We use the multipole acceptance criterion of Dehnen, 2002, JCoPh, Volume 179,
 * Issue 1, pp.27-42, equation 10.
 *
 * @param r_max2 The square of the size of the multipole.
 * @param theta_crit2 The square of the critical opening angle.
 * @param r2 Square of the distance (periodically wrapped) between the
 * particle and the multipole.
 */
__attribute__((always_inline, const)) INLINE static int gravity_M2P_accept(
    const float r_max2, const float theta_crit2, const float r2) {

  // MATTHIEU: Make this mass-dependent ?

  /* Multipole acceptance criterion (Dehnen 2002, eq.10) */
  return (r2 * theta_crit2 > r_max2);
}

/**
 * @brief Checks whether a particle-cell interaction can be appromixated by a
 * M2P interaction using an advanced opening criteria.
 *
 * We use the multipole acceptance criterion of Dehnen, 2014, Computational 
 * Astrophysics and Cosmology, Volume 1, article id.1 24pp equations 12->16.
 *
 * Always use the classical criteria on the 0th timestep (as we have no 
 * accelerations yet). 
 *
 * If ADVANCED_OPENING_CRITERIA is not defined, always use the classical case.
 *
 * @param gp_a The reference gpart.
 * @param m_b The multipole of the cell we are interacting with.
 * @param r_max2 The square of the size of the multipole.
 * @param r2 Square of the distance (periodically wrapped) between the
 * particle and the multipole.
 * @param gravity_props The gravity properties struct, for some constants.
 * @param step The current timestep.
 */
__attribute__((always_inline, const)) INLINE static int gravity_M2P_accept_advanced(
    const struct gpart *gp_a, const struct multipole *m_b, const double r_max2,
    const double r2, const struct gravity_props *gravity_props, const int step) {

  if (step == 0) {
#ifdef ADVANCED_OPENING_CRITERIA
    /* Never accept M2P on the first timestep */
    return 0;
#else
    /* Classic criteria */
    const double theta_crit2 = gravity_props->theta_crit2;
    return gravity_M2P_accept(r_max2, theta_crit2, r2);
#endif

  } else {
#ifdef ADVANCED_OPENING_CRITERIA
    /* gpart and multipole information */
    const double min_a_grav_norm = m_b->min_a_grav_norm;
    const double M_a = gp_a->mass;

    /* Advanced acceptance criteria */
    const double r_p = r2 * pow(r2, SELF_GRAVITY_MULTIPOLE_ORDER-2);
    const double r_max_p = pow(r_max2, SELF_GRAVITY_MULTIPOLE_ORDER-2);

    const double E = r_max_p / r_p;
    const double E_bar = 8 * E;

    const double theta2 = r_max2 / r2;
    const double rel_force_error = gravity_props->rel_force_error;

    return ((E_bar * M_a / r2) < (rel_force_error * min_a_grav_norm) && (theta2 < 1.0));
#else
    /* Classic criteria */
    const double theta_crit2 = gravity_props->theta_crit2;
    return gravity_M2P_accept(r_max2, theta_crit2, r2);
#endif
  }
}

/**
 * @brief Checks whether a cell-cell interaction can be appromixated by a M-M
 * interaction using an advanced opening criteria.
 *
 * We use the multipole acceptance criterion of Dehnen, 2014, Computational 
 * Astrophysics and Cosmology, Volume 1, article id.1 24pp equations 12->16.
 *
 * Rather than the pure geometrical version of Dehnen 2002, this attempts to 
 * estimate the opening criteria based upon the relative force error, using the 
 * accelerations from the previous timestep.
 *
 * We therefore cannot use this criteria on the 0th timestep, and revert to the
 * 'classical' criteria to obtain the initial accelerations.
 *
 * If ADVANCED_OPENING_CRITERIA is not defined, always use the classical case.
 *
 * @param m_a The m_pole of cell A.
 * @param m_b The m_pole of cell B.
 * @param r_crit_a The size of the multipole A.
 * @param r_crit_b The size of the multipole B.
 * @param r2 Square of the distance (periodically wrapped) between the
 * multipoles.
 * @param gravity_props The gravity_properties struct, for some constants.
 * @param step The current timestep.
 */
__attribute__((always_inline, const)) INLINE static int gravity_M2L_accept_advanced(
    const struct multipole *m_a, const struct multipole *m_b,
    const double r_crit_a, const double r_crit_b, const double r2,
    const struct gravity_props *gravity_props, const int step) {

  /* Get some constants from the engine */
  const double theta_crit2 = gravity_props->theta_crit2;

  if (step == 0) {
    /* Need to use classical critetia on the first timestep
     * as we have no accelerations yet */
    return gravity_M2L_accept(r_crit_a, r_crit_b, theta_crit2, r2);

  } else {
#ifdef ADVANCED_OPENING_CRITERIA
    /* Some constatnts for the calculation */
    const double M_a = m_a->M_000;
    const double min_a_grav_norm = m_b->min_a_grav_norm;
    const double size = r_crit_a + r_crit_b;
    const double size2 = size * size;
    const double r_p = r2 * pow(r2, SELF_GRAVITY_MULTIPOLE_ORDER-2);
    const double theta2 = size2 / r2;
    const double rel_force_error = gravity_props->rel_force_error;

    /* Binomial coefficients for m=0-5 */
    const int binom_coefficient[6][6] = {{1,0,0,0,0,0},
                                         {1,1,0,0,0,0},
                                         {1,2,1,0,0,0},
                                         {1,3,3,1,0,0},
                                         {1,4,6,4,1,0},
                                         {1,5,10,10,5,1}};
  
    /* Sum over each multipole order to evaluate E_a->b (eq. 13) */
    double E = 0.0;
    for (int i=0; i <= SELF_GRAVITY_MULTIPOLE_ORDER; i++) {
      E += binom_coefficient[SELF_GRAVITY_MULTIPOLE_ORDER][i] *
              sqrt(m_a->power[i]) *
              pow(r_crit_b, SELF_GRAVITY_MULTIPOLE_ORDER-i);   
    }
    E /= M_a;
    E /= r_p;
  
    /* Compute E_a->b_bar (eq. 15)
     * This predicts the maximum(ish) relative force error (da/a) that you would expect
     * between the two multipoles via an M2L interaction i.e., da/a <= E_bar */
    const double E_bar = E * ((8 * max(r_crit_a, r_crit_b)) / (r_crit_a + r_crit_b));

    return ((E_bar * M_a / r2) < (rel_force_error * min_a_grav_norm) && (theta2 < 1.0));
#else
    return gravity_M2L_accept(r_crit_a, r_crit_b, theta_crit2, r2);
#endif
  }
}
#endif /* SWIFT_MULTIPOLE_H */
