/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *               2015 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_VECTOR_H
#define SWIFT_VECTOR_H
#include <stdio.h>

/* Config parameters. */
#include "../config.h"

/* Local headers */
#include "inline.h"

#ifdef WITH_VECTORIZATION

/* Need to check whether compiler supports this (IBM does not)
   This will prevent the macros to be defined and switch off
   explicit vectorization if the compiled does not support it */
#ifdef HAVE_IMMINTRIN_H
/* Include the header file with the intrinsics for Intel architecture. */
#include <immintrin.h>
#endif

/* Define the vector macro. */
#define VEC_MACRO(elcount, type) \
  __attribute__((vector_size((elcount) * sizeof(type)))) type

/* So what will the vector size be? */

#include <arm_sve.h>
#define VEC_SIZE 4
static const int zero_int = 0;
static const float zero_float = 0.f;
static const int one_int = 1;
static const float one_float = 1.f;
#define PRED_TRUE_32 svptrue_b32()
#define VEC_FLOAT svfloat32_t
#define VEC_DBL svfloat64_t
#define VEC_INT svint32_t
#define VEC_UINT svuint32_t
#define mask_t svbool_t
#define vec_load(a) svld1(PRED_TRUE_32,a)
#define vec_store(a, addr) svst1(PRED_TRUE_32,addr,a) 
#define vec_setzero() svld1(PRED_TRUE_32, &zero_float)
#define vec_setintzero() svld1(PRED_TRUE_32, &zero_int)
#define vec_set1() svld1(PRED_TRUE_32, &one_float)
#define vec_setint1() svld1(PRED_TRUE_32, &one_int)
#define vec_add(a, b) svadd_x(PRED_TRUE_32, a, b)
#define vec_mask_add(a, b, mask) svadd_x(mask, a, b)
#define vec_sub(a, b) svsub_x(PRED_TRUE_32, a, b)
#define vec_mask_sub(a, b, mask) svsub_x(mask, a, b)
#define vec_fma(a, b, c) svmad_x(PRED_TRUE_32, a, b, c)
#define vec_fnma(a, b, c) svmsb_x(PRED_TRUE_32, a, b, c)

/* Define the composite types for element access. */
typedef union {
  float f[VEC_SIZE];
  double d[VEC_SIZE / 2];
  int i[VEC_SIZE];
  int ui[VEC_SIZE];
} vector_array;

//inline void print_vector(vector v) {
//  printf("vector %.5e %.5e %.5e %.5e\n", v.f[0], v.f[1], v.f[2], v.f[3]);
//  fflush(stdout);
//}
//
///* Define the mask type depending on the instruction set used. */
//#ifdef HAVE_AVX512_F
//typedef __mmask16 mask_t;
//#else
//typedef vector mask_t;
//#endif
//
///**
// * @brief Calculates the inverse ($1/x$) of a vector using intrinsics and a
// * Newton iteration to obtain the correct level of accuracy.
// *
// * @param x #vector to be inverted.
// * @return x_inv #vector inverted x.
// */
//__attribute__((always_inline)) INLINE vector vec_reciprocal(vector x) {
//
//  vector x_inv;
//
//  x_inv.v = vec_rcp(x.v);
//  x_inv.v = vec_sub(x_inv.v,
//                    vec_mul(x_inv.v, (vec_fma(x.v, x_inv.v, vec_set1(-1.0f)))));
//
//  return x_inv;
//}
//
///**
// * @brief Calculates the inverse and square root (\f$1/\sqrt{x}\f$) of a vector
// * using intrinsics and a Newton iteration to obtain the correct level of
// * accuracy.
// *
// * @param x #vector to be inverted.
// * @return x_inv #vector inverted x.
// */
//__attribute__((always_inline)) INLINE vector vec_reciprocal_sqrt(vector x) {
//
//  vector x_inv;
//
//  x_inv.v = vec_rsqrt(x.v);
//  x_inv.v = vec_sub(
//      x_inv.v,
//      vec_mul(vec_mul(vec_set1(0.5f), x_inv.v),
//              (vec_fma(x.v, vec_mul(x_inv.v, x_inv.v), vec_set1(-1.0f)))));
//
//  return x_inv;
//}
//
///**
// * @brief Returns a new vector with data loaded from a memory address.
// *
// * @param x memory address to load from.
// * @return Loaded #vector.
// */
//__attribute__((always_inline)) INLINE vector vector_load(float *const x) {
//
//  vector temp;
//  temp.v = vec_load(x);
//  return temp;
//}
//
///**
// * @brief Returns a vector filled with one value.
// *
// * @param x value to set each element.
// * @return A #vector filled with a given constant.
// */
//__attribute__((always_inline)) INLINE vector vector_set1(const float x) {
//
//  vector temp;
//  temp.v = vec_set1(x);
//  return temp;
//}
//
///**
// * @brief Returns a new vector filled with zeros.
// *
// * @return temp set #vector.
// * @return A #vector filled with zeros.
// */
//__attribute__((always_inline)) INLINE vector vector_setzero(void) {
//
//  vector temp;
//  temp.v = vec_setzero();
//  return temp;
//}

#else
/* Needed for cache alignment. */
#define VEC_SIZE 8
#endif /* WITH_VECTORIZATION */

#endif /* SWIFT_VECTOR_H */
