#ifndef SWIFT_VECTOR_TESTS
#define SWIFT_VECTOR_TESTS
// #define __ARM_NEON 1
#define WITH_VECTORIZATION 1
#define __ARM_FEATURE_SVE 1
#include "../../config.h"
// #include "../../src/swift.h"
#undef __ARM_NEON
#include "../../src/minmax.h"
#include "../../src/vector_sve.h"
#include <stdio.h>
#include <math.h>

#endif
