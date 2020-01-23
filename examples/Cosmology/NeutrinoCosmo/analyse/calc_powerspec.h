/**
 * @calc_powerspec.h
 * @author  Willem Elbers <whe@willemelbers.com>
 * @version 1.0
 *
 * @section DESCRIPTION
 *
 * Calculates the 3D power spectrum from a grid using DFT.
 */

#include <iostream>
#include <vector>

#ifndef FFTW3_H_FUN
	#include <fftw3.h>
#endif

#ifndef FIELD_IO_H
	#define FIELD_IO_H
	#include "field_io.h"
#endif

#define CLOUD_IN_CELL 0

#include "calc_powerspec.cpp"
