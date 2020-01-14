#ifndef CONFIG_H
#define CONFIG_H

// User defined parameters
#define OUTPUT_DIR "output/"
#define GRID_WIDTH 64   // cells
#define BOX_WIDTH 378.944875363  // Mpc

// Spectral index and power spectrum normalization
#define N_S 0.9619
#define A_S 2.215e-9
#define PIVOT_SCALE 0.05 // in 1/Mpc

// Present-day fluctuation scale
#define SIGMA_8 0.8282

// Whether to normalize using sigma_8 (NORM_SIGMA) or A_S (NORM_CMB)
#define NORM_SIGMA 0
#define NORM_CMB 1
#define NORMALIZATION_METHOD    NORM_CMB

// Initial velocity methods
#define VEL_ZELDOVICH 0
#define VEL_TRANSFER 1
#define VELOCITY_METHOD     VEL_ZELDOVICH

// Whether to set the initial displacements with Zel'dovich or Infinite-LPT
#define FALSE 0
#define TRUE 1
#define USE_INFINI_LPT TRUE

// Whether to infer the velocity of the particles from the original grid/glass
// position (FALSE) or from the velocity field at the displaced position (TRUE)
#define VELOCITY_AT_DISPLACED_POS TRUE

// Hubble constant
#define H_0 67.556  // km/s/Mpc

// If particles are generated from a grid, the grid should be size NP^3
#define NP 16
#define PARTICLE_NUM NP* NP* NP

// If neutrino particles are generated from a grid, the grid should be size
// NNUP^3
#define NNUP 16
#define NEUTRINO_NUM NNUP* NNUP* NNUP

// Starting redshift of the simulation
#define Z_START 40;

// Comsological parameters are in cosmo.h

#endif
