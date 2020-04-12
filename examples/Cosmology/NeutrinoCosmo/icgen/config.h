#ifndef CONFIG_H
#define CONFIG_H

// User defined parameters
#define OUTPUT_DIR "output/"
#define GRID_WIDTH 48            // cells
#define BOX_WIDTH 256  // Mpc

// Spectral index and power spectrum normalization
#define N_S 0.954
#define A_S 2.09052e-9
#define PIVOT_SCALE 0.05  // in 1/Mpc

// Present-day fluctuation scale
#define SIGMA_8 0.0 //not used

// Whether to normalize using sigma_8 (NORM_SIGMA) or A_S (NORM_CMB)
#define NORM_SIGMA 0
#define NORM_CMB 1
#define NORMALIZATION_METHOD NORM_CMB

// Initial velocity methods
#define VEL_ZELDOVICH 0
#define VEL_TRANSFER 1
#define VELOCITY_METHOD VEL_TRANSFER

// Whether to set the initial displacements with Zel'dovich or Infinite-LPT
#define FALSE 0
#define TRUE 1
#define USE_INFINI_LPT FALSE

// Whether to infer the velocity of the particles from the original grid/glass
// position (FALSE) or from the velocity field at the displaced position (TRUE)
#define VELOCITY_AT_DISPLACED_POS TRUE

// Whether to output the generalized velocity coordinate p = a^2(dx/dt) or the
// peculiar velocity coordinate v = a (dx/dt) = p/a. HeWon and Swift use p
// internally, but Swift wants v in the initial condition files.
#define GENERALIZED_VELOCITY 0
#define PECULIAR_VELOCITY 1
#define OUTPUT_VELOCITY_FORMAT PECULIAR_VELOCITY

// Whether to output the velocity with a relativistic correction. If yes,
// (dx/dt) is replaced with (dx/ds), i.e. the spatial part of the 4-velocity.
// This differs from the normal output by a factor gamma = 1/sqrt(1-v^2), with
// v=a*(dx/dt) the peculiar velocity in units of c.
// Depending on OUTPUT_VELOCITY_FORMAT, you either get a^2*(dx/ds) or a*(dx/ds)
#define OUTPUT_RELATIVISTIC_VELOCITY TRUE

// Whether to sample the neutrino momenta from a homogenous background model
// or from a linearly perturbed temperature field
#define NU_TEMPERATURE_HOMOGENEOUS 0
#define NU_TEMPERATURE_LINEAR 1
#define NU_TEMPERATURE_MODE NU_TEMPERATURE_LINEAR

// Hubble constant
#define H_0 67.9  // km/s/Mpc

// If particles are generated from a grid, the grid should be size NP^3
#define NP 32
#define PARTICLE_NUM NP* NP* NP

// If neutrino particles are generated from a grid, the grid should be size
// NNUP^3
#define NNUP 1
#define NEUTRINO_NUM NNUP* NNUP* NNUP

// Starting redshift of the simulation
#define Z_START 40;

// Comsological parameters are in cosmo.h

#endif
