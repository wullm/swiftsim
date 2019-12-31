#ifndef CONFIG_H
#define CONFIG_H

// User defined parameters
#define OUTPUT_DIR "output/"
#define GRID_WIDTH 64  // cells
#define BOX_WIDTH 1000   // Mpc

// Spectral index
#define N_S 0.9619

// Hubble constant
#define H_0 67.556 // km/s/Mpc

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
