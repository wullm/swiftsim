#ifndef CONFIG_H
#define CONFIG_H

// User defined parameters
#define OUTPUT_DIR "output/"
#define GRID_WIDTH 64  // cells
#define BOX_WIDTH 64   // Mpc

// If particles are generated from a grid, the grid should be size NP^3
#define NP 32
#define PARTICLE_NUM NP* NP* NP

#define NEUTRINO_NUM 10

// Starting redshift of the simulation
#define Z_START 40;

// Comsological parameters are in cosmo.h

#endif
