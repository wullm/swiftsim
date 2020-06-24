/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Willem Elbers (whe@willemelbers.com)
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

/* Config parameters. */
#include "../../config.h"

/* This object's header. */
#include "firebolt_interface.h"

/* Firebolt structures */
struct params firebolt_params;
struct units firebolt_units;
struct perturb_data firebolt_ptdat;

struct multipoles firebolt_mL; //multipoles in standard Legendre basis
struct multipoles firebolt_mmono; //multipoles in monomial basis
struct multipoles firebolt_mgauge; //Legendre multipole gauge transformations
struct grids firebolt_grs;

/* SWIFT structure reference */
struct cosmology *c;

/* Redshift as a function of the logarithm of conformal time */
double z_at_log_tau(double log_tau) {
    double conformal_time = exp(log_tau);
    double a = cosmology_get_scale_factor_from_conformal_time(c, conformal_time);
    return 1./a - 1.;
}

int firebolt_init(struct swift_params *params, struct renderer *rend, const struct engine *e) {
    char settingsFileName[200] = "";
    parser_get_param_string(params, "Boltzmann:firebolt_settings_file",
                            settingsFileName);

    c = e->cosmology;

    /* When is the simulation supposed to start? */
    double a_begin = c->a_begin;
    double tau_begin_sim = cosmology_get_conformal_time(c, a_begin);

    printf("\n");
    printf("Running Firebolt.\n\n");

    /* Read the Firebolt settings and units */
    readParams(&firebolt_params, settingsFileName);
    readUnits(&firebolt_units, settingsFileName);

    /* Read the perturbation vector and initialize the interpolation spline */
    readPerturb(&firebolt_params, &firebolt_units, &firebolt_ptdat);
    initPerturbInterp(&firebolt_ptdat);

    /* Find indices corresponding to specific functions */
    int h_prime_index = 0, eta_prime_index = 0;
    for (int i=0; i<firebolt_ptdat.n_functions; i++) {
      if (strcmp(firebolt_ptdat.titles[i], "h_prime") == 0) {
        h_prime_index = i;
        printf("Found '%s', index = %d\n", firebolt_ptdat.titles[i], i);
      } else if (strcmp(firebolt_ptdat.titles[i], "eta_prime") == 0) {
        eta_prime_index = i;
        printf("Found '%s', index = %d\n", firebolt_ptdat.titles[i], i);
      }
    }

    /* Make sure that the interpolation functions correspond to h' and eta' */
    switchPerturbInterp(&firebolt_ptdat, h_prime_index, 0);
    switchPerturbInterp(&firebolt_ptdat, eta_prime_index, 1);

    /* Read the gaussian random field */
    int N = rend->primordial_box_small_N;
    double box_len = rend->primordial_box_len;

    /* Determine the maximum and minimum wavenumbers */
    double dk = 2*M_PI/box_len;
    double k_max = sqrt(3)*dk*N/2;
    double k_min = dk;

    /* Ensure a safe error margin */
    k_max *= 1.2;
    k_min /= 1.2;

    /* Propagate the grid dimensions */
    firebolt_params.GridSize = N;
    firebolt_params.BoxLen = box_len;

    /* The system to solve */
    double tau_ini = exp(firebolt_ptdat.log_tau[0]);
    double tau_fin = tau_begin_sim; //integrate up to the beginning of the sim
    double a_fin = a_begin; //integrate up to the beginning of the sim

    const struct phys_const *phys_const = e->physical_constants;
    const double kb = phys_const->const_boltzmann_k;
    const double eV = phys_const->const_electron_volt;
    const double M_nu = c->M_nu[0];
    const double M = M_nu * eV / (kb * c->T_nu);
    const double c_vel = e->physical_constants->const_speed_light_c;

    /* Size of the problem */
    int l_max = firebolt_params.MaxMultipole;
    int l_max_convert = firebolt_params.MaxMultipoleConvert;
    int k_size = firebolt_params.NumberWavenumbers;
    double q_min = firebolt_params.MinMomentum;
    double q_max = firebolt_params.MaxMomentum;
    int q_steps = firebolt_params.NumberMomentumBins;
    double tol = firebolt_params.Tolerance;

    /* Store the momentum range in the renderer struct for later use */
    rend->firebolt_q_size = q_steps;
    rend->firebolt_log_q_min = log(q_min);
    rend->firebolt_log_q_max = log(q_max);

    printf("\n");
    printf("[k_min, k_max, k_size] = [%f, %f, %d]\n", k_min, k_max, k_size);
    printf("[l_max, q_steps, q_max, tol] = [%d, %d, %.1f, %.3e]\n", l_max, q_steps, q_max, tol);
    printf("[l_max_convert] = %d\n", l_max_convert);
    printf("\n");

    printf("The initial time is %f\n", exp(firebolt_ptdat.log_tau[0]));
    printf("Speed of light c = %f\n", c_vel);
    printf("Neutrino mass M_nu = %f (%f eV)\n", M, M_nu);

    /* Initialize the multipoles */
    initMultipoles(&firebolt_mL, k_size, q_steps, l_max, q_min, q_max, k_min, k_max);

    /* Also initialize the multipoles in monomial basis (with much lower l_max) */
    initMultipoles(&firebolt_mmono, k_size, q_steps, l_max_convert+1, q_min, q_max, k_min, k_max);

    /* Initialize gauge transforms (only Psi_0 and Psi_1 are gauge dependent) */
    int l_size_gauge = 2;
    initMultipoles(&firebolt_mgauge, k_size, q_steps, l_size_gauge, q_min, q_max, k_min, k_max);

    /* Calculate the multipoles in Legendre basis */
    evolveMultipoles(&firebolt_mL, &firebolt_ptdat, tau_ini, tau_fin, tol, M, c_vel, 0);

    /* Compute the gauge transforms in a separate struct (only Psi_0, Psi_1) */
    convertMultipoleGauge_Nb(&firebolt_mgauge, &firebolt_ptdat, log(tau_fin), a_fin, M, c_vel);

    /* Convert from Legendre basis to monomial basis */
    convertMultipoleBasis_L2m(&firebolt_mL, &firebolt_mmono, l_max_convert);

    /* Convert the gauge transformations to monomial base and add it on top */
    convertMultipoleBasis_L2m(&firebolt_mgauge, &firebolt_mmono, 1);

    printf("Done with integrating. Processing the moments.\n");

    /* Initialize the multipole interpolation splines */
    initMultipoleInterp(&firebolt_mmono);

    /* Generate grids with the monomial multipoles */
    initGrids(&firebolt_params, &firebolt_mmono, &firebolt_grs);
    rend->firebolt_grids = &firebolt_grs;

    generateGrids(&firebolt_params, &firebolt_units, &firebolt_mmono,
                  rend->primordial_phases_small, &firebolt_grs);

    return 0;
}

int firebolt_update(const struct renderer *rend, const struct engine *e) {

    /* Integrate from where to where? */
    double a_old = c->a_old;
    double a_new = c->a;
    double tau_ini = cosmology_get_conformal_time(c, a_old);
    double tau_fin = cosmology_get_conformal_time(c, a_new);

    if (tau_fin == tau_ini) return 0;

    const struct phys_const *phys_const = e->physical_constants;
    const double kb = phys_const->const_boltzmann_k;
    const double eV = phys_const->const_electron_volt;
    const double M_nu = c->M_nu[0];
    const double M = M_nu * eV / (kb * c->T_nu);
    const double c_vel = e->physical_constants->const_speed_light_c;
    const double tol = firebolt_params.Tolerance;

    /* Evolve the Legendre multipoles to the current time */
    evolveMultipoles(&firebolt_mL, &firebolt_ptdat, tau_ini, tau_fin, tol, M, c_vel, 0);

    /* Reset the monomial multipoles and gauge transformations */
    resetMultipoles(&firebolt_mmono);
    resetMultipoles(&firebolt_mgauge);

    /* Compute the gauge transformations (only Psi_0 and Psi_1) */
    convertMultipoleGauge_Nb(&firebolt_mgauge, &firebolt_ptdat, log(tau_fin), a_new, M, c_vel);

    /* Convert Legendre multipoles to monomial basis */
    convertMultipoleBasis_L2m(&firebolt_mL, &firebolt_mmono, firebolt_mmono.l_size - 1);

    /* Convert the gauge transformations to monomial basis and add it on top */
    convertMultipoleBasis_L2m(&firebolt_mgauge, &firebolt_mmono, 1);

    /* Initialize the multipole interpolation splines */
    cleanMultipoleInterp();
    initMultipoleInterp(&firebolt_mmono);

    /* Read the gaussian random field */
    int N = rend->primordial_box_small_N;
    double box_len = firebolt_params.BoxLen;

    /* Generate grids with the monomial multipoles */
    generateGrids(&firebolt_params, &firebolt_units, &firebolt_mmono,
                  rend->primordial_phases_small, &firebolt_grs);

    /* Store the grids on master node */
    if (e->nodeID == 0) {
        for (int index_q=0; index_q<firebolt_mmono.q_size; index_q++) {
            for (int index_l=0; index_l<firebolt_mmono.l_size; index_l++) {
                char boxname[40];
                sprintf(boxname, "grid_l%d_q%d.hdf5", index_l, index_q);
                double *the_grid = firebolt_grs.grids + index_l * (N * N * N)
                                 * firebolt_mmono.q_size + index_q * (N * N * N);
                writeGRF_H5(the_grid, N, box_len, boxname);
            }
        }
    }

    return 0;
}

int firebolt_free(void) {

    /* Clean up the Firebolt structures */
    cleanParams(&firebolt_params);
    cleanPerturbInterp(&firebolt_ptdat);
    cleanPerturb(&firebolt_ptdat);
    cleanMultipoles(&firebolt_mmono);
    cleanMultipoles(&firebolt_mL);
    cleanMultipoles(&firebolt_mgauge);

    cleanMultipoleInterp();

    return 0;
}
