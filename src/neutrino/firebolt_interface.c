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
struct multipoles firebolt_mL; //multipoles in standard Legendre basis

/* SWIFT structure reference */
struct cosmology *c;
struct renderer *r;

/* Redshift as a function of the logarithm of conformal time */
double redshift_func(double log_tau) {
    return rend_perturb_z_at_log_tau(r, log_tau);
}

/* Metric derivative functions */
double h_prime_func(double k, double log_tau) {
    /* Bilinear interpolation indices in (log_tau, k) space */
    int tau_index = 0, k_index = 0;
    double u_tau = 0.f, u_k = 0.f;
    /* Find the time index */
    rend_interp_locate_tau(r, log_tau, &tau_index, &u_tau);
    /* Find the k-space interpolation index */
    rend_interp_locate_k(r, k, &k_index, &u_k);

    /* Interpolate the transfer function */
    double Tr = rend_custom_interp(r, k_index, tau_index, u_tau, u_k,
                                   r->index_transfer_h_prime);
    return Tr;
}

double eta_prime_func(double k, double log_tau) {
    /* Bilinear interpolation indices in (log_tau, k) space */
    int tau_index = 0, k_index = 0;
    double u_tau = 0.f, u_k = 0.f;
    /* Find the time index */
    rend_interp_locate_tau(r, log_tau, &tau_index, &u_tau);
    /* Find the k-space interpolation index */
    rend_interp_locate_k(r, k, &k_index, &u_k);

    /* Interpolate the transfer function */
    double Tr = rend_custom_interp(r, k_index, tau_index, u_tau, u_k,
                                   r->index_transfer_eta_prime);
    return Tr;
}

int firebolt_init(struct swift_params *params, struct renderer *rend, const struct engine *e) {
    c = e->cosmology;
    r = rend;

    /* When is the simulation supposed to start? */
    double a_begin = c->a_begin;
    double tau_begin_sim = cosmology_get_conformal_time(c, a_begin);

    printf("\n");
    printf("Running Firebolt.\n\n");

    /* Read the gaussian random field */
    int N = rend->primordial_grid_N;
    double box_len = rend->primordial_dims[0];

    /* Determine the maximum and minimum wavenumbers */
    double dk = 2*M_PI/box_len;
    double k_max = sqrt(3)*dk*N/2;
    double k_min = dk;

    /* Ensure a safe error margin */
    k_max *= 1.2;
    k_min /= 1.2;

    /* The system to solve */
    double tau_ini = exp(rend->transfer.log_tau[0]);
    double tau_fin = tau_begin_sim; //integrate up to the beginning of the sim
    double a_fin = a_begin; //integrate up to the beginning of the sim

    const struct phys_const *phys_const = e->physical_constants;
    const double kb = phys_const->const_boltzmann_k;
    const double eV = phys_const->const_electron_volt;
    const double M_nu = c->M_nu[0];
    const double M = M_nu * eV / (kb * c->T_nu);
    const double c_vel = e->physical_constants->const_speed_light_c;
    const double rho_crit = c->critical_density;
    const double h_planck = e->physical_constants->const_planck_h;

    /* Size of the problem */
    int l_max = parser_get_param_int(params, "Boltzmann:l_max");
    int k_size = parser_get_param_int(params, "Boltzmann:k_size");
    double q_min = parser_get_param_double(params, "Boltzmann:q_min");
    double q_max = parser_get_param_double(params, "Boltzmann:q_max");
    int q_size = parser_get_param_int(params, "Boltzmann:q_size");
    int verbose = parser_get_opt_param_int(params, "Boltzmann:verbose", 0);
    double tol = parser_get_opt_param_double(params, "Boltzmann:tol", 1e-6);

    rend->boltzmann.P_cdm_prev = malloc(k_size * sizeof(double));
    rend->boltzmann.P_cdm_prime = malloc(k_size * sizeof(double));

    /* Store the momentum range in the renderer struct for later use */
    // rend->firebolt_q_size = q_steps;
    // rend->firebolt_log_q_min = log(q_min);
    // rend->firebolt_log_q_max = log(q_max);

    printf("\n");
    printf("[k_min, k_max, k_size] = [%f, %f, %d]\n", k_min, k_max, k_size);
    printf("[l_max, q_size, q_min, q_max, tol] = [%d, %d, %.1f, %.1f, %.3e]\n", l_max, q_size, q_min, q_max, tol);
    printf("\n");

    printf("The initial time is %f\n", tau_ini);
    printf("Speed of light c = %f\n", c_vel);
    printf("Neutrino mass M_nu = %f (%f eV)\n", M, M_nu);

    /* Initialize the multipoles */
    initMultipoles(&firebolt_mL, k_size, q_size, l_max + 1, q_min, q_max,
                   k_min, k_max);

    /* Store the k dimensions */
    rend->boltzmann.tol = tol;
    rend->boltzmann.k_size = k_size;
    rend->boltzmann.k_min = k_min;
    rend->boltzmann.k_max = k_max;
    rend->boltzmann.verbose = verbose;
    rend->boltzmann.k = firebolt_mL.k;

    /* Calculate the multipoles in Legendre basis */
    evolveMultipoles(&firebolt_mL, tau_ini, tau_fin, tol, M, c_vel,
                     redshift_func, h_prime_func, eta_prime_func, verbose);

    printf("Done with integrating. Processing the moments.\n");

    /* Retrieve the background ncdm density */

    /* Indices in the tau directions */
    int tau_index = 0;
    /* Spacing (0 <= u <= 1) between subsequent indices */
    double u_tau;

    /* Find the index and spacing */
    rend_interp_locate_tau(rend, log(tau_fin), &tau_index, &u_tau);

    double Onu = (1 - u_tau) * rend->transfer.Omegas[rend->transfer.tau_size * rend->index_transfer_delta_ncdm + tau_index]
                     + u_tau * rend->transfer.Omegas[rend->transfer.tau_size * rend->index_transfer_delta_ncdm + tau_index + 1];

    double rho_nu = Onu * rho_crit;
    printf("Onu = %e, rho_crit = %e\n", Onu, rho_nu);

    double degeneracy = 2;
    double factor_ncdm = degeneracy * 4 * M_PI * pow(c->T_nu * kb / (h_planck/(2*M_PI)) / c_vel, 3) * c->T_nu * kb / pow(c_vel, 2);

    /* Calculate the density */
    double *delta_nu = malloc(k_size * sizeof(double));
    computeDensity(&firebolt_mL, delta_nu, M, a_fin, rho_nu, factor_ncdm);
    free(delta_nu);

    return 0;
}

int firebolt_update(struct renderer *rend, const struct engine *e) {

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
    const double rho_crit = c->critical_density;
    const double h_planck = e->physical_constants->const_planck_h;

    /* Size of the problem */
    int k_size = firebolt_mL.k_size;
    double tol = rend->boltzmann.tol;
    int verbose = rend->boltzmann.verbose;

    /* Evolve the Legendre multipoles to the current time */
    evolveMultipoles(&firebolt_mL, tau_ini, tau_fin, tol, M, c_vel,
                     redshift_func, h_prime_func, eta_prime_func, verbose);

    printf("Done with integrating. Processing the moments.\n");

    /* Retrieve the background ncdm density */

    /* Indices in the tau directions */
    int tau_index = 0;
    /* Spacing (0 <= u <= 1) between subsequent indices */
    double u_tau;

    /* Find the index and spacing */
    rend_interp_locate_tau(rend, log(tau_fin), &tau_index, &u_tau);

    double Onu = (1 - u_tau) * rend->transfer.Omegas[rend->transfer.tau_size * rend->index_transfer_delta_ncdm + tau_index]
                     + u_tau * rend->transfer.Omegas[rend->transfer.tau_size * rend->index_transfer_delta_ncdm + tau_index + 1];

    double rho_nu = Onu * rho_crit;
    printf("Onu = %e, rho_crit = %e\n", Onu, rho_nu);

    double degeneracy = 2;
    double factor_ncdm = degeneracy * 4 * M_PI * pow(c->T_nu * kb / (h_planck/(2*M_PI)) / c_vel, 3) * c->T_nu * kb / pow(c_vel, 2);

    printf("%e %e %e %e %e %e\n", factor_ncdm, rho_nu, a_new, c->T_nu * kb, 1.0 / pow(c_vel, 7), pow(c->T_nu * kb / (h_planck/(2*M_PI)) / c_vel, 3) * c->T_nu * kb / pow(c_vel, 4));

    /* Calculate the density */
    double *delta_nu = malloc(k_size * sizeof(double));
    computeDensity(&firebolt_mL, delta_nu, M, a_new, rho_nu, factor_ncdm);
    free(delta_nu);

    return 0;
}

int firebolt_free(void) {

    /* Clean up the Firebolt structures */
    cleanMultipoles(&firebolt_mL);
    free(r->boltzmann.P_cdm_prev);
    free(r->boltzmann.P_cdm_prime);

    return 0;
}
