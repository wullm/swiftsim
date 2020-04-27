/**
 * @cosmo.h
 * @author  Willem Elbers <whe@willemelbers.com>
 * @version 1.0
 * @license MIT (https://opensource.org/licenses/MIT)
 *
 * @section DESCRIPTION
 *
 * Calculate the transfer function per Eisenstein & Hu (1997) for LCDM,
 * as well as other cosmological quantities such as time(redshift) and
 * growth factor D(redshift).
 */

 #ifndef COSMO_H
 #define COSMO_H

#include <math.h>

// Spectral index of primordial fluctuations
const double n_s = N_S;

// Normalization of the CMB power spectrum
const double A_s = A_S;

// Pivot scale in 1/Mpc
const double k_pivot = PIVOT_SCALE;

// Matter power spectrum flucutations on scale 8 Mpc/h
const double sigma_8 = SIGMA_8;

float powerSpectrum(float k) {
  return pow(k, n_s);  // normalization is fixed by sigma_8 in main code, so
                       // A_s/k_star are irrelevant
}

// User-defined paramters
// const double H_0 = 67.556;  // km/s/Mpc (now defined in config.h)
const double h = H_0 / 100;
const double Theta_2p7 = 1;  // CMB T = 2.7 K * Theta_2p7

// Omega paramters
const double Omega_b = 0.0482386599;
const double Omega_m = 0.31346452918;
const double Omega_nu = 0.0013978431;
const double Omega_c = Omega_m - Omega_b;
const double Omega_r = 5.373809e-05;  // this is Omega_g
const double Omega_total = 1.0;
const double Omega_L = Omega_total - Omega_m - Omega_r - Omega_nu;

// Numerical constant [(2/3)*(pi^2/zeta(3))]^(1/3) appearing in neutrino
// temperature-density relation
const double T_nu_const = 1.762359429;
const double T_nu = 1.952784;  // Kelvin (today)
const double mu_nu = 0.0;      // chemical potential

// Neutrino related constants
const double eV = 1.602176634e-19;    // J
const double eV_mass = 1.782662e-36;  // kg (this is eV/c^2)
// Boltzmann's constant in units of eV/K
const double k_b = 8.617333262145e-5;
// The mass of the neutrino
const double M_nu = 0.02;               // eV_mass (per species)
const double M_nu_kg = M_nu * eV_mass;  // kg

// Warm Dark Matter power spectrum, via  Bode, Turok and Ostriker (2001)
float powerSpectrumWDM(float k) {
  float nu = 1.12;
  float alpha = 0.25 / h;
  float damping = pow(1.0 + pow(alpha * k, 2 * nu), -10.0 / nu);
  return pow(k, n_s) * damping;
}

// Calculate matter-radiation equality epoch
const double z_eq = 2.5e4 * Omega_m * h * h * pow(Theta_2p7, -4);
const double k_eq = 7.46e-2 * Omega_m * h * h * pow(Theta_2p7, -2);  // Mpc^-1

// Calculate drag epoch (when the baryons are released from the Compton drag)
const double b_1 = 0.313 * pow(Omega_m * h * h, -0.419) *
                   (1 + 0.607 * pow(Omega_m * h * h, 0.674));
const double b_2 = 0.238 * pow(Omega_m * h * h, 0.223);
const double z_d = 1291 * pow(Omega_m * h * h, 0.251) /
                   (1 + 0.659 * pow(Omega_m * h * h, 0.828)) *
                   (1 + b_1 * pow(Omega_b * h * h, b_2));

// Calculate the sound horizon at the drag epoch
const double R_0 =
    31.5 * Omega_b * h * h * pow(Theta_2p7, -4);  // dimensionless
double R_bpm(double z) {  // The baryon to photon momentum ratio
  return R_0 / (z / 1000);
}
const double R_eq = R_bpm(z_eq);
const double R_d = R_bpm(z_d);
const double s =
    (2.0 / (3 * k_eq)) * sqrt(6.0 / R_eq) *
    log((sqrt(1 + R_d) + sqrt(R_d + R_eq)) / (1 + sqrt(R_eq)));  // Mpc

// Calculate the silk damping scale
const double k_silk = 1.6 * pow(Omega_b * h * h, 0.52) *
                      pow(Omega_m * h * h, 0.73) *
                      (1 + pow(10.4 * Omega_m * h * h, -0.95));  // Mpc^-1

// Function 15 in Eisenstein & Hu
double G(double y) {
  return y * (-6 * sqrt(1 + y) +
              (2 + 3 * y) * log((sqrt(1 + y) + 1) / (sqrt(1 + y) - 1)));
}

// The basic shape ln(k)/k^2 of the transfer function
double T_0_tilde(double k, double alpha_c, double beta_c) {
  double q = k / (13.41 * k_eq);
  double C = 14.2 / alpha_c + 386.0 / (1 + 69.9 * pow(q, 1.08));
  double output =
      log(M_E + 1.8 * beta_c * q) / (log(M_E + 1.8 * beta_c * q) + C * q * q);
  return output;
}

// Spherical Bessel function
double j_0(double x) { return sin(x) / x; }
double sinc(double x) { return sin(x) / x; }

// Now calculate the CDM & baryon transfer function
double Transfer(double k) {
  double q = k / (13.41 * k_eq);

  // CDM transfer function
  double a_1 = pow(46.9 * Omega_m * h * h, 0.670) *
               (1 + pow(32.1 * Omega_m * h * h, -0.532));
  double a_2 = pow(12.0 * Omega_m * h * h, 0.424) *
               (1 + pow(45.0 * Omega_m * h * h, -0.582));
  double alpha_c =
      pow(a_1, -Omega_b / Omega_m) * pow(a_2, -pow(Omega_b / Omega_m, 3));
  double bb_1 = 0.944 / (1 + pow(458 * Omega_m * h * h, -0.708));
  double bb_2 = pow(0.395 * Omega_m * h * h, -0.0266);
  double beta_c = 1.0 / (1 + bb_1 * (pow(Omega_c / Omega_m, bb_2) - 1));
  double f = 1.0 / (1 + pow(k * s / 5.4, 4));  // interpolation factor
  double T_c =
      f * T_0_tilde(k, 1, beta_c) + (1 - f) * T_0_tilde(k, alpha_c, beta_c);

  // Baryon transfer function
  double beta_node = 8.41 * pow(Omega_m * h * h, 0.435);
  double s_tilde = s / pow(1 + pow(beta_node / (k * s), 3),
                           0.33333);  // phenomenological sound horizon
  double alpha_b =
      2.07 * k_eq * s * pow(1 + R_d, -0.75) * G((1 + z_eq) / (1 + z_d));
  double beta_b =
      0.5 + Omega_b / Omega_m +
      (3 - 2 * Omega_b / Omega_m) * sqrt(pow(17.2 * Omega_m * h * h, 2) + 1);
  double T_b =
      (T_0_tilde(k, 1, 1) / (1 + pow(k * s / 5.2, 2)) +
       alpha_b / (1 + pow(beta_b / (k * s), 3)) * exp(-pow(k / k_silk, 1.4))) *
      j_0(k * s_tilde);

  // The final result is a weighted average
  double T = (Omega_b / Omega_m) * T_b + (Omega_c / Omega_m) * T_c;
  return T;
}

// Time in Gyr
double cosmic_time(
    double z) {  // valid for flat LCDM (source: arXiv:gr-qc/0606038v1)
  double conversion_factor = 978.5706;  //(km/s/Mpc)^-1 to Gyr
  return (2 / H_0) / (3 * sqrt(Omega_L)) *
         asinh(sqrt(Omega_L / Omega_m) / pow(z + 1, 1.5)) * conversion_factor;
}

// The scale factor in any cosmology, by definition
double a_scale_factor_of_z(double z) { return 1.0 / (1 + z); }
// The Hubble constant in LCDM, by Friedman's equation
double H_hubble_of_z(double z) {
  double conversion_factor = 978.5706;  //(km/s/Mpc)^-1 to Gyr
  double Hubble = H_0 * sqrt(Omega_r * pow(1 + z, 4) + Omega_m * pow(1 + z, 3) +
                             (1 - Omega_m - Omega_L) * pow(1 + z, 2) + Omega_L);
  return Hubble / conversion_factor;  // in units of Gyr^-1
}

// Approximation for the scale factor in flat LCDM
// (https://physics.stackexchange.com/q/263725)
double a_scale_factor(double t) {       // time in Gyr
  double conversion_factor = 978.5706;  //(km/s/Mpc)^-1 to Gyr
  return pow(
      Omega_m / Omega_L *
          pow(sinh(1.5 * sqrt(Omega_L) * H_0 * t / conversion_factor), 2),
      0.3333);
}
// Approximation for flat LCDM (same source)
double coth(double x) { return cosh(x) / sinh(x); }
double H_hubble(double t) {             // time in Gyr
  double conversion_factor = 978.5706;  //(km/s/Mpc)^-1 to Gyr
  double Hubble = sqrt(Omega_L) * H_0 *
                  coth(1.5 * sqrt(Omega_L) * H_0 * t /
                       conversion_factor);  // in units of (km/s/Mpc)^-1
  return Hubble / conversion_factor;        // in units of Gyr^-1
}

//The fractional contribution of cdm at redshift z
double Omega_cdm_at_z(double z) {
    double conversion_factor = 978.5706;  //(km/s/Mpc)^-1 to Gyr
    double Hubble_ratio = (H_0 / conversion_factor) / H_hubble_of_z(z);
    return Omega_c * pow(1 + z, 3) * pow(Hubble_ratio, 2);
}
//The fractional contribution of baryons at redshift z
double Omega_b_at_z(double z) {
    double conversion_factor = 978.5706;  //(km/s/Mpc)^-1 to Gyr
    double Hubble_ratio = (H_0 / conversion_factor) / H_hubble_of_z(z);
    return Omega_b * pow(1 + z, 3) * pow(Hubble_ratio, 2);
}

double D_growth_factor(double z) {
  double conversion_factor = 978.5706;  //(km/s/Mpc)^-1 to Gyr
  double Hubble_ratio = (H_0 / conversion_factor) / H_hubble_of_z(z);

  double Omega_L_of_z = Omega_L * pow(Hubble_ratio, 2);
  double Omega_m_of_z = Omega_m * pow(1 + z, 3) * pow(Hubble_ratio, 2);
  double D = 2.5 / (1 + z) * Omega_m_of_z /
             (pow(Omega_m_of_z, 4. / 7) - Omega_L_of_z +
              (1 + Omega_m_of_z / 2) * (1 + Omega_L_of_z / 70));
  return D;
}

double logarithmic_derivative_f_1(double z) {
  double conversion_factor = 978.5706;  //(km/s/Mpc)^-1 to Gyr
  double Hubble_ratio = (H_0 / conversion_factor) / H_hubble_of_z(z);

  double Omega_L_of_z = Omega_L * pow(Hubble_ratio, 2);
  double Omega_m_of_z = Omega_m * pow(1 + z, 3) * pow(Hubble_ratio, 2);

  // Constant to a good approximation, see eq. (3) in Jenkins (2009)
  // return pow(Omega_m, 5.0/9.0);
  // Better value given by Hamilton (2001) astro-ph/0006089
  return pow(Omega_m_of_z, 4.0 / 7.0) +
         (1 + Omega_m_of_z / 2) * Omega_L_of_z / 70;
}

double logarithmic_derivative_f_1_simple(double z) {
  double conversion_factor = 978.5706;  //(km/s/Mpc)^-1 to Gyr
  double Hubble_ratio = (H_0 / conversion_factor) / H_hubble_of_z(z);

  double Omega_L_of_z = Omega_L * pow(Hubble_ratio, 2);
  double Omega_m_of_z = Omega_m * pow(1 + z, 3) * pow(Hubble_ratio, 2);

  // Constant to a good approximation, see eq. (3) in Jenkins (2009)
  return pow(Omega_m_of_z, 5.0 / 9.0);
}

// This equals \int_t^{t+\Delta t} dt/a(t).
double kick_step_size(double z1, double z2) {
  double a1 = a_scale_factor_of_z(z1);
  double a2 = a_scale_factor_of_z(z2);
  double a_avg = 0.5 * (a1 + a2);
  double H_avg = 0.5 * (H_hubble_of_z(z1) + H_hubble_of_z(z2));

  // Use dt = da/(aH)
  double delta_t = (a2 - a1) / (a_avg * H_avg);
  double inv_a_avg = 0.5 * (pow(a1, -1) + pow(a2, -1));

  return inv_a_avg * delta_t;

  // Below is an analytical value, but I think it may be incorrect
  // See appendix A of Quinn et al. (1997) astro-ph/9710043 for the kick
  // operator K(t).
  // double conversion_factor = 978.5706; //(km/s/Mpc)^-1 to Gyr
  // return 2.0/(H_0/conversion_factor) * (sqrt(a_scale_factor_of_z(z2)) -
  // sqrt(a_scale_factor_of_z(z1)));
}
// This equals \int_t^{t+\Delta t} dt/a(t)^2.
double drift_step_size(double z1, double z2) {
  double a1 = a_scale_factor_of_z(z1);
  double a2 = a_scale_factor_of_z(z2);
  double a_avg = 0.5 * (a1 + a2);
  double H_avg = 0.5 * (H_hubble_of_z(z1) + H_hubble_of_z(z2));

  // Use dt = da/(aH)
  double delta_t = (a2 - a1) / (a_avg * H_avg);
  double inv_aa_avg = 0.5 * (pow(a1, -2) + pow(a2, -2));

  return inv_aa_avg * delta_t;

  // Below is an analytical value, but I think it may be incorrect
  // See appendix A of Quinn et al. (1997) astro-ph/9710043 for the drift
  // operator D(t).
  // double conversion_factor = 978.5706; //(km/s/Mpc)^-1 to Gyr
  // return 2.0/(H_0/conversion_factor) * (1.0/sqrt(a_scale_factor_of_z(z1))
  // - 1.0/sqrt(a_scale_factor_of_z(z2)));
}

// The main units used in the HeWon simulation code
const double Mpc = 3.086e22;  // m
const double Gyr = 3.154e16;  // s
const double km = 1000;       // m
const double kg = 1000;       // g
const double cm = 0.01;       // m

//Mass unit used by Swift is 10^10 M_sol
const double M_sol = 1.98848e+33; //g
const double M_swift = 1e10 * M_sol; //g

// The speed of light in units of Mpc/Gyr
const double c_vel = 299792458 * Gyr / Mpc;

// One unit of critical density
const double H_0_Hz = H_0 / (Mpc / km);  // s^-1
const double G_newton = 6.67430e-11;     // m^3/kg/s^2
const double rho_crit = 3 * pow(H_0_Hz, 2) / (8 * M_PI * G_newton);  // kg m^-3

// Initial WDM velocity dispersion
double rms_vel_km = 1000;                 // m/sf
double rms_vel = rms_vel_km / Mpc * Gyr;  // Mpc/Gyr

// A test of the various functions introduced in cosmo.h
// Compare e.g. with Weinberg et al. (2013) arXiv:1201.2434
inline void test_cosmology(std::string fname) {
  std::ofstream of(fname);

  of << "t z a_1 a_2 H_1 H_2 D f_1 f_2" << std::endl;
  double z_scaling = 1.01;
  double z_start = Z_START;
  double z_rs = z_start * z_scaling;
  double z_penultimate = 0.001;
  while (z_rs >= z_penultimate) {
    // Jump down in redshift (include z=0 as final row)
    z_rs = (z_rs > z_penultimate * z_scaling) ? z_rs / z_scaling : 0;

    // Test our methods
    double t = cosmic_time(z_rs);
    double a1 = a_scale_factor(t);
    double a2 = a_scale_factor_of_z(z_rs);
    double H1 = H_hubble(t);
    double H2 = H_hubble_of_z(z_rs);
    double D = D_growth_factor(z_rs);
    double f1 = logarithmic_derivative_f_1(z_rs);
    double f2 = logarithmic_derivative_f_1_simple(z_rs);

    // Output a row
    of << t << " " << z_rs << " " << a1 << " " << a2 << " " << H1 << " " << H2
       << " " << D << " " << f1 << " " << f2 << std::endl;
  }

  of.close();

  if (Omega_total != 1) {
    std::cout << "Warning: Omega != 1, while this is assumed for some methods."
              << std::endl;
  }
}

#endif
