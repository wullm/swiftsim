/**
 * @read_transfer.cpp
 * @author  Willem Elbers <whe@willemelbers.com>
 * @version 1.0
 * @license MIT (https://opensource.org/licenses/MIT)
 *
 * @section DESCRIPTION
 *
 * Read CLASS format transfer function files.
 */

#include "read_transfer.h"

void read_transfer(std::vector<double>& ks, std::vector<double>& T_rho, std::vector<double>& T_theta) {
    std::ifstream f("transfer/class_transfer_z40.dat");

    //Ignore lines starting with # until we find a line that doesn't
    std::string dummyLine;
    bool done = false;
    while (!done) {
        std::getline(f, dummyLine);
        if (dummyLine[0] != '#') {
            done = true;
        }
    }

    //Number of columns in the CLASS transfer function output file
    int ncol = 11;
    //Column that contains the wavenumber k (h/Mpc)
    int kcol = 0;
    //Column that contains the cdm transfer function d_cdm
    int d_cdm_col = 3;

    //Get the values of k (h/Mpc) in column 0
    // std::vector<double> ks;
    //Get the neutrino density transfer function from column 5
    // std::vector<double> T_rho_nu;
    //Get the neutrino theta (velocity divergence) transfer function from column 16
    // std::vector<double> T_theta_nu;

    int i = 0;
    double dbl;
    while( f >> dbl ) {
        if (i%ncol == kcol) {
            ks.push_back(dbl);
        } else if (i%ncol == d_cdm_col) {
            T_rho.push_back(dbl);
        }
        // else if (i%20 == 16) {
        //     T_theta_nu.push_back(dbl);
        // }
        i++;
    }

    //Number of rows
    int kvals = ks.size();

    //Multiply the transfer function by -1/k^2 to get the same format as
    //the Eisenstein-Hu transfer functions (also the format of CAMB)
    for (int j=0; j<kvals; j++) {
        T_rho[j] *= -pow(ks[j], -2);
    }

    //Normalize the power spectrum such that T(k) -> 1 as k -> 0
    double amplitude_k_min = T_rho[0]; //Amplitude at the largest scale
    for (int j=0; j<kvals; j++) {
        T_rho[j] /= amplitude_k_min;
    }


    // for (int i = 0; i < ks.size(); i++) {
    // 	std::cout << ks[i] << " " << T_rho[i] << " " << T_theta_nu[i] << std::endl;
    // }
}
