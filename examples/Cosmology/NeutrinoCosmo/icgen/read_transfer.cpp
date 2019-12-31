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
#include "config.h"

void read_transfer(std::vector<double>& ks, std::vector<double>& T_rho, std::vector<double>& T_rho_nu) {
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
    int ncol = 15;
    //Column that contains the wavenumber k (h/Mpc)
    int kcol = 0;
    //Column that contains the cdm transfer function d_cdm
    int d_cdm_col = 3;
    //Column that contains the neutrino transfer function d_cdm
    int d_ncdm_col = 5;

    int i = 0;
    double dbl;
    while( f >> dbl ) {
        if (i%ncol == kcol) {
            ks.push_back(dbl);
        } else if (i%ncol == d_cdm_col) {
            T_rho.push_back(dbl);
        } else if (i%ncol == d_ncdm_col) {
            T_rho_nu.push_back(dbl);
        }
        i++;
    }

    //Number of rows
    int kvals = ks.size();

    //Multiply the transfer functions by -1/k^2 to get the same format as
    //the Eisenstein-Hu transfer functions (also the format of CAMB)
    for (int j=0; j<kvals; j++) {
        T_rho[j] *= -pow(ks[j], -2);
        T_rho_nu[j] *= -pow(ks[j], -2);
    }

    //Normalize the power spectrum such that T(k) -> 1 as k -> 0
    double amplitude_k_min = T_rho[0]; //Amplitude of CDM at the largest scale
    for (int j=0; j<kvals; j++) {
        T_rho[j] /= amplitude_k_min;
        T_rho_nu[j] /= amplitude_k_min;
    }

    //Convert from h/Mpc to 1/Mpc
    for (int j=0; j<kvals; j++) {
        ks[j] /= (H_0 / 100.0);
    }
}

//Make an indexed search table for the transfer function interpolation
void make_index_table(int TF_I_max, float *TF_index, std::vector<double>& TF_ks, float *log_k_min, float *log_k_max) {
    //The smallest and largest k-values in the transfer fucntion tables
    float k_min = TF_ks[0];
    float k_max = TF_ks[TF_ks.size() - 1];
    *log_k_min = log(k_min);
    *log_k_max = log(k_max);

    //Make the index table
    for (int i=0; i<TF_I_max; i++) {
        float u = (float) i/TF_I_max;
        float v = *log_k_min + u * (*log_k_max - *log_k_min);
        float w = exp(v);

        //Find the largest bin such that w > k
        float maxJ = 0;
        int int_id = 0;
        for(int j=0; j<TF_ks.size(); j++) {
            if (TF_ks[j] < w) {
                maxJ = j;
            }
        }
        TF_index[i] = maxJ;
    }
}
