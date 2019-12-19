/**
 * @sampler.cpp
 * @author  Willem Elbers <whe@willemelbers.com>
 * @version 1.0
 * @license MIT (https://opensource.org/licenses/MIT)
 *
 * @section DESCRIPTION
 *
 * Generates random numbers according to a specified probability density
 * function.
 */

#include "sampler.h"

#include <iostream>
#include <random>
#include <math.h>
#include <chrono>
#include <algorithm>

//The unnormalized pdf (Fermi-Dirac version)
double pdf(double x, double T, double mu) {
    return x*x/(exp((x-mu)/T) + 1);
}

//Numerically integrate the unnormalized pdf to find the cdf and the normalization
double cdf(double x, double T, double mu, int N) {
    //Trapezoid rule integration
    double out = 0;
    double step = x/N;
    for (int i=0; i<N; i++) {
        double mult = (i==0 || i==N-1) ? 0.5 : 1.0;
        out += mult * step * pdf(step*i,T,mu);
    }

    return out;
}


void prepare_intervals(struct sampler *s, double T, double mu = 0) {
    //Setup
    double bl = 0.001;
    double br = 10*T;
    int numIntegrateN = 1000;

    double err = 1e-10;

    //Iteratively determine a good lower bound
    while(pdf(bl,T,mu) > err) {
        bl /= 1.1;
    }

    //Iteratively determine a good upper bound
    while(pdf(br,T,mu) > err) {
        br *= 1.1;
    }

    //Normalize the pdf
    double norm = 1.0/cdf(br, T, mu , 1000);
    s->norm = norm;

    //Create the intervals
    interval first = {0,bl,br,norm*cdf(bl,T,mu,numIntegrateN),norm*cdf(br,T,mu,numIntegrateN)};
    s->intervals.push_back(first);

    int nIntervals = 1;
    interval iv = first;

    bool done = false;
    while(!done) {
        //Check if the interval is too big
        if (iv.Fr - iv.Fl > 0.05) {
            //Split the interval in half
            double m = iv.l + 0.5*(iv.r - iv.l);
            double Fm = norm*cdf(m,T,mu,numIntegrateN);

            //Insert new interval
            int id = nIntervals;
            interval niv = {id, m, iv.r, Fm, iv.Fr};
            niv.nid = iv.nid;
            s->intervals.resize(nIntervals+1);
            s->intervals[niv.id] = niv;
            nIntervals++;

            //Update old interval
            iv.r = m;
            iv.Fr = Fm;
            iv.nid = id;
            s->intervals[iv.id] = iv;
        } else if (iv.r == br) {
            //Stop if we are at the end
            done = true;
        } else {
            //Move on to the next interval
            iv = s->intervals[iv.nid];
        }
    }

    done = false;
    iv = s->intervals[0];
    while(!done) {
        //Evaluate the pdf at the endpoints
        double fl = norm * pdf(iv.l,T,mu);
        double fr = norm * pdf(iv.r,T,mu);

        //Calculate the cubic Hermite approximation
        iv.a0 = iv.l;
        iv.a1 = (iv.Fr - iv.Fl)/fl;
        iv.a2 = 3*(iv.r - iv.l) - (iv.Fr - iv.Fl)*(2./fl + 1./fr);
        iv.a3 = 2*(iv.l - iv.r) + (iv.Fr - iv.Fl)*(1./fl + 1./fr);

        //Evaluate the error at the midpoint
        double u = 0.5*(iv.Fr + iv.Fl);
        double H = iv.a0 + iv.a1*0.5 + iv.a2*pow(0.5,2) + iv.a3*pow(0.5,3);
        iv.error = abs(norm * cdf(H,T,mu,numIntegrateN) - u);
        s->intervals[iv.id] = iv;

        //Monotonicity check
        double delta = (iv.Fr - iv.Fl)/(iv.r - iv.l);
        bool monotonic = (delta <= 3*fl) && (delta <= 3*fr);

        //If the error is too big or if the polynomial is not monotonic
        if (iv.error > 1e-6 || !monotonic) {
            //Split the interval in half
            double m = iv.l + 0.5*(iv.r - iv.l);
            double Fm = norm * cdf(m,T,mu,numIntegrateN);

            //Insert new interval
            int id = nIntervals;
            interval niv = {id, m, iv.r, Fm, iv.Fr};
            niv.nid = iv.nid;
            s->intervals.resize(nIntervals+1);
            s->intervals[niv.id] = niv;
            nIntervals++;

            //Update old interval
            iv.r = m;
            iv.Fr = Fm;
            iv.nid = id;
        } else if (iv.r == br) {
            //Stop if we are at the end
            done = true;
        } else {
            //Move on to the next interval
            iv = s->intervals[iv.nid];
        }
    }

    //The number of intervals
    s->intervalNum = s->intervals.size();
    //Sort the intervals
    std::sort(s->intervals.begin(), s->intervals.end(), compareByLeft);
    //Reset the ids
    for(int i=0; i<s->intervalNum; i++) {
        s->intervals[i].id = i;
    }

    //Allocate memory for the search index
    s->index = (float*) malloc(s->I_max * sizeof(float));

    //Generate the search index array
    for (int i=0; i<s->I_max; i++) {
        float u = (float) i/s->I_max;

        //Find the largest interval such that u > F(p)
        float maxJ = 0;
        int int_id = 0;
        for(auto iv : s->intervals) {
            if (iv.Fr < u && iv.r > maxJ) {
                maxJ = iv.r;
                int_id = iv.id;
            }
        }
        s->index[i] = int_id;
    }
}

double sampler_draw(struct sampler *s) {
    float u = (float) rand() / RAND_MAX;
    const int I_max = s->I_max;
    int I = floor(u*I_max);
    int idx = s->index[I<I_max ? I : I_max-1];

    //Find the largest interval such that u > F(p)
    float maxJ = 0;
    int int_id = idx;
    for(int j=idx; j<s->intervalNum; j++) {
        interval iv = s->intervals[j];
        if (iv.Fr < u && iv.r > maxJ) {
            maxJ = iv.r;
            int_id = iv.id;
        } else {
            break;
        }
    }

    interval iv = s->intervals[int_id];

    float u_tilde = (u - iv.Fl)/(iv.Fr - iv.Fl);
    float H = iv.a0 + iv.a1*u_tilde + iv.a2*u_tilde*u_tilde + iv.a3*u_tilde*u_tilde*u_tilde;

    return H;
}
