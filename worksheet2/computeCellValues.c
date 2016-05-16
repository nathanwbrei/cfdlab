#include "LBDefinitions.h"
#include "computeCellValues.h"
#include <stdio.h>

void computeDensity(const double *const currentCell, double *density){
    /*
     * Computes the macroscopic density within currentCell
     * rho(x,t) = sum(f_i, i=0:Q-1)
     */

    int i;

    *density = 0;
    for (i=0; i<Q; i++) {
        *density += currentCell[i];
    }
}

void computeVelocity(const double * const currentCell, const double * const density, double *velocity){
    /*
     * Computes the velocity within currentCell 
     * u(x,t) = sum(f[i] * c[i] for i in [0:Q-1]) / rho
     */

    int i;

    velocity[0] = 0;
    velocity[1] = 0;
    velocity[2] = 0;

    for (i=0; i<Q; i++) {
        velocity[0] += LATTICEVELOCITIES[i][0] * currentCell[i];
        velocity[1] += LATTICEVELOCITIES[i][1] * currentCell[i];
        velocity[2] += LATTICEVELOCITIES[i][2] * currentCell[i];
    }

    velocity[0] /= *density;
    velocity[1] /= *density;
    velocity[2] /= *density;
}


static double cs2 = C_S * C_S;
void computeFeq(const double * const density, const double * const velocity, double *feq){
    /*
     * feq[i] = w[i] * rho * (1 + c_i*u / (c_s^2) + ... 
     */

    int i;
    double ci0, ci1, ci2, u0, u1, u2;

    u0 = velocity[0];
    u1 = velocity[1];
    u2 = velocity[2];

    double t1 = (u0*u0 + u1*u1 + u2*u2 ) / (2 * cs2);

    for (i=0; i<Q; i++) {
        
        ci0 = LATTICEVELOCITIES[i][0];
        ci1 = LATTICEVELOCITIES[i][1];
        ci2 = LATTICEVELOCITIES[i][2];

        double t2 = (ci0*u0 + ci1*u1 + ci2*u2) / cs2;

        feq[i] = (1 + t2 + (t2 * t2 / 2) - t1) * (*density) * LATTICEWEIGHTS[i];

    }
}



