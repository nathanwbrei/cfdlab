#include "helper.h"
#include "collision.h"

void computePostCollisionDistributions(double *currentCell, const double * const tau, const double *const feq){
    /* Compute the post-collision distribution f*[i] according to BGK update rule. See Eq 14.  */

    int i;
    double fi;

    for (i=0; i<19; i++) {
        fi = currentCell[i];
        currentCell[i] = fi - (fi - feq[i]) / *tau;
    }
}

void doCollision(double *collideField, int *flagField,const double * const tau,int xlength){
    /* 
     * For each inner grid cell in collideField, compute the post-collide
     * distribution
     */
 
    double density;
    double velocity[3];
    double feq[19];

    double * currentCell;

    int x,y,z;
    int n = xlength + 2;

    // Loop over inner cells: compare to streaming.c
    for (x=1; x<= xlength; x++) {
        for (y=1; y<=xlength; y++) {
            for (z=1; z<=xlength; z++) {
                currentCell = getEl(collideField, x, y, z, 0, n);
                computeDensity(currentCell, &density);
                computeVelocity(currentCell, &density, velocity);
                computeFeq(&density, velocity, feq);
                computePostCollisionDistributions(currentCell, tau, feq);
            }
        }
    }
    // TODO: Consider iterating over non-flagged boundaries-- How
    // general do we want this?
}

