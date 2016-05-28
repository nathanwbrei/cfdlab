#include "helper.h"
#include "collision.h"

void computePostCollisionDistributions(double *currentCell, const double * const tau, const double *const feq){
    /* Compute the post-collision distribution f*[i] according to BGK update rule. See Eq 14.  */

    int i;
    double fi;

    for (i=0; i<Q; i++) {
        fi = currentCell[i];
        currentCell[i] = fi - (fi - feq[i]) / *tau;
    }
}

void doCollision(double *collideField, int *flagField, const double * const tau, int * length){
    /* 
     * For each inner grid cell in collideField, compute the post-collide
     * distribution
     */
 
    double density;
    double velocity[D];
    double feq[Q];

    double * currentCell;

    int x,y,z;
    int n[3] = { length[0] + 2, length[1] + 2, length[2] + 2 };

    // Loop over inner cells: compare to streaming.c
    for (z=1; z<=length[0]; z++) {
        for (y=1; y<=length[1]; y++) {
            for (x=1; x<=length[2]; x++) {
                currentCell = getEl(collideField, x, y, z, 0, n);
                computeDensity(currentCell, &density);
                computeVelocity(currentCell, &density, velocity);
                computeFeq(&density, velocity, feq);
                computePostCollisionDistributions(currentCell, tau, feq);
            }
        }
    }
}

