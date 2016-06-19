#include "helper.h"
#include "collision.h"

double computeExternal (int i, double density, double * extForces){
    /* Computes the influence of external forces, like gravity */
    int j = 0;
    double dprod = 0;
    for (j = 0; j < D; ++j)
        dprod += LATTICEVELOCITIES[i][j] * extForces[j];
    return dprod * density * LATTICEWEIGHTS[i];
}

void computePostCollisionDistributions(double *currentCell, const double * const tau, const double *const feq, double density, double * extForces){
    /* Compute the post-collision distribution f*[i] according to BGK update rule. See Eq 14.  */

    int i;
    double fi;

    for (i=0; i<Q; i++) {
        fi = currentCell[i];
        currentCell[i] = fi - (fi - feq[i]) / *tau + computeExternal(i, density, extForces);
    }
}

void doCollision(double *collideField, int *flagField, double * massField, double * fractionField, const double * const tau, int * length, double * extForces){
    /* 
     * For each inner grid cell in collideField, compute the post-collide
     * distribution
     */
 
    double density;
    double velocity[D];
    double feq[Q];

    double * currentCell;

    int x, y, z, node[3], flag;
    int n[3] = { length[0] + 2, length[1] + 2, length[2] + 2 };

    // Loop over inner cells: compare to streaming.c
    for (z = 1; z <= length[2]; z++) {
        node[2] = z;
        for (y = 1; y <= length[1]; y++) {
            node[1] = y;
            for (x = 1; x <= length[0]; x++) {
                node[0] = x;
                flag = *getFlag(flagField, node, n);
                if (flag == FLUID || flag == INTERFACE || flag == GAS) {
                    currentCell = getEl(collideField, node, 0, n);
                    computeDensity(currentCell, &density);
                    computeVelocity(currentCell, &density, velocity);
                    computeFeq(&density, velocity, feq);
                    computePostCollisionDistributions(currentCell, tau, feq, density, extForces);

                    /* Update fluid fraction */
                    *getFraction(fractionField, node, n) = *getMass(massField, node, n) / density;
                }
            }
        }
    }
}

