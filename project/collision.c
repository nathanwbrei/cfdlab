#include "helper.h"
#include "collision.h"

double dotProd (const int * array1, double * array2, int length){
    int i;
    double product =0;
    for (i = 0; i < length; ++i){
        product += array1[i] * array2[i];
    }
    return product;
}

void computeNormal(int* node, int* length, double * fractionField, double* normal_dir){
    int n[3] = { length[0] + 2, length[1] + 2, length[2] + 2 };
    int node_back[3] = { node[0], node[1], node[2] };
    int node_front[3] = { node[0], node[1], node[2] };

    for (int i = 0; i < D; ++i){
        node_back[i] --;
        node_front[i] ++;
        normal_dir[i] = (*getFraction(fractionField, node_back, n) - *getFraction(fractionField, node_front, n))/2;
        node_back[i] = node[i];
        node_front[i] = node[i];
    }
}

double computeExternal (int i, double density, double * extForces){
    /* Computes the influence of external forces, like gravity */

    return dotProd(LATTICEVELOCITIES[i], extForces, D) * density * LATTICEWEIGHTS[i];
}

void computePostCollisionDistributions(int *node, double * currentCell, int* flagField, double* fractionField, const double * const tau, const double *const feq, const double *const feqAtm, double density, double * extForces, int * length){
    /* Compute the post-collision distribution f*[i] according to BGK update rule. See Eq 14.  */

    int i;
    double fi;
    // int n[3] = { length[0] + 2, length[1] + 2, length[2] + 2 };

    // if (*getFlag(flagField, node, n) == INTERFACE){
    //     int source_node[D];
    //     double normal_dir[3];

    //     for (i=0; i<Q; i++) {
    //         source_node[0] = node[0] - LATTICEVELOCITIES[i][0];
    //         source_node[1] = node[1] - LATTICEVELOCITIES[i][1];
    //         source_node[2] = node[2] - LATTICEVELOCITIES[i][2];
    //         computeNormal(node, length, fractionField, normal_dir);

    //         if (*getFlag(flagField, source_node, n) == GAS || dotProd(LATTICEVELOCITIES[i],normal_dir,D) >0 ){
    //             currentCell[Q-i] = feqAtm[i] + feqAtm[Q-i] - currentCell[i];
    //         }
    //         else{
    //             fi = currentCell[i];
    //             currentCell[i] = fi - (fi - feq[i]) / *tau + computeExternal(i, density, extForces);
    //         }
    //     }            
    // }
    // else{
        for (i=0; i<Q; i++) {
            fi = currentCell[i];
            currentCell[i] = fi - (fi - feq[i]) / *tau + computeExternal(i, density, extForces);
        }
    // }
}

void doCollision(double *collideField, int *flagField, double * massField, double * fractionField, const double * const tau, int * length, double * extForces){
    /* 
     * For each inner grid cell in collideField, compute the post-collide
     * distribution
     */
 
    double density, densityAtm =1;
    double velocity[D];
    double feq[Q], feqAtm[Q];

    double * currentCell;

    int x, y, z, node[3], flag, isFluid;
    int n[3] = { length[0] + 2, length[1] + 2, length[2] + 2 };

    // Loop over inner cells: compare to streaming.c
    for (z = 1; z <= length[2]; z++) {
        node[2] = z;
        for (y = 1; y <= length[1]; y++) {
            node[1] = y;
            for (x = 1; x <= length[0]; x++) {
                node[0] = x;
                flag = *getFlag(flagField, node, n);
                isFluid = flag == FLUID;
                if (isFluid || flag == INTERFACE) {
                    currentCell = getEl(collideField, node, 0, n);
                    computeDensity(currentCell, &density);
                    computeVelocity(currentCell, &density, velocity);
                    computeFeq(&density, velocity, feq);
                    computeFeq(&densityAtm, velocity, feqAtm);
                    computePostCollisionDistributions(node, currentCell, flagField, fractionField, tau, feq, feqAtm, density, extForces, n);

                    /* Update fluid fraction */
                    *getFraction(fractionField, node, n) = *getMass(massField, node, n) / density;
                    if (isFluid) {
                        *getMass(massField, node, n) = density;
                    }
                }
            }
        }
    }
}

