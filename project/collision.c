#include "helper.h"
#include "collision.h"
#include <omp.h>

float dotProd (const float * array1, float * array2, int length){
    int i;
    float product =0;
    for (i = 0; i < length; ++i){
        product += array1[i] * array2[i];
    }
    return product;
}

void computeNormal(int* node, int* n, float * fractionField, float* normal_dir){
    int node_back[3] = { node[0], node[1], node[2] };
    int node_front[3] = { node[0], node[1], node[2] };

    for (int i = 0; i < D; ++i){
        node_back[i] --;
        node_front[i] ++;
        normal_dir[i] = (*getFraction(fractionField, node_back, n) - *getFraction(fractionField, node_front, n))/2;
        //printf("%e %e , ", *getFraction(fractionField, node_back, n) , *getFraction(fractionField, node_front, n));
        //printf ("%d %d %d , %d %d %d , %d %d %d \n", node_back[0],node_back[1],node_back[2],node[0],node[1],node[2],node_front[0],node_front[1],node_front[2]);
        node_back[i] = node[i];
        node_front[i] = node[i];
    }
    //printf("\n");
}

float computeExternal (int i, float density, float * extForces){
    /* Computes the influence of external forces, like gravity */

    return dotProd(LATTICEVELOCITIES[i], extForces, D) * density * LATTICEWEIGHTS[i];
}

void computePostCollisionDistributions(int *node, float * currentCell, int* flagField, float* fractionField, const float * const tau_inv, const float *const feq, const float *const feqAtm, float density, float * extForces, int * n){
    /* Compute the post-collision distribution f*[i] according to BGK update rule. See Eq 14.  */

    int i;
    float fi;

    if (*getFlag(flagField, node, n) == INTERFACE){
        // int source_node[D];
        float normal_dir[3];
        // float epsilon = 0.00000000001;

        for (i=0; i<Q; i++) {
            /* Make this part with the normak work, folowing eq 4.5 in PHD */
            // source_node[0] = node[0] + LATTICEVELOCITIES[i][0];
            // source_node[1] = node[1] + LATTICEVELOCITIES[i][1];
            // source_node[2] = node[2] + LATTICEVELOCITIES[i][2];
            computeNormal(node, n, fractionField, normal_dir);

            // if (*getFlag(flagField, source_node, n) == GAS || dotProd(LATTICEVELOCITIES[i],normal_dir,D) > epsilon ){
            //     printf("Node %d %d %d, latttice %d %d %d, GasDest %d, normal: %e %e %e , diotprod %e vs %e \n",node[0],node[1],node[2],LATTICEVELOCITIES[i][0],LATTICEVELOCITIES[i][1],LATTICEVELOCITIES[i][2], *getFlag(flagField, source_node, n) ==GAS, normal_dir[0],normal_dir[1],normal_dir[2], dotProd(LATTICEVELOCITIES[i],normal_dir,D), epsilon );
            //     currentCell[Q-i] = feqAtm[i] + feqAtm[Q-i] - currentCell[i];
            // }
            // else{
               fi = currentCell[i];
               currentCell[i] = fi - (fi - feq[i]) * (*tau_inv) + computeExternal(i, density, extForces);
            // }
        }            
    }
    else{
        for (i=0; i<Q; i++) {
            fi = currentCell[i];
            currentCell[i] = fi - (fi - feq[i]) * (*tau_inv) + computeExternal(i, density, extForces);
        }
    }
}

void doCollision(float *collideField, int *flagField, float * massField, float * fractionField, const float * const tau, int * length, float * extForces, int n_threads){
    /* 
     * For each inner grid cell in collideField, compute the post-collide
     * distribution
     */
 
    float density, densityAtm =1;
    float velocity[D];
    float feq[Q], feqAtm[Q];

    float * currentCell;

    int x, y, z, node[3], flag, isFluid;
    int n[3] = { length[0] + 2, length[1] + 2, length[2] + 2 };
    const float tau_inv = 1 / *tau;

#pragma omp parallel for schedule(dynamic) private(x, y, node, density, feq, velocity, densityAtm, feqAtm, isFluid, flag, currentCell) num_threads(n_threads)
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
                    computePostCollisionDistributions(node, currentCell, flagField, fractionField, &tau_inv, feq, feqAtm, density, extForces, n);

                    /* Update fluid fraction */
                    if (isFluid) {
                        *getMass(massField, node, n) = density;
                        *getFraction(fractionField, node, n) = 1; 
                    } else {
                        *getFraction(fractionField, node, n) = *getMass(massField, node, n) / density;
                    }
                }
            }
        }
    }
}

