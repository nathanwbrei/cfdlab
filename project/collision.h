#ifndef _COLLISION_H_
#define _COLLISION_H_

#include "computeCellValues.h"

/** computes the post-collision distribution functions according to the BGK update rule and
 *  stores the results again at the same position.
 */
void computePostCollisionDistributions(int *node, float * currentCell, int* flagField, float* fractionField, const float * const tau, const float *const feq, const float *const feqAtm, float density, float * extForces, int * n);

/** carries out the whole local collision process. Computes density and velocity and
 *  equilibrium distributions. Carries out BGK update.
 */
void doCollision(float *collideField, int *flagField, float * massField, float * fractionField, const float * const tau,int * length, float * extForces);
#endif

