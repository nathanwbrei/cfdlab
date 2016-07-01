#ifndef _COLLISION_H_
#define _COLLISION_H_

#include "computeCellValues.h"

/** computes the post-collision distribution functions according to the BGK update rule and
 *  stores the results again at the same position.
 */
void computePostCollisionDistributions(int *node, float * currentCell, int* flagField, double* fractionField, const double * const tau, const float *const feq, const float *const feqAtm, float density, double * extForces, int * n);

/** carries out the whole local collision process. Computes density and velocity and
 *  equilibrium distributions. Carries out BGK update.
 */
void doCollision(float *collideField, int *flagField, double * massField, double * fractionField, const double * const tau,int * length, double * extForces);
#endif

