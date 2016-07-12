#ifndef _COLLISION_H_
#define _COLLISION_H_

#include "computeCellValues.h"

/** carries out the whole local collision process. Computes density and velocity and
 *  equilibrium distributions. Carries out BGK update.
 */
void doCollision(float *collideField, int *flagField, float * massField, float * fractionField, const float * const tau,int * length, float * extForces, int n_threads);
#endif
