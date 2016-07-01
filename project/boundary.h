#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

/** handles the boundaries in our simulation setup */
void treatBoundary(float * collideField,
                   int * flagField,
                   const char * const problem,
                   const double * const Re,
                   const double * const ro_ref,
                   const double * const ro_in,
                   const double * const velocity,
                   int * length);

#endif

