#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

/** handles the boundaries in our simulation setup */
void treatBoundary(float * collideField,
                   int * flagField,
                   const char * const problem,
                   const float * const Re,
                   const float * const ro_ref,
                   const float * const ro_in,
                   const float * const velocity,
                   int * length);

#endif

