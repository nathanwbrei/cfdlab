#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

/** Function treatBoundary handles the boundaries in the simulation setup for the lattice Boltzmann model
    This function has no return values and it writes directly in the fields */
void treatBoundary(float * collideField,            /*Pointer to the start of the collide field, with all the lattices */
                   int * flagField,                 /*Pointer to the start of the collide field*/
                   const char * const scenatio,     /*Scenario used to indicate parabolic inflow */
                   const float * const Re,          /*Is the inverse of the viscosity */
                   const float * const ro_ref,      /*Reference density */
                   const float * const ro_in,       /*Inflow density */
                   const float * const velocity,    /*Vector with the velocity of the flow*/
                   int * length,                    /*Length of the cavity on each dimension*/
                   int n_threads);                  /*Number of threads available to parallelize the code*/

#endif