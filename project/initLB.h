#ifndef _INITLB_H_
#define _INITLB_H_
#include "helper.h"


/* reads the parameters for the lid driven cavity scenario from a config file */
int readParameters(
    int *length,                       /* reads domain size. Parameter name: "xlength" */
    float *tau,                        /* relaxation parameter tau. Parameter name: "tau" */
    float *velocity,               /* velocity of the lid. Parameter name: "characteristicvelocity" */
    float *extForces,                  /* External force, like gravity or electromagnetic */
    int *timesteps,                     /* number of timesteps. Parameter name: "timesteps" */
    int *timestepsPerPlotting,          /* timesteps between subsequent VTK plots. Parameter name: "vtkoutput" */
    int argc,                           /* number of arguments. Should equal 2 (program + name of config file */
    char *argv[],                       /* argv[1] shall contain the path to the config file */
    char *problem,                      /* specifies the considered problem scenario (parabolic or constant inflow )*/
    float *ro_ref,                     /* reference density nomally set to 1 */
    float *ro_in,                       /* density of inflow/outflow */
    int *boundaries,                     /* definition of the type of boundaries on each one of the walls, for definitions see LBDefinitios.h*/
    int * r,                             /* radius of the drop */
    int * n_threads, 
    float * exchange
);


/* initialises the particle distribution functions and the flagfield */
void initialiseFields(float *collideField, float *streamField,int *flagField, float * massField, float * fractionField, int *length, int * boundaries, int r, char *argv[], double * num_fluid_cell);

/* Initializes oen cell in the fields, this is called by initialiseFields*/
void initialiseCell(float *collideField, float *streamField, int *flagField, int *n, int * node, int flag);
#endif

