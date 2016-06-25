#ifndef _INITLB_H_
#define _INITLB_H_
#include "helper.h"


/* reads the parameters for the lid driven cavity scenario from a config file */
int readParameters(
    int *length,                       /* reads domain size. Parameter name: "xlength" */
    double *tau,                        /* relaxation parameter tau. Parameter name: "tau" */
    double *velocity,               /* velocity of the lid. Parameter name: "characteristicvelocity" */
    double *extForces,                  /* External force, like gravity or electromagnetic */
    int *timesteps,                     /* number of timesteps. Parameter name: "timesteps" */
    int *timestepsPerPlotting,          /* timesteps between subsequent VTK plots. Parameter name: "vtkoutput" */
    int argc,                           /* number of arguments. Should equal 2 (program + name of config file */
    char *argv[],                       /* argv[1] shall contain the path to the config file */
    char *problem,                      /* specifies the considered problem scenario (parabolic or constant inflow )*/
    double *ro_ref,                     /* reference density nomally set to 1 */
    double *ro_in,                       /* density of inflow/outflow */
    int *boundaries,                     /* definition of the type of boundaries on each one of the walls, for definitions see LBDefinitios.h*/
    int *Proc,                          /* Array whit the number of processors per dimention */
    int my_rank,                         /* Indicates the process number to allow print, if is not pararlel here goes a 0*/
    int * r                             /* radius of the drop */
);


/* initialises the particle distribution functions and the flagfield */
void initialiseFields(double *collideField, double *streamField,int *flagField, double * massField, double * fractionField, int *length, int * boundaries, int r, char *argv[], int * my_pos, int* Proc);

/* Initializes oen cell in the fields, this is called by initialiseFields*/
void initialiseCell(double *collideField, double *streamField, int *flagField, int *n, int * node, int flag);

/* Allocate memory for sendBuffer and readBuffer */
void initBuffers(double ** readBuffer, double ** sendBuffer, int * length);
#endif

