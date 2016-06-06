#ifndef _INITLB_H_
#define _INITLB_H_
#include "helper.h"


/* reads the parameters for the lid driven cavity scenario from a config file */
int readParameters(
    int *xlength,                       /* reads domain size. Parameter name: "xlength" */
    double *tau,                        /* relaxation parameter tau. Parameter name: "tau" */
    double *velocityWall,               /* velocity of the lid. Parameter name: "characteristicvelocity" */
    int *timesteps,                     /* number of timesteps. Parameter name: "timesteps" */
    int *timestepsPerPlotting,          /* timesteps between subsequent VTK plots. Parameter name: "vtkoutput" */
    int argc,                           /* number of arguments. Should equal 2 (program + name of config file */
    char *argv[],                       /* argv[1] shall contain the path to the config file */
    int *Proc                           /* Array whit the number of processors per dimention */
);


/* initialises the particle distribution functions and the flagfield */
void initialiseFields(double *collideField, double *streamField,int *flagField, int xlength);

/* Initializes oen cell in the fields, this is called by initialiseFields*/
void initialiseCell(double *collideField, double *streamField, int *flagField, int length_tot, int x, int y, int z, int flag);
#endif

