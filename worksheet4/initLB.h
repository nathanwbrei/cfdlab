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
    int *Proc,                          /* Array whit the number of processors per dimention */
    int my_rank
);


/* initialises the particle distribution functions and the flagfield */
void initialiseFields(double *collideField, double *streamField,int *flagField, int * length, int * my_pos, int * Proc);

/* Initializes oen cell in the fields, this is called by initialiseFields*/
void initialiseCell(double *collideField, double *streamField, int *flagField, int * n, int * node, int flag);

/* Obtain the position in the D dimentional space my_pos, given the process number rank amd the number of process per dimention (Proc)*/
void get_rank_pos(int * my_pos, int rank, int *Proc);

/* Obtain the size of the prtion assigned to the procces my_lengths, using the length of the cavity xlength, the position on each dimetion y the total number of processes per dimention*/
void get_my_lengths(int* my_pos, int xlength, int* my_lengths, int * Proc);

#endif

