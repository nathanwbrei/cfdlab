#include "helper.h"
#include "initLB.h"
#include "LBDefinitions.h"
#include "mpi.h"

/**
 * Reads the parameters for the lid driven cavity scenario from a config file.
 * Throws an error if number of program arguments does not equal 2.
 **/
int readParameters(
    int *xlength,                       /* reads domain size. Parameter name: "xlength" */
    double *tau,                        /* relaxation parameter tau. Parameter name: "tau" */
    double *velocityWall,               /* velocity of the lid. Parameter name: "characteristicvelocity" */
    int *timesteps,                     /* number of timesteps. Parameter name: "timesteps" */
    int *timestepsPerPlotting,          /* timesteps between subsequent VTK plots. Parameter name: "vtkoutput" */
    int argc,                           /* number of arguments. Should equal 2 (program + name of config file */
    char *argv[],                       /* argv[1] shall contain the path to the config file */
    int *Proc,                          /* Array whit the number of processors per dimention */
    int my_rank                         /* Indicates the process number to allow print, if is not pararlel here goes a 0*/
    ){

    const char *szFileName = argv[1];

    if (argc != 2) {
        ERROR("number of arguments is incorrect");
    }

    read_int(szFileName, "xlength", xlength, my_rank);
    read_double(szFileName,"tau", tau, my_rank);
    read_double(szFileName, "characteristicvelocity_x", &velocityWall[0], my_rank);
    read_double(szFileName, "characteristicvelocity_y", &velocityWall[1], my_rank);
    read_double(szFileName, "characteristicvelocity_z", &velocityWall[2], my_rank);
    read_int(szFileName, "timesteps", timesteps, my_rank);
    read_int(szFileName, "vtkoutput", timestepsPerPlotting, my_rank);
    read_int(szFileName, "iProc", &Proc[0], my_rank);
    read_int(szFileName, "jProc", &Proc[1], my_rank);
    read_int(szFileName, "kProc", &Proc[2], my_rank);

    return 0;
}

void initialiseCell(double *collideField, double *streamField, int *flagField, int * n, int * node, int flag) {	
    int i;

    *getFlag(flagField, node, n) = flag;

    for (i = 0; i < Q; ++i){
        *getEl(streamField, node, i, n) = LATTICEWEIGHTS[i];
        *getEl(collideField, node, i, n) = LATTICEWEIGHTS[i];
    }
}

void initialiseFields(double *collideField, double *streamField, int *flagField, int * length){
    int x, y, z;
    int node[3];

    int n[3] = { length[0] + 2, length[1] + 2, length[2] + 2 };

    /* Definition of the fields */
    /* z = 0 */
    node[2] = 0;
    for (y = 0; y < n[1]; ++y){
        node[1] = y;
        for (x = 0; x < n[0]; ++x) {
            node[0] = x;
            initialiseCell(collideField, streamField, flagField, n, node, NOSLIP);
        }
    }

    for (z = 1; z <= length[2]; ++z){
        node[2] = z;

        /* y = 0 */
        node[1] = 0;

        for (x = 0; x < n[0]; ++x) {
            node[0] = x;
            initialiseCell(collideField, streamField, flagField, n, node, NOSLIP);
        }

        for (y = 1; y <= length; ++y){
            node[1] = y;

            /* x = 0 */
            node[0] = 0;
            initialiseCell(collideField, streamField, flagField, n, node, NOSLIP);

            for (x = 1; x <= length[0]; ++x){
                node[0] = x;
                initialiseCell(collideField, streamField, flagField, n, node, FLUID);
            }
            
            /* x = length+1 */
            node[0] = length[0] + 1;;
            initialiseCell(collideField, streamField, flagField, n, node, NOSLIP);
        }
    
        /* y = length+1; */
        node[1] = length[1] + 1;
        for (x = 0; x < n[0]; ++x){
            node[0] = x;
            initialiseCell(collideField, streamField, flagField, n, node, NOSLIP);
        }
    }

    /* z = length + 1 */
    node[2] = length[2] + 1;
    for (y = 0; y < n[1]; ++y){
        node[1] = y;
        for (x = 0; x < n[0]; ++x){
            node[0] = x;
            initialiseCell(collideField, streamField, flagField, n, node, MOVING_WALL);
        }
    }
}

void get_rank_pos(int * my_pos, int rank, int *Proc){
    int i;
    int aux_rank = rank;
    for ( i = 0; i < D; ++i){
        my_pos[i] = aux_rank % Proc[i];
        aux_rank = (aux_rank- my_pos[i]) / Proc[i] ;
    }
}

void get_my_lengths(int* my_pos, int xlength, int* my_lengths, int * Proc){
    int i;
    for (i = 0; i < D; ++i){
        my_lengths[i] = (int) (1.0*xlength / Proc[i]);  /*Divides the length of the cavity betwen the number of sections in that dimention (takes the floor)*/
        if ((xlength % Proc[i]) > my_pos[i])    /*If the number of sections is not a divisor of the number of cells then an extra cell is assigned to the first processes*/
            my_lengths[i] ++;
    }
}
