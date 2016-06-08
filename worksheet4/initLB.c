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

void initialiseCell(double *collideField, double *streamField, int *flagField, int length_tot, int x, int y, int z, int flag) {	
    int i;

    flagField[z * length_tot * length_tot + y * length_tot + x] = flag;

        for (i = 0; i < Q; ++i){
            *getEl(streamField, x, y, z, i, length_tot) = LATTICEWEIGHTS[i];
            *getEl(collideField, x, y, z, i, length_tot) = LATTICEWEIGHTS[i];
        }
}

void initialiseFields(double *collideField, double *streamField, int *flagField, int xlength){
    int x, y, z, length_tot;
    length_tot = xlength + 2;

    /* Definition of the fields */
    z = 0;
    for (y = 0; y < length_tot; ++y){
        for (x = 0; x < length_tot; ++x){			
            initialiseCell(collideField, streamField, flagField, length_tot, x, y, z, NOSLIP);
        }
    }
    for (z = 1; z <= xlength; ++z){
        y = 0;
        for (x = 0; x < length_tot; ++x){
            initialiseCell(collideField, streamField, flagField, length_tot, x, y, z, NOSLIP);
        }

        for (y = 1; y <= xlength; ++y){
            x = 0;
            initialiseCell(collideField, streamField, flagField, length_tot, x, y, z, NOSLIP);
            for (x = 1; x <= xlength; ++x){				
                initialiseCell(collideField, streamField, flagField, length_tot, x, y, z, FLUID);
            }
            x = xlength+1;
            initialiseCell(collideField, streamField, flagField, length_tot, x, y, z, NOSLIP);
        }
    
        y = xlength+1;
        for (x = 0; x < length_tot; ++x){
            initialiseCell(collideField, streamField, flagField, length_tot, x, y, z, NOSLIP);
        }
    }
    z = xlength+1;
    for (y = 0; y < length_tot; ++y){
        for (x = 0; x < length_tot; ++x){
            initialiseCell(collideField, streamField, flagField, length_tot, x, y, z, MOVING_WALL);
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
