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
    ) {

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

/* Initialize flagField, streamField and collideField for passed node */
void initialiseCell(double *collideField, double *streamField, int *flagField, int * n, int * node, int flag) {	
    int i;

    *getFlag(flagField, node, n) = flag;

    for (i = 0; i < Q; ++i){
        *getEl(streamField, node, i, n) = LATTICEWEIGHTS[i];
        *getEl(collideField, node, i, n) = LATTICEWEIGHTS[i];
    }
}

void initialiseFields(double *collideField, double *streamField, int *flagField, int * length, int * my_pos, int* Proc) {
    int x, y, z, i;
    int node[3], walls[6]; /* WALLS : LEFT RIGHT FRONT BACK BOTTOM TOP */

    int n[3] = { length[0] + 2, length[1] + 2, length[2] + 2 };

    /* Determine which type each boundary has */
    for (i = 0; i < D; ++i) {
        if (my_pos[i] == 0)
            walls[2 * i] = NOSLIP;
        else {
            walls[2 * i] = PARALLEL;
        }
        if (my_pos[i] == (Proc[i] - 1))
            walls[2 * i + 1] = NOSLIP;
        else {
            walls[2 * i + 1] = PARALLEL;
        }
    }

    /* Set TOP boundary to MOVING_WALL */
    if(my_pos[2] == (Proc[2] - 1)) {
        walls[5] = MOVING_WALL;
    }

    /* Definition of the fields */
    /* z = 0 */
    node[2] = 0;
    for (y = 0; y < n[1]; ++y){
        node[1] = y;
        for (x = 0; x < n[0]; ++x) {
            node[0] = x;
            initialiseCell(collideField, streamField, flagField, n, node, walls[4]);
        }
    }

    /* Z inner */
    for (z = 1; z <= length[2]; ++z){
        node[2] = z;

        /* y = 0 */
        node[1] = 0;

        for (x = 0; x < n[0]; ++x) {
            node[0] = x;
            initialiseCell(collideField, streamField, flagField, n, node, walls[2]);
        }

        /* Y inner */
        for (y = 1; y <= length[1]; ++y){
            node[1] = y;

            /* x = 0 */
            node[0] = 0;
            initialiseCell(collideField, streamField, flagField, n, node, walls[0]);

            /* X inner */
            for (x = 1; x <= length[0]; ++x){
                node[0] = x;
                initialiseCell(collideField, streamField, flagField, n, node, FLUID);
            }
            
            /* x = length+1 */ 
            node[0] = length[0] + 1;;
            initialiseCell(collideField, streamField, flagField, n, node, walls[1]);
        }
    
        /* y = length+1; */
        node[1] = length[1] + 1;
        for (x = 0; x < n[0]; ++x){
            node[0] = x;
            initialiseCell(collideField, streamField, flagField, n, node, walls[3]);
        }
    }

    /* z = length + 1 */
    node[2] = length[2] + 1;
    for (y = 0; y < n[1]; ++y){
        node[1] = y;
        for (x = 0; x < n[0]; ++x){
            node[0] = x;
            initialiseCell(collideField, streamField, flagField, n, node, walls[5]);
        }
    }
}

/* Allocate memory for sendBuffer and readBuffer */
void initBuffers(double ** readBuffer, double ** sendBuffer, int * length) {
    int n[3] = { length[0] + 2, length[1] + 2, length[2] + 2 };

    /* TODO think about sizes carefully */
    /* LEFT: YZ */
    readBuffer[0] = (double *)calloc(n[1] * n[2] * q, sizeof(double));
    sendBuffer[0] = (double *)calloc(n[1] * n[2] * q, sizeof(double));
    /* RIGHT: YZ */
    readBuffer[1] = (double *)calloc(n[1] * n[2] * q, sizeof(double));
    sendBuffer[1] = (double *)calloc(n[1] * n[2] * q, sizeof(double));

    /* TOP: XY */
    readBuffer[2] = (double *)calloc(n[0] * n[1] * q, sizeof(double));
    sendBuffer[2] = (double *)calloc(n[0] * n[1] * q, sizeof(double));
    /* BOTTOM: XY */
    readBuffer[3] = (double *)calloc(n[0] * n[1] * q, sizeof(double));
    sendBuffer[3] = (double *)calloc(n[0] * n[1] * q, sizeof(double));

    /* FRONT: XZ */
    readBuffer[4] = (double *)calloc(n[0] * n[2] * q, sizeof(double));
    sendBuffer[4] = (double *)calloc(n[0] * n[2] * q, sizeof(double));
    /* BACK: XZ */
    readBuffer[5] = (double *)calloc(n[0] * n[2] * q, sizeof(double));
    sendBuffer[5] = (double *)calloc(n[0] * n[2] * q, sizeof(double));
}
