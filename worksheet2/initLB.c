#include "helper.h"
#include "initLB.h"
#include "LBDefinitions.h"


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
    char *argv[]                        /* argv[1] shall contain the path to the config file */
    ){

    const char *szFileName = argv[1];

    if (argc != 2) {
        ERROR("number of arguments is incorrect");
    }

    read_int(szFileName, "xlength", xlength);
    read_double(szFileName,"tau", tau);
    read_double(szFileName, "characteristicvelocity_x", &velocityWall[0]);
    read_double(szFileName, "characteristicvelocity_y", &velocityWall[1]);
    read_double(szFileName, "characteristicvelocity_z", &velocityWall[2]);
    read_int(szFileName, "timesteps", timesteps);
    read_int(szFileName, "vtkoutput", timestepsPerPlotting);

    return 0;
}

void initialiseCell(double *collideField, double *streamField, int *flagField, int length_tot, int x, int y, int z, int flag){  
>>>>>>> 633c556... Branch ready to merge

    flagField[z * length_tot * length_tot + y * length_tot + x] = flag;

    for (int i = 0; i < 19; ++i){
        *getEl(collideField, x, y, z, i, length_tot) = LATTICEWEIGHTS[i];
    }
}

void initialiseFields(double *collideField, double *streamField, int *flagField, int xlength){

    int x, y, z, length_tot;
    length_tot = xlength +2 ;

    z = 0;
    for (y = 0; y <= xlength+1; ++y){
        for (x = 0; x <= xlength+1; ++x){
            initialiseCell(collideField, streamField, flagField, length_tot, x, y, z, 1);
        }
    }
    for (z = 1; z <= xlength; ++z){
        y = 0;
        for (x = 0; x <= xlength +1; ++x){
            initialiseCell(collideField, streamField, flagField, length_tot, x, y, z, 1);
        }       
        for (y = 1; y <= xlength; ++y){
            x = 0;
            initialiseCell(collideField, streamField, flagField, length_tot, x, y, z, 1);
            for (x = 1; x <= xlength; ++x){
                initialiseCell(collideField, streamField, flagField, length_tot, x, y, z, 0);
            }
            x = xlength+1;
            initialiseCell(collideField, streamField, flagField, length_tot, x, y, z, 1);
        }
        y = xlength+1;
        for (x = 0; x <= xlength+1; ++x){
            initialiseCell(collideField, streamField, flagField, length_tot, x, y, z, 1);
        }
    }
    z = xlength+1;
    for (y = 0; y <= xlength+1; ++y){
        for (x = 0; x <= xlength+1; ++x){
            initialiseCell(collideField, streamField, flagField, length_tot, x, y, z, 2);
        }
    }
}

