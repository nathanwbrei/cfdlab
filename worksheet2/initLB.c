#include "helper.h"
#include "initLB.h"


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


void initialiseFields(double *collideField, double *streamField, int *flagField, int xlength){
  /* TODO */
}

