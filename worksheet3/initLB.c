#include "helper.h"
#include "initLB.h"
#include "LBDefinitions.h"


/**
 * Reads the parameters for the lid driven cavity scenario from a config file.
 * Throws an error if number of program arguments does not equal 2.
 **/
int readParameters(
    int *length,                        /* reads domain size. Parameter name: "xlength" */
    double *tau,                        /* relaxation parameter tau. Parameter name: "tau" */
    double *velocityWall,               /* velocity of the lid. Parameter name: "characteristicvelocity" */
    int *timesteps,                     /* number of timesteps. Parameter name: "timesteps" */
    int *timestepsPerPlotting,          /* timesteps between subsequent VTK plots. Parameter name: "vtkoutput" */
    int argc,                           /* number of arguments. Should equal 2 (program + name of config file */
    char *argv[],                       /* argv[1] shall contain the path to the config file */
    char *problem,                      /* specifies the considered problem scenario (parabolic or constant inflow )*/
    double *ro_ref,                     /* reference density nomally set to 1 */
    double *ro_in,                      /* density of inflow/outflow */
    int *boundaries                     /* definition of the type of boundaries on each one of the walls, for definitions see LBDefinitios.h*/
    ){

    if (argc != 2) {
        ERROR("number of arguments is incorrect");
    }
    char filename[20] ;
    strcpy(filename,argv[1]);
    const char *szFileName = strcat( filename,".dat"); /* Concatenate .dat so the imput is independent of extesion*/

    read_int(szFileName, "zlength", &length[0]);
    read_int(szFileName, "ylength", &length[1]);
    read_int(szFileName, "xlength", &length[2]);
    read_double(szFileName,"tau", tau);
    read_double(szFileName, "characteristicvelocity_x", &velocityWall[0]);
    read_double(szFileName, "characteristicvelocity_y", &velocityWall[1]);
    read_double(szFileName, "characteristicvelocity_z", &velocityWall[2]);
    read_int(szFileName, "timesteps", timesteps);
    read_int(szFileName, "vtkoutput", timestepsPerPlotting);
    read_string(szFileName, "problem", problem);
    read_double(szFileName,"ro_ref", ro_ref);
    read_double(szFileName,"ro_in", ro_in);
    read_int(szFileName, "wallup", &boundaries[0]);
    read_int(szFileName, "walldown", &boundaries[1]);
    read_int(szFileName, "wallleft", &boundaries[2]);
    read_int(szFileName, "wallright", &boundaries[3]);
    read_int(szFileName, "wallback", &boundaries[4]);
    read_int(szFileName, "wallfront", &boundaries[5]);


    return 0;
}

void initialiseCell(double *collideField, double *streamField, int *flagField, int *length, int x, int y, int z, int flag) {	
    int i;

    flagField[z * (length[0]+2) * (length[0]+2) + y * (length[1]+2) + x] = flag;    /*asigns the flag to the specified cell */

    for (i = 0; i < Q; ++i){
        *getEl(streamField, x, y, z, i, length) = LATTICEWEIGHTS[i];    /*Insert on each cell the initial value */
        *getEl(collideField, x, y, z, i, length) = LATTICEWEIGHTS[i];
    }
}

void initialiseFields(double *collideField, double *streamField, int *flagField, int * length, int * boundaries, char *argv[]){
    int x, y, z;
    int ** image;
    char filename[20] ;    

    strcpy(filename,argv[1]);                   /*Copy the value from argv to filename to not modify the original one*/ 
    strcat( filename,".pgm");
    image = read_pgm(filename);/* Concatenate .pgm so the imput is independent of extesion*/

    /* 
    *   Definition of the fields 
    *   The structure of the looping makes possible to define the boudaries
    */
    z = 0;
    for (y = 0; y <= length[1]+1; ++y){
        for (x = 0; x <= length[2]+1; ++x){			
            initialiseCell(collideField, streamField, flagField, length, x, y, z, boundaries[1]);   /* Loop for the down wall of the cavity*/
        }
    }
    for (z = 1; z <= length[0]; ++z){
        y = 0;
        for (x = 0; x <= length[2]+1; ++x){
            initialiseCell(collideField, streamField, flagField, length, x, y, z, boundaries[5]);   /* Loop for the front wall of the cavity*/
        }
        for (y = 1; y <= length[1]; ++y){
            x = 0;
            initialiseCell(collideField, streamField, flagField, length, x, y, z, boundaries[2]);   /* Loop for the left wall of the cavity*/
            for (x = 1; x <= length[2]; ++x){				
                initialiseCell(collideField, streamField, flagField, length, x, y, z, image[x][z]);           /* Loop for the interior points*/
            }
            x = length[2]+1;
            initialiseCell(collideField, streamField, flagField, length, x, y, z, boundaries[3]);   /* Loop for the right wall of the cavity*/
        }    
        y = length[1]+1;
        for (x = 0; x <= length[2]+1; ++x){
            initialiseCell(collideField, streamField, flagField, length, x, y, z, boundaries[4]);   /* Loop for the back wall of the cavity*/
        }
    }
    z = length[0]+1;
    for (y = 0; y <= length[1]+1; ++y){
        for (x = 0; x < length[2]+2; ++x){
            initialiseCell(collideField, streamField, flagField, length, x, y, z, boundaries[0]);   /* Loop for the up wall of the cavity*/
        }
    }
    free_imatrix( image, 0, length[0]+2,0,length[2]+2);
}
