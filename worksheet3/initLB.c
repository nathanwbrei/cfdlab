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
    double *velocity,               /* velocity of the lid. Parameter name: "characteristicvelocity" */
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

    read_double(szFileName, "velocity_x", &velocity[0]);
    read_double(szFileName, "velocity_y", &velocity[1]);
    read_double(szFileName, "velocity_z", &velocity[2]);
    
    read_int(szFileName, "timesteps", timesteps);
    read_int(szFileName, "vtkoutput", timestepsPerPlotting);
    read_string(szFileName, "problem", problem);
    read_double(szFileName,"ro_ref", ro_ref);
    read_double(szFileName,"ro_in", ro_in);

    read_int(szFileName, "wall_x0", &boundaries[0]);
    read_int(szFileName, "wall_xmax", &boundaries[1]);
    read_int(szFileName, "wall_y0", &boundaries[2]);
    read_int(szFileName, "wall_ymax", &boundaries[3]);
    read_int(szFileName, "wall_z0", &boundaries[4]);
    read_int(szFileName, "wall_zmax", &boundaries[5]);

    if (strcmp(problem, PARABOLIC_SCENARIO) != 0 && strcmp(problem, CONSTANT_SCENARIO) != 0) {
        ERROR("Unrecognized scenario");
    }

    return 0;
}

void initialiseCell(double *collideField, double *streamField, int *flagField, int *length, int x, int y, int z, int flag) {	
    int i;

    int n[3] = { length[0] + 2, length[1] + 2, length[2] + 2 };

    *getFlag(flagField, x, y, z, n) = flag;    /*asigns the flag to the specified cell */

    for (i = 0; i < Q; ++i){
        *getEl(streamField, x, y, z, i, n) = LATTICEWEIGHTS[i];    /*Insert on each cell the initial value */
        *getEl(collideField, x, y, z, i, n) = LATTICEWEIGHTS[i];
    }
}

/**
 * Checks that in input image there was no too thin boundaries
 */
void checkForbiddenPatterns(int ** image, int * length) {
    int x, z;
    int c1, c2, c3, c4, result;

    for (x = 0; x <= length[2]; x++) {
        for (z = 0; z <= length[0]; z++) {
            c1 = image[x][z] == FLUID;
            c2 = image[x + 1][z] == FLUID;
            c3 = image[x][z + 1] == FLUID;
            c4 = image[x + 1][z + 1] == FLUID;
            result = c1 << 3 | c2 << 2 | c3 << 1 | c4;

            if (result == 0b0110 || result == 0b1001) {
                ERROR("forbidden boundary");
            }
        }
    }
}

void initialiseFields(double *collideField, double *streamField, int *flagField, int * length, int * boundaries, char *argv[]){
    int x, y, z;
    int ** image;
    char filename[20] ;

    strcpy(filename,argv[1]);                   /*Copy the value from argv to filename to not modify the original one*/ 
    strcat( filename,".pgm");
    image = read_pgm(filename);/* Concatenate .pgm so the imput is independent of extension*/
    checkForbiddenPatterns(image, length);

    for (z = 0; z <= length[0]; z++) {
        for (x = 0; x <= length[2]; x++) {
            printf("%d ", image[x][z]);
        }
        printf("\n");
    }

    /* 
    *   Definition of the fields 
    *   The structure of the looping makes possible to define the boudaries
    */
    z = 0;
    for (y = 0; y <= length[1]+1; ++y){
        for (x = 0; x <= length[2]+1; ++x){			
            initialiseCell(collideField, streamField, flagField, length, x, y, z, boundaries[4]);   /* Loop for the down wall of the cavity*/
        }
    }
    
    for (z = 1; z <= length[0]; ++z){
        y = 0;
        for (x = 0; x <= length[2]+1; ++x){
            initialiseCell(collideField, streamField, flagField, length, x, y, z, boundaries[2]);   /* Loop for the front wall of the cavity*/
        }

        for (y = 1; y <= length[1]; ++y){
            x = 0;
            initialiseCell(collideField, streamField, flagField, length, x, y, z, boundaries[0]);   /* Loop for the left wall of the cavity*/
            for (x = 1; x <= length[2]; ++x){				
                initialiseCell(collideField, streamField, flagField, length, x, y, z, image[x][z]);           /* Loop for the interior points*/
            }
            x = length[2]+1;
            initialiseCell(collideField, streamField, flagField, length, x, y, z, boundaries[1]);   /* Loop for the right wall of the cavity*/
        }

        y = length[1]+1;
        for (x = 0; x <= length[2]+1; ++x){
            initialiseCell(collideField, streamField, flagField, length, x, y, z, boundaries[3]);   /* Loop for the back wall of the cavity*/
        }
    }

    z = length[0]+1;
    for (y = 0; y <= length[1]+1; ++y){
        for (x = 0; x < length[2]+2; ++x){
            initialiseCell(collideField, streamField, flagField, length, x, y, z, boundaries[5]);   /* Loop for the up wall of the cavity*/
        }
    }
    free_imatrix( image, 0, length[0]+2,0,length[2]+2);
}

