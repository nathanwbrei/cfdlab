#ifndef _MAIN_C_
#define _MAIN_C_

#include <time.h>

#include "collision.h"
#include "streaming.h"
#include "initLB.h"
#include "visualLB.h"
#include "boundary.h"
#include "LBDefinitions.h"

int main(int argc, char *argv[]){
    int xlength, timesteps, timestepsPerPlotting;
    double tau, velocityWall[3];    
    int t;
    double *swap=NULL;
    clock_t start_time, total_time = 0;

    readParameters(&xlength, &tau, velocityWall, &timesteps, &timestepsPerPlotting, argc, argv);

    double  *collideField = (double *)  malloc((size_t)( Q*(xlength+2)*(xlength+2)*(xlength+2) ) * sizeof( double ));
    double  *streamField  = (double *)  malloc((size_t)( Q*(xlength+2)*(xlength+2)*(xlength+2) ) * sizeof( double ));
    int  *flagField = (int *)  malloc((size_t)( (xlength+2)*(xlength+2)*(xlength+2) ) * sizeof( int ));

    if (collideField == NULL || streamField == NULL || flagField == NULL) {
        ERROR("Unable to allocate matrices.");
    } 


    initialiseFields(collideField, streamField, flagField, xlength);

    treatBoundary(collideField, flagField, velocityWall, xlength);

    for (t = 0; t < timesteps; t++) {

        start_time = clock();  // Start the timer for the lattice updates

        doStreaming(collideField, streamField, flagField, xlength);
        swap = collideField;
        collideField = streamField;
        streamField = swap;
        doCollision(collideField,flagField,&tau,xlength);
        treatBoundary(collideField, flagField, velocityWall, xlength);

        total_time += clock() - start_time; // Add elapsed ticks to total_time

        if (t % timestepsPerPlotting == 0) {
            writeVtkOutput(collideField, flagField, argv[1], t, xlength);
        }
    }

    // Compute average mega-lattice-updates-per-second in order to judge performance
    float elapsed_time = total_time/((float)CLOCKS_PER_SEC);
    float mlups = ((xlength+2) * (xlength+2) * (xlength+2) * timesteps) / (elapsed_time * 1000000);
    printf("Elapsed time (excluding vtk writes) = %f\nAverage MLUPS = %f\n", elapsed_time, mlups);

    free(collideField);
    free(streamField);
    free(flagField);
    
    return 0;
}

#endif

