#ifndef _MAIN_C_
#define _MAIN_C_

#include <time.h>

#include "collision.h"
#include "streaming.h"
#include "initLB.h"
#include "visualLB.h"
#include "boundary.h"
#include "flag.h"
#include "checks.h"

int main(int argc, char *argv[]){
    int length[3], timesteps, timestepsPerPlotting, boundaries[6], r;
    /* TODO do we need inVelocity? */
    double tau, velocity[3], ro_in, ro_ref, extForces[3];
    int t;
    double *swap=NULL;
    clock_t start_time, total_time = 0;
    char problem [10];
    /* Read the file of parameters */
    readParameters(length, &tau, velocity, extForces, &timesteps, &timestepsPerPlotting, argc, argv, problem, &ro_ref, &ro_in, boundaries, &r);
    /* Allocate memory */
    double *collideField = (double *) malloc((size_t)( Q*(length[0]+2)*(length[1]+2)*(length[2]+2)) * sizeof( double ));
    double *streamField  = (double *) malloc((size_t)( Q*(length[0]+2)*(length[1]+2)*(length[2]+2)) * sizeof( double ));
    int *flagField = (int *) malloc((size_t)( (length[0]+2)*(length[1]+2)*(length[2]+2) ) * sizeof( int ));
    double * massField = (double *) malloc((size_t)( (length[0]+2)*(length[1]+2)*(length[2]+2) ) * sizeof(double));
    /* fluid fraction field */
    double * fractionField = (double *) malloc((size_t)( (length[0]+2)*(length[1]+2)*(length[2]+2) ) * sizeof(double));

    double viscosity = C_S * C_S * (tau - 0.5);
    double Re = 1 / viscosity;

    if (collideField == NULL || streamField == NULL || flagField == NULL) {
        ERROR("Unable to allocate matrices.");
    } 
    
    initialiseFields(collideField, streamField, flagField, massField, fractionField, length, boundaries, r, argv);

    // /* TODO This is only for debugging reasons, errase (comment) it when you don't need it */
        // int a=0;
    // for (int z = 0; z <= length[0]+1; ++z)
    // {
    //     for (int y = 0; y <= length[1]+1; ++y)
    //     {
    //         for (int x = 0; x <= length[2]+1; ++x)
    //         {
    //             printf("%d ", flagField[a]);
    //             a++;
    //         }
    //         printf("\n");
    //     }
    //     printf("\n");
    // }

    treatBoundary(collideField, flagField, problem, &Re, &ro_ref, &ro_in, velocity, length);

    for (t = 0; t < timesteps; t++) {

        start_time = clock();  // Start the timer for the lattice updates

        doStreaming(collideField, streamField, flagField, massField, fractionField, length);
        swap = collideField;
        collideField = streamField;
        streamField = swap;
        doCollision(collideField, flagField, massField, fractionField, &tau, length, extForces);
        updateFlagField(collideField, flagField, fractionField, length);
        treatBoundary(collideField, flagField, problem, &Re, &ro_ref, &ro_in, velocity, length);

        total_time += clock() - start_time; // Add elapsed ticks to total_time

        if (t % timestepsPerPlotting == 0) {
            run_checks(collideField, massField, flagField, length, t );
            writeVtkOutput(collideField, flagField, argv[1], t, length);
            printf("Time step %i finished, vtk file was created\n", t);
        }
    }

    // Compute average mega-lattice-updates-per-second in order to judge performance
    float elapsed_time = total_time/((float)CLOCKS_PER_SEC);
    float mlups = ((length[0]+2) * (length[1]+2) * (length[2]+2) * timesteps) / (elapsed_time * 1000000);
    printf("Elapsed time (excluding vtk writes) = %f\nAverage MLUPS = %f\n", elapsed_time, mlups);

    free(collideField);
    free(streamField);
    free(flagField);
    free(massField);
    free(fractionField);
    
    return 0;
}

#endif

