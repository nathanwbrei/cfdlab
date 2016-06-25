#ifndef _MAIN_C_
#define _MAIN_C_

#include <time.h>
#include <mpi.h>

#include "parallel.h"
#include "collision.h"
#include "streaming.h"
#include "initLB.h"
#include "visualLB.h"
#include "boundary.h"
#include "flag.h"
#include "checks.h"

int main(int argc, char *argv[]){
    int length[3], timesteps, timestepsPerPlotting, boundaries[6], r, t;
    
    /* Parallel parameters */
    int my_rank, number_of_ranks, Proc[3], my_pos[3], my_lengths[3], my_origin[3];

    /* TODO do we need inVelocity? */
    double tau, velocity[3], ro_in, ro_ref, extForces[3];
    double *swap=NULL;
    
    double * sendBuffer[6]; /* [0:left,1:right,2:top,3:bottom,4:front,5:back] */
    double * readBuffer[6];

    clock_t start_time, total_time = 0;

    char problem [10];

    /* Initialize MPI */
    initializeMPI(&my_rank, &number_of_ranks, argc, argv);

    /* Read the file of parameters */
    readParameters(length, &tau, velocity, extForces, &timesteps, &timestepsPerPlotting, argc, argv, problem, &ro_ref, &ro_in, boundaries, Proc, my_rank, &r);
    
    /* If the number of process in the imput file and the arguments does not coicide, launch an error */
    if (number_of_ranks != (Proc[0]*Proc[1]*Proc[2])){
        printf("Proc %d : Number of processes does not match %d vs %d \n", my_rank, number_of_ranks,(Proc[0]*Proc[1]*Proc[2]));
    }

    /* Get position in the topology */
    // TODO MPI topology could be used
    get_rank_pos(my_pos, my_rank, Proc);

    /* Get local lengths */
    get_my_lengths(my_pos, length, Proc, my_lengths, my_origin);

    initBuffers(readBuffer, sendBuffer, my_lengths);

    /* Allocate memory */
    double *collideField = (double *) malloc((size_t)( Q*(my_lengths[0]+2)*(length[1]+2)*(length[2]+2)) * sizeof( double ));
    double *streamField  = (double *) malloc((size_t)( Q*(my_lengths[0]+2)*(length[1]+2)*(length[2]+2)) * sizeof( double ));
    int *flagField = (int *) malloc((size_t)( (my_lengths[0]+2)*(length[1]+2)*(length[2]+2) ) * sizeof( int ));
    double * massField = (double *) malloc((size_t)( (my_lengths[0]+2)*(length[1]+2)*(length[2]+2) ) * sizeof(double));
    /* fluid fraction field */
    double * fractionField = (double *) malloc((size_t)( (my_lengths[0]+2)*(length[1]+2)*(length[2]+2) ) * sizeof(double));

    double viscosity = C_S * C_S * (tau - 0.5);
    double Re = 1 / viscosity;

    if (collideField == NULL || streamField == NULL || flagField == NULL) {
        ERROR("Unable to allocate matrices.");
    } 
    
    initialiseFields(collideField, streamField, flagField, massField, fractionField, my_lengths, boundaries, r, argv, my_pos, Proc);

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

        exchange(LEFT, collideField, sendBuffer[0], readBuffer[0], my_lengths, my_pos, Proc, my_rank);
        exchange(RIGHT, collideField, sendBuffer[1], readBuffer[1], my_lengths, my_pos, Proc, my_rank);
        exchange(TOP, collideField, sendBuffer[2], readBuffer[2], my_lengths, my_pos, Proc, my_rank);
        exchange(BOTTOM, collideField, sendBuffer[3], readBuffer[3], my_lengths, my_pos, Proc, my_rank);
        exchange(BACK, collideField, sendBuffer[4], readBuffer[4], my_lengths, my_pos, Proc, my_rank);
        exchange(FRONT, collideField, sendBuffer[5], readBuffer[5], my_lengths, my_pos, Proc, my_rank);
 
        doStreaming(collideField, streamField, flagField, massField, fractionField, my_lengths);
        swap = collideField;
        collideField = streamField;
        streamField = swap;
        doCollision(collideField, flagField, massField, fractionField, &tau, my_lengths, extForces);
        updateFlagField(collideField, flagField, fractionField, my_lengths);
        treatBoundary(collideField, flagField, problem, &Re, &ro_ref, &ro_in, velocity, my_lengths);

        total_time += clock() - start_time; // Add elapsed ticks to total_time

        if (t % timestepsPerPlotting == 0) {
            run_checks(collideField, massField, flagField, my_lengths, t );
            writeVtkOutput(collideField, flagField, argv[1], t, my_lengths, my_origin, my_rank);
            printf("Time step %i finished, vtk file was created\n", t);
        }
    }

    /* Compute average mega-lattice-updates-per-second in order to judge performance */
    if (my_rank==0)
    {
        float elapsed_time = total_time/((float)CLOCKS_PER_SEC);
        float mlups = ((length[0]+2) * (length[1]+2) * (length[2]+2) * timesteps) / (elapsed_time * 1000000);
        printf("Elapsed time (excluding vtk writes) = %f\nAverage MLUPS = %f\n", elapsed_time, mlups);
    }

    free(collideField);
    free(streamField);
    free(flagField);
    free(massField);
    free(fractionField);
    
    for (int i = 0; i < 6; i++) {
        free(sendBuffer[i]);
        free(readBuffer[i]);
    }

    Programm_Stop("End");
 
    return 0;
}

#endif

