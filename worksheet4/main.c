#ifndef _MAIN_C_
#define _MAIN_C_

#include "parallel.h"
#include <time.h>
#include <mpi.h>

#include "collision.h"
#include "streaming.h"
#include "initLB.h"
#include "visualLB.h"
#include "boundary.h"
#include "LBDefinitions.h"


int main(int argc, char *argv[]){
    int xlength, timesteps, timestepsPerPlotting;
    int my_rank, number_of_ranks, Proc[3], my_pos[3], my_lengths[3], my_origin[3];
    double tau, velocityWall[3];    
    int t;
    clock_t start_time, total_time = 0;
    double * sendBuffer[6]; /* [0:left,1:right,2:top,3:bottom,4:front,5:back] */
    double * readBuffer[6];
    double * s;

    initializeMPI(&my_rank, &number_of_ranks, argc, argv);

    readParameters(&xlength, &tau, velocityWall, &timesteps, &timestepsPerPlotting, argc, argv, Proc, my_rank);

    /* If the number of process in the imput file and the arguments does not coicide, launch an error */
    if (number_of_ranks != (Proc[0]*Proc[1]*Proc[2])){
        printf("Proc %d : Number of processes does not match %d vs %d \n", my_rank, number_of_ranks,(Proc[0]*Proc[1]*Proc[2]));
    }

    get_rank_pos(my_pos, my_rank, Proc);
    get_my_lengths(my_pos, xlength, Proc, my_lengths, my_origin);

    initBuffers(readBuffer, sendBuffer, my_lengths);


    double *collideField = (double *) malloc((size_t)(Q*(my_lengths[0]+2)*(my_lengths[1]+2)*(my_lengths[2]+2)) * sizeof(double));
    double *streamField = (double *) malloc((size_t)(Q*(my_lengths[0]+2)*(my_lengths[1]+2)*(my_lengths[2]+2)) * sizeof(double));
    int *flagField = (int *) malloc((size_t)((my_lengths[0]+2)*(my_lengths[1]+2)*(my_lengths[2]+2) ) * sizeof(int));

    if (collideField == NULL || streamField == NULL || flagField == NULL) {
        ERROR("Unable to allocate matrices.");
    }

    initialiseFields(collideField, streamField, flagField, my_lengths, my_pos, Proc);

    treatBoundary(collideField, flagField, velocityWall, my_lengths);

    for (t = 0; t < timesteps; t++) {
        start_time = clock();  // Start the timer for the lattice updates

        exchange(LEFT, collideField, sendBuffer[0], readBuffer[0], my_lengths, my_pos, Proc, my_rank);
        exchange(RIGHT, collideField, sendBuffer[1], readBuffer[1], my_lengths, my_pos, Proc, my_rank);
        exchange(TOP, collideField, sendBuffer[2], readBuffer[2], my_lengths, my_pos, Proc, my_rank);
        exchange(BOTTOM, collideField, sendBuffer[3], readBuffer[3], my_lengths, my_pos, Proc, my_rank);
        exchange(BACK, collideField, sendBuffer[4], readBuffer[4], my_lengths, my_pos, Proc, my_rank);
        exchange(FRONT, collideField, sendBuffer[5], readBuffer[5], my_lengths, my_pos, Proc, my_rank);
 
        // TODO: Is this ordering correct? Remember that we send the RIGHT face to the right

        doStreaming(collideField, streamField, flagField, my_lengths);

        s = collideField;
        collideField = streamField;
        streamField = s;

        doCollision(collideField,flagField,&tau,my_lengths);
        treatBoundary(collideField, flagField, velocityWall, my_lengths);

        total_time += clock() - start_time; /* Add elapsed ticks to total_time */

        if (t % timestepsPerPlotting == 0) {
            printf("Process %i finished step %i\n", my_rank, t);
            writeVtkOutput(collideField, flagField, argv[1], t, my_lengths, my_origin, my_pos, my_rank);
        }
    }

    /* Compute average mega-lattice-updates-per-second in order to judge performance */
    if (my_rank==0)
    {
        float elapsed_time = total_time/((float)CLOCKS_PER_SEC);
        float mlups = ((xlength+2) * (xlength+2) * (xlength+2) * timesteps) / (elapsed_time * 1000000);
        printf("Elapsed time (excluding vtk writes) = %f\nAverage MLUPS = %f\n", elapsed_time, mlups);
    }

    free(collideField);
    free(streamField);
    free(flagField);

    /* TODO free buffers */

    Programm_Stop("End");
    
    return 0;
}

#endif

