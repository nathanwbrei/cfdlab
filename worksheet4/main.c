#ifndef _MAIN_C_
#define _MAIN_C_

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
    int my_rank, number_of_ranks, Proc[3], my_pos[3], my_lengths[3];
    double tau, velocityWall[3];    
    int t;
    double *swap=NULL;
    clock_t start_time, total_time = 0;
    double *sendBuffer[6];// [0:left,1:right,2:top,3:bottom,4:front,5:back]
    double *readBuffer[6];

    
    //initializeMPI(&my_rank,&number_of_ranks,argc,argv);
    /* initializeMPI */
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &number_of_ranks );
    MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );

    readParameters(&xlength, &tau, velocityWall, &timesteps, &timestepsPerPlotting, argc, argv, Proc, my_rank);

    /* If the number of process in the imput file and the arguments does not coicide, launch an error */
    if (number_of_ranks != (Proc[0]*Proc[1]*Proc[2])){
        printf("Proc %d : Nuber of processes does not match %d vs %d \n", my_rank, number_of_ranks,(Proc[0]*Proc[1]*Proc[2]));
    }

    get_rank_pos( my_pos, my_rank, Proc);
    get_my_lengths( my_pos, xlength, my_lengths, Proc);    

    MPI_Barrier( MPI_COMM_WORLD ); /* Barrier to get order in the output, just for thebug */
    printf("Debug: Process %d: Position x, y z %d %d %d ; lenghts %d %d %d \n", my_rank, my_pos[0], my_pos[1], my_pos[2], my_lengths[0], my_lengths[1], my_lengths[2]);

    double  *collideField   = (double *)    malloc((size_t)( Q*(my_lengths[0]+2)*(my_lengths[1]+2)*(my_lengths[2]+2) ) * sizeof( double ));
    double  *streamField    = (double *)    malloc((size_t)( Q*(my_lengths[0]+2)*(my_lengths[1]+2)*(my_lengths[2]+2) ) * sizeof( double ));
    int     *flagField      = (int *)       malloc((size_t)(   (my_lengths[0]+2)*(my_lengths[1]+2)*(my_lengths[2]+2) ) * sizeof( int ));

    if (collideField == NULL || streamField == NULL || flagField == NULL) {
        ERROR("Unable to allocate matrices.");
    } 


/*
    initialiseFields(collideField, streamField, flagField, my_lengths);

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

*/
    free(collideField);
    free(streamField);
    free(flagField);

    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize();
    return 0;
}

#endif

