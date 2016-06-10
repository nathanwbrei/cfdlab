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

// TODO: Move these to their own file

typedef enum {LEFT, RIGHT, TOP, BOTTOM, FRONT, BACK} face_t;

double * getBufferEl(double * buffer, int m, int n, int i, int * length) {
    return buffer + n * length[0] * Q + m * Q + i;
}

void extract_z(double * field,  double ** sendBuffer, int * node, int * length, int * n) {
    int x, y, i, bufferLength[2];

    bufferLength[0] = length[0];
    bufferLength[1] = length[1];

    for (y = 1; y <= length[1]; y++) {
        node[1] = y;
        for (x = 0; x < n[0]; x++) {
            node[0] = x;
            for (i = 1; i < Q; i++) {
                *getBufferEl(sendBuffer[0], x, y, i, bufferLength) = *getEl(field, node, i, n);
            }
        }
    }
}

void extract_y(double * field,  double ** sendBuffer, int * node, int * length, int * n) {
    int x, z, i, bufferLength[2];

    bufferLength[0] = length[0];
    bufferLength[1] = length[2];
 
    for (z = 0; z < n[2]; z++) {
        node[2] = z;
        for (x = 0; x < n[0]; x++) {
            node[0] = x;
            for (i = 1; i < Q; i++) {
                *getBufferEl(sendBuffer[0], x, z, i, bufferLength) = *getEl(field, node, i, n);
            }
        }
    }
}

void extract_x(double * field,  double ** sendBuffer, int * node, int * length, int * n) {
    int y, z, i, bufferLength[2];

    bufferLength[0] = length[1];
    bufferLength[1] = length[2];
 
    for (z = 1; z <= length[2]; z++) {
        node[2] = z;
        for (y = 1; y <= length[1]; y++) {
            node[1] = y;
            for (i = 1; i < Q; i++) {
                *getBufferEl(sendBuffer[0], y, z, i, bufferLength) = *getEl(field, node, i, n);
            }
        }
    }
}

void extract(face_t face, double * field, double ** sendBuffer, int * length, int * my_pos, int * Proc) {
    int n[3] = { length[0] + 2, length[1] + 2, length[2] + 2 };
    int node[3];

    if (face == LEFT) {
        if (my_pos[0] != 0) {
            node[0] = 0;
            extract_x(field, sendBuffer, node, length, n);
       }
    } else if (face == RIGHT) {
        if (my_pos[0] != Proc[0] - 1) {
            node[0] = length[0] + 1;
            extract_x(field, sendBuffer, node, length, n);
        }
    } else if (face == TOP) {
        if (my_pos[2] != Proc[2] - 1) {
            node[2] = length[2] + 1;
            extract_z(field, sendBuffer, node, length, n);
        }
    } else if (face == BOTTOM) {
        if (my_pos[2] != 0) {
            node[2] = 0;
            extract_z(field, sendBuffer, node, length, n);
        }
    } else if (face == FRONT) {
        if (my_pos[1] != 0) {
            node[1] = 0;
            extract_y(field, sendBuffer, node, length, n);
        }
    } else if (face == BACK) {
        if (my_pos[1] != Proc[1] - 1) {
            node[1] = length[1] + 1;
            extract_y(field, sendBuffer, node, length, n);
        }
    }

    // Copy the cells of `field` along `face` into corresponding slot of `sendBuffer`
}

void swap(face_t face, double * sendBuffer, double * receiveBuffer) {
    // MPI send to neighboring process, block until receives 
    // TODO: Needs to know who its neighbors are / whether they exist
}

void inject_z(double * field,  double ** readBuffer, int * node, int * length, int * n) {
    int x, y, i, bufferLength[2];

    bufferLength[0] = length[0];
    bufferLength[1] = length[1];
 
    for (y = 1; y <= length[1]; y++) {
        node[1] = y;
        for (x = 0; x < n[0]; x++) {
            node[0] = x;
            for (i = 1; i < Q; i++) {
                *getEl(field, node, i, n) = *getBufferEl(readBuffer[0], x, y, i, bufferLength);
            }
        }
    }
}

void inject_y(double * field,  double ** readBuffer, int * node, int * length, int * n) {
    int x, z, i, bufferLength[2];

    bufferLength[0] = length[0];
    bufferLength[1] = length[2];
 
    for (z = 0; z < n[2]; z++) {
        node[2] = z;
        for (x = 0; x < n[0]; x++) {
            node[0] = x;
            for (i = 1; i < Q; i++) {
                *getEl(field, node, i, n) = *getBufferEl(readBuffer[0], x, z, i, bufferLength);
            }
        }
    }
}

void inject_x(double * field,  double ** readBuffer, int * node, int * length, int * n) {
    int y, z, i, bufferLength[2];

    bufferLength[0] = length[1];
    bufferLength[1] = length[2];
 
    for (z = 1; z <= length[2]; z++) {
        node[2] = z;
        for (y = 1; y <= length[1]; y++) {
            node[1] = y;
            for (i = 1; i < Q; i++) {
                *getEl(field, node, i, n) = *getBufferEl(readBuffer[0], y, z, i, bufferLength);
            }
        }
    }
}

void inject(face_t face, double * readBuffer, double * field, int * length, int * my_pos, int * Proc) {
    // Copy the contents of the received buffer into the field along the corresponding face
    int n[3] = { length[0] + 2, length[1] + 2, length[2] + 2 };
    int node[3];

    if (face == LEFT) {
        if (my_pos[0] != 0) {
            node[0] = 0;
            inject_x(field, readBuffer, node, length, n);
       }
    } else if (face == RIGHT) {
        if (my_pos[0] != Proc[0] - 1) {
            node[0] = length[0] + 1;
            inject_x(field, readBuffer, node, length, n);
        }
    } else if (face == TOP) {
        if (my_pos[2] != Proc[2] - 1) {
            node[2] = length[2] + 1;
            inject_z(field, readBuffer, node, length, n);
        }
    } else if (face == BOTTOM) {
        if (my_pos[2] != 0) {
            node[2] = 0;
            inject_z(field, readBuffer, node, length, n);
        }
    } else if (face == FRONT) {
        if (my_pos[1] != 0) {
            node[1] = 0;
            inject_y(field, readBuffer, node, length, n);
        }
    } else if (face == BACK) {
        if (my_pos[1] != Proc[1] - 1) {
            node[1] = length[1] + 1;
            inject_y(field, readBuffer, node, length, n);
        }
    }
}

void exchange(face_t face,
              double * field,
              double ** sendBuffer,
              double ** readBuffer,
              int * length,
              int * my_pos,
              int * Proc) {
    extract(face, field, sendBuffer, length, my_pos, Proc);  
    swap(face, sendBuffer, readBuffer);
    inject(face, field, readBuffer, length, my_pos, Proc);
}

int main(int argc, char *argv[]){
    int xlength, timesteps, timestepsPerPlotting;
    int my_rank, number_of_ranks, Proc[3], my_pos[3], my_lengths[3];
    double tau, velocityWall[3];    
    int t;
    double *swap=NULL;
    clock_t start_time, total_time = 0;
    double * sendBuffer[6]; /* [0:left,1:right,2:top,3:bottom,4:front,5:back] */
    double * readBuffer[6];

    initializeMPI(&my_rank, &number_of_ranks, argc, argv);
    /* initializeMPI */
    //MPI_Init( &argc, &argv );
    //MPI_Comm_size( MPI_COMM_WORLD, &number_of_ranks );
    //MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );

    readParameters(&xlength, &tau, velocityWall, &timesteps, &timestepsPerPlotting, argc, argv, Proc, my_rank);

    /* If the number of process in the imput file and the arguments does not coicide, launch an error */
    if (number_of_ranks != (Proc[0]*Proc[1]*Proc[2])){
        printf("Proc %d : Number of processes does not match %d vs %d \n", my_rank, number_of_ranks,(Proc[0]*Proc[1]*Proc[2]));
    }

    get_rank_pos(my_pos, my_rank, Proc);
    get_my_lengths(my_pos, xlength, my_lengths, Proc);
    initBuffers(readBuffer, sendBuffer, my_lengths);

    MPI_Barrier(MPI_COMM_WORLD); /* Barrier to get order in the output, just for thebug */
    printf("Debug: Process %2d: Position x, y, z %d %d %d ; lenghts %d %d %d \n", my_rank, my_pos[0], my_pos[1], my_pos[2], my_lengths[0], my_lengths[1], my_lengths[2]);

    double *collideField = (double *) malloc((size_t)(Q*(my_lengths[0]+2)*(my_lengths[1]+2)*(my_lengths[2]+2)) * sizeof(double));
    double *streamField = (double *) malloc((size_t)(Q*(my_lengths[0]+2)*(my_lengths[1]+2)*(my_lengths[2]+2)) * sizeof(double));
    int *flagField = (int *) malloc((size_t)((my_lengths[0]+2)*(my_lengths[1]+2)*(my_lengths[2]+2) ) * sizeof(int));

    if (collideField == NULL || streamField == NULL || flagField == NULL) {
        ERROR("Unable to allocate matrices.");
    }

    initialiseFields(collideField, streamField, flagField, my_lengths, my_pos, Proc);

    treatBoundary(collideField, flagField, velocityWall, my_lengths);
    
    writeVtkOutput(collideField, flagField, argv[1], t, my_lengths, my_pos, my_rank);
    /*
    for (t = 0; t < timesteps; t++) {

        start_time = clock();  // Start the timer for the lattice updates

        exchange(RIGHT, collideField, sendBuffer, receiveBuffer);
        exchange(LEFT, collideField, sendBuffer, receiveBuffer);
        exchange(BACK, collideField, sendBuffer, receiveBuffer);
        exchange(FRONT, collideField, sendBuffer, receiveBuffer);
        exchange(BOTTOM, collideField, sendBuffer, receiveBuffer);
        exchange(TOP, collideField, sendBuffer, receiveBuffer);
        // TODO: Is this ordering correct? Remember that we send the RIGHT face to the right

        doStreaming(collideField, streamField, flagField, xlength);

        double * swap = collideField;
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

    Programm_Stop("End");
    
    return 0;
}

#endif

