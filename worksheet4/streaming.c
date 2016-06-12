#include "streaming.h"
#include "helper.h"
#include "LBDefinitions.h"
#include "parallel.h"
#include <mpi.h>

void doStreaming(double * collideField, double * streamField, int * flagField, int * length){
    int x, y, z, i;
    int node[3], source_node[3];
    int n[3] = { length[0] + 2, length[1] + 2, length[2] + 2 };
    int myrank;

    printXZvelocities(collideField, 2, n);
    
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    /* Loop for inner cells */
     for (z = 1; z <= length[2]; z++) {
         node[2] = z;
         for (y = 1; y <= length[1]; y++) {
             node[1] = y;
             for (x = 1; x <= length[0]; x++) {
                 node[0] = x;
                 for (i = 0; i < Q; i++) {
                     source_node[0] = x - LATTICEVELOCITIES[i][0];
                     source_node[1] = y - LATTICEVELOCITIES[i][1];
                     source_node[2] = z - LATTICEVELOCITIES[i][2];
                     if (myrank == 0 && *getFlag(flagField, source_node, n) == PARALLEL) {
                         printf("(%i %i %i) %f -> (%i %i %i) %f\n",
                                node[0], node[1], node[2],
                                *getEl(streamField, node, i, n),
                                source_node[0], source_node[1], source_node[2],
                                *getEl(collideField, source_node, i, n)
                             );
                     }

                     *getEl(streamField, node, i, n) = *getEl(collideField, source_node, i, n);
                 }
             }
         }
     }
}

