#include "streaming.h"
#include "helper.h"
#include "LBDefinitions.h"
#include "parallel.h"

void doStreaming(double * collideField, double * streamField, int * flagField, int * length){
    int x, y, z, i;
    int node[3], source_node[3];
    int n[3] = { length[0] + 2, length[1] + 2, length[2] + 2 };

    //printXZvelocities(collideField, 2, n);

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
                     if (*getFlag(flagField, source_node, n) == PARALLEL) {
                         // printf("streaming from parallel boundary\n");
                     }

                     *getEl(streamField, node, i, n) = *getEl(collideField, source_node, i, n);
                 }
             }
         }
     }
}

