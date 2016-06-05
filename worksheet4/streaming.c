#include "streaming.h"
#include "helper.h"
#include "LBDefinitions.h"

void doStreaming(double * collideField, double * streamField, int * flagField, int xlength){
    int x, y, z, i;
    int n = xlength + 2;

    /* Loop for inner cells */
    for (z = 1; z <= xlength; z++) {
        for (y = 1; y <= xlength; y++) {
            for (x = 1; x <= xlength; x++) {
                for (i = 0; i < Q; i++) {
                    *getEl(streamField, x, y, z, i, n) = *getEl(collideField, x - LATTICEVELOCITIES[i][0], y - LATTICEVELOCITIES[i][1], z - LATTICEVELOCITIES[i][2], i, n);
                }
            }
        }
    }
}

