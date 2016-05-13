#include "streaming.h"
#include "helper.h"
#include "LBDefinitions.h"

void doStreaming(double * collideField, double * streamField, int * flagField, int xlength){
    int x, y, z, i;
    int n = xlength + 2;

    /* Loop for inner cells */
    for (x = 1; x < n; x++) {
        for (y = 1; y < n; y++) {
            for (z = 1; z < n; z++) {
                for (i = 1; i < 19; i++) {
                    *getEl(collideField, x, y, z, i, xlength) = *getEl(streamField, x - LATTICEVELOCITIES[i][0], y - LATTICEVELOCITIES[i][1], z - LATTICEVELOCITIES[i][2], i, xlength);
                }
            }
        }
    }
}

