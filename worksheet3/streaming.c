#include "streaming.h"
#include "helper.h"
#include "LBDefinitions.h"

void doStreaming(double * collideField, double * streamField, int * flagField, int * length){
    int x, y, z, i;
    int n[3] = { length[0] + 2, length[1] + 2, length[2] + 2 };

    /* Loop for inner cells */
    for (z = 1; z <= length[0]; z++) {
        for (y = 1; y <= length[1]; y++) {
            for (x = 1; x <= length[2]; x++) {
                if (*getFlag(flagField, x, y, z, n) == FLUID) {
                    for (i = 0; i < Q; i++) {
                        *getEl(streamField, x, y, z, i, n) = *getEl(collideField, x - LATTICEVELOCITIES[i][0], y - LATTICEVELOCITIES[i][1], z - LATTICEVELOCITIES[i][2], i, n);
                    }
                }
            }
        }
    }
}

