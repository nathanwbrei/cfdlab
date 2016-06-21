#include "LBDefinitions.h"
#include "helper.h"
#include "flag.h"

void updateFlagField(double * collideField, int * flagField, double * fractionField, int * length) {
    int x, y, z, flag, nFilled = 0, nEmptied = 0, k;
    int node[3];
    double fraction;
    int n[3] = { length[0] + 2, length[1] + 2, length[2] + 2 };
    
    int filledCells[n[0] * n[1] * n[2]][3];
    int emptiedCells[n[0] * n[1] * n[2]][3];

    /*
      Updating flags for INTERFACE cells:
         if fraction > 1 set flag to FLUID;
         if fraction < 0 set flag to GAS.
      Saving all emptied and filled cells to emptiedCells and filledCells arrays.
    */
    /* TODO fill filledCells and emptiedCells on collide step */
    for (z = 1; z <= length[2]; z++) {
        node[2] = z;
        for (y = 1; y <= length[1]; y++) {
            node[1] = y;
            for (x = 1; x <= length[0]; x++) {
                node[0] = x;
                flag = *getFlag(flagField, node, n);

                /* We are interested only in INTERFACE cells now */
                if (flag == INTERFACE) {
                    fraction = *getFraction(fractionField, node, n);

                    if (fraction > 1) {
                        filledCells[nFilled][0] = node[0];
                        filledCells[nFilled][1] = node[1];
                        filledCells[nFilled][2] = node[2];
                        nFilled++;

                        *getFlag(flagField, node, n) = FLUID;
                    } else if (fraction < 0) {
                        emptiedCells[nEmptied][0] = node[0];
                        emptiedCells[nEmptied][1] = node[1];
                        emptiedCells[nEmptied][2] = node[2];
                        nEmptied++;

                        *getFlag(flagField, node, n) = GAS;
                    }
                }
            }
        }
    }

    /*
      On this step we need to update neighbors of filled and emptied cells
      in order to have closed interface layer
     */

    /*
      First - neighbors of the filled cells.
      We need to convert all empty neighbors of filled cells to INTERFACE.
      Distribution of these cells will be initialized with f_eq(rho_avg, v_avg),
      where v_avg and rho_avg is average characteristics of surrounding FLUID and INTERFACE cells.
      All such neighbors must be deleted from emptiedCells array.
     */
    for (k = 0; k < nFilled; k++) {
        /* TODO */
    }

    /*
      Second - neighbors of the emptied cells.
     */
    updateEmptiedNeighbors(collideField, flagField, emptiedCells, nEmptied, n);
}

void updateEmptiedNeighbors(double * collideField, int * flagField, int emptiedCells[][3], int nEmptied, int * n) {
    int i, k, neighbor_node[3], * flag;

    for (k = 0; k < nEmptied; k++) {
        for (i = 0; i < Q; i++) {
            neighbor_node[0] = emptiedCells[k][0] + LATTICEVELOCITIES[i][0];
            neighbor_node[1] = emptiedCells[k][1] + LATTICEVELOCITIES[i][1];
            neighbor_node[2] = emptiedCells[k][2] + LATTICEVELOCITIES[i][2];
            //         printf("%i %i %i\n", neighbor_node[0], neighbor_node[1], neighbor_node[2]);
            if (neighbor_node[0] == 0 || neighbor_node[0] == n[0] ||
                neighbor_node[1] == 0 || neighbor_node[1] == n[1] ||
                neighbor_node[2] == 0 || neighbor_node[2] == n[2]) {
                continue;
            }

            flag = getFlag(flagField, neighbor_node, n);

            if (*flag == FLUID) {
                *flag = INTERFACE;
            }
        }
    }

}
