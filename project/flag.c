#include "LBDefinitions.h"
#include "helper.h"

void updateFlagField(double * collideField, int * flagField, double * fractionField, int * length) {
    int x, y, z, flag, nFilled = 0, nEmptied = 0, k, i;
    int node[3];
    double fraction, * cell;
    int n[3] = { length[0] + 2, length[1] + 2, length[2] + 2 };
    
    double * filledCells[n[0] * n[1] * n[2]];
    double * emptiedCells[n[0] * n[1] * n[2]];

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
                    cell = getEl(collideField, node, 0, n);
                    fraction = *getFraction(fractionField, node, n);

                    if (fraction > 1) {
                        filledCells[nFilled++] = cell;
                        *getFlag(flagField, node, n) = FLUID;
                    } else if (fraction < 0) {
                        emptiedCells[nEmptied++] = cell;
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
    for (k = 0; k < nEmptied; k++) {
        /* TODO */
    }
}
