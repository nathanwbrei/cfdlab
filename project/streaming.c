#include "streaming.h"
#include "helper.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"
#include <omp.h>

void doStremingCell(double * collideField, double * streamField, int * flagField, double * massField, double * fractionField, int * node, double * el, int * n, int isInterface) {
    int i, flag;
    int source_node[3];
    double fi_nb, se;

    for (i = 0; i < Q; i++) {
        /* neighboring cell from which particles are obtained */
        source_node[0] = node[0] - LATTICEVELOCITIES[i][0];
        source_node[1] = node[1] - LATTICEVELOCITIES[i][1];
        source_node[2] = node[2] - LATTICEVELOCITIES[i][2];

        /* Amount of particles that goes to this cell */
        fi_nb = *getEl(collideField, source_node, i, n);
        *(el + i) = fi_nb;

        if (isInterface) {
        /*
          Update mass field according formula (4.3) (PhD thesis):
          dm_i(x, t + dt) = se * (E(x + dt*e_i, t) + E(x, t)) / 2
          se = f_i_inv(x + dt*e_i, t) - f_i(x,t)

          m(x, t + dt) = m(x, t) + sum(dm_i(x, t + dt));
          i in formula corresponds to Q - 1 - i in code and vice versa
        */

            flag = *getFlag(flagField, source_node, n); 
            if (flag == FLUID || flag == INTERFACE) {
                se = fi_nb - *getEl(collideField, source_node, Q - 1 - i, n);
                *getMass(massField, node, n) += se * (*getFraction(fractionField, node, n) + *getFraction(fractionField, source_node, n)) / 2;
            }
        }
    }

}

void doStreaming(double * collideField, double * streamField, int * flagField, double * massField, double * fractionField, int * length){
    int x, y, z, *flag, isFluid, isInterface;
    int node[3];
    double * el;

    int n[3] = { length[0] + 2, length[1] + 2, length[2] + 2 };

    /* Loop for inner cells */
//#pragma omp parallel for schedule(dynamic) private(y, x, node, isFluid, flag, isInterface, el) num_threads(6)
    for (z = 1; z <= length[2]; z++) {
        node[2] = z;
         
        for (y = 1; y <= length[1]; y++) {
            node[1] = y;
            for (x = 1; x <= length[0]; x++) {
                el = getEl(streamField, node, 0, n);
                flag = getFlag(flagField, node, n);
                node[0] = x;
                isFluid = *flag == FLUID;
                isInterface = *flag == INTERFACE;

                if (isFluid || isInterface) {
                    doStremingCell(collideField, streamField, flagField, massField, fractionField, node, el, n, isInterface);
                }


            }
        }
    }
}

