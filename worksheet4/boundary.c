#include "boundary.h"
#include "LBDefinitions.h"
#include "helper.h"
#include "computeCellValues.h"

void boundaryCell(double * collideField, int* flagField, const double * const wallVelocity, int * n, int * node){
    int i, coord_dest[3];
    double * cell_ptr;
    double dotProd;
    double density;

    for (i = 0; i < Q; i++) {
        coord_dest[0] = node[0] + LATTICEVELOCITIES[i][0];
        coord_dest[1] = node[1] + LATTICEVELOCITIES[i][1];
        coord_dest[2] = node[2] + LATTICEVELOCITIES[i][2];
    
        if (coord_dest[0] > 0 && coord_dest[0] < n[0] &&
            coord_dest[1] > 0 && coord_dest[1] < n[1] &&
            coord_dest[2] > 0 && coord_dest[2] < n[2]) {
      
            cell_ptr = getEl(collideField, node , i, n);
            *cell_ptr= *getEl(collideField, coord_dest, 18-i, n);

            if (*getFlag(flagField, node, n) == MOVING_WALL) {
                dotProd = 0;
                
                for (int j = 0; j < 3; ++j)	{
                    dotProd += LATTICEVELOCITIES[i][j] * wallVelocity[j];
                }	
        
                computeDensity(getEl(collideField, coord_dest, 0, n), &density);
                *cell_ptr += 2.0 * LATTICEWEIGHTS[i] * dotProd * density / C_S_2;
            }
        }
    }
}

void treatBoundary(double *collideField, int *flagField, const double * const wallVelocity, int * length){
    int x, y, z;
    int node[3];
    int n[3] = { length[0] + 2, length[1] + 2, length[2] + 2 };

    /* z = 0 */
    node[2] = 0;
    for (y = 0; y < n[1]; ++y){
        node[1] = y;
        for (x = 0; x < n[0]; ++x){           	        
            node[0] = x;
            boundaryCell(collideField, flagField, wallVelocity, n, node);
        }
    }

    for (z = 1; z <= length[2]; ++z){
        /* y = 0 */
        node[1] = 0;
        for (x = 0; x < n[0]; ++x){
            node[0] = x;
            boundaryCell(collideField, flagField, wallVelocity, n, node);
        }

        for (y = 1; y <= length[1]; ++y){
            node[1] = y;
            // x = 0;
            node[0] = 0;
            boundaryCell(collideField, flagField, wallVelocity, n, node);
            // x = length[0] + 1;
            node[0] = length[0] + 1;
            boundaryCell(collideField, flagField, wallVelocity, n, node);
        }

        y = length[1] + 1;
        node[1] = y;
        for (x = 0; x < n[0]; ++x){
            node[0] = x;
            boundaryCell(collideField, flagField, wallVelocity, n, node);
        }
    }

    /* z = length[2] + 1; */
    node[2] = length[2] + 1;
    for (y = 0; y < n[1]; ++y){
        node[1] = y;
        for (x = 0; x < n[0]; ++x){
            node[0] = x;
            boundaryCell(collideField, flagField, wallVelocity, n, node);
        }
    }
}
