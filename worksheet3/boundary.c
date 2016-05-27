#include "boundary.h"
#include "LBDefinitions.h"
#include "helper.h"
#include "computeCellValues.h"

void boundaryCell(double * collideField, int* flagField, const double * const wallVelocity, int xlength, int x, int y, int z){
    int i, coord_dest[3];
    double * cell_ptr;
    double dotProd;
    double density;
    int length_tot = xlength + 2;

    /* for each lattice */
    for (i = 0; i < Q; i++) {
        /* compute a cell where lattice is pointing*/
        coord_dest[0] = x + LATTICEVELOCITIES[i][0];
        coord_dest[1] = y + LATTICEVELOCITIES[i][1];
        coord_dest[2] = z + LATTICEVELOCITIES[i][2];

        /* if the cell is inner */
        if (coord_dest[0] > 0 && coord_dest[0] <= xlength &&
            coord_dest[1] > 0 && coord_dest[1] <= xlength &&
            coord_dest[2] > 0 && coord_dest[2] <= xlength) {

            /* get pointer to the i-th lattice of boundary cell */
            cell_ptr = getEl(collideField, x, y, z, i, length_tot);

            /* set i-th lattice to inverse lattice of the computed inner cell */
            *cell_ptr= *getEl(collideField, coord_dest[0], coord_dest[1], coord_dest[2], 18-i, length_tot);

            /* if boundary has type MOVING_WALL */
            if (flagField[length_tot*length_tot*z+length_tot*y+x] == MOVING_WALL) {
                dotProd = 0;

                /* compute inner product of wall velocity and i-th lattice velocity */
                for (int j = 0; j < 3; ++j)	{
                    dotProd += LATTICEVELOCITIES[i][j] * wallVelocity[j];
                }	
        
                computeDensity(getEl(collideField, coord_dest[0], coord_dest[1], coord_dest[2], 0, length_tot), &density);

                /* Set boundary i-th lattice with respect to the formula */
                *cell_ptr += 2.0 * LATTICEWEIGHTS[i] * dotProd * density / (C_S*C_S);
            }
        }
    }
}

void treatBoundary(double *collideField, int *flagField, const double * const wallVelocity, int xlength){
    int x, y, z;
    int n = xlength + 2;

    /* Left boundary for Y */
    z = 0;
    /* All cells on XY plane */
    for (y = 0; y < n; ++y){
        for (x = 0; x < n; ++x){           	        
            boundaryCell(collideField, flagField, wallVelocity, xlength, x, y, z);
        }
    }

    /* Inner cells for Z-axis */
    for (z = 1; z <= xlength; ++z) {
        /* Left boundary for Y */
        y = 0;
        /* All cells on X axis */
        for (x = 0; x < n; ++x) {
            boundaryCell(collideField, flagField, wallVelocity, xlength, x, y, z);
        }

        /* Inner cells for Y-axis */
        for (y = 1; y <= xlength; ++y){
            /* Left boundary for X */
            x = 0;
            boundaryCell(collideField, flagField, wallVelocity, xlength, x, y, z);
            /* Right boundary for X */
            x = xlength+1;
            boundaryCell(collideField, flagField, wallVelocity, xlength, x, y, z);
        }
        
        /* Right boundary for Y */
        y = xlength+1;
        /* All cells on X axis */
        for (x = 0; x < n; ++x){
            boundaryCell(collideField, flagField, wallVelocity, xlength, x, y, z);
        }
    }
    
    /* Right boundary for Z */
    z = xlength+1;
    /* All cells on XY plane */
    for (y = 0; y < n; ++y){
        for (x = 0; x < n; ++x){
            boundaryCell(collideField, flagField, wallVelocity, xlength, x, y, z);
        }
    }
}
