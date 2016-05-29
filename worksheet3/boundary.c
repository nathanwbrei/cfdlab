#include "boundary.h"
#include "LBDefinitions.h"
#include "helper.h"
#include "computeCellValues.h"

/**
 * Set NOSLIP condition
 */
void setNoSlip(double * collideField, int * flagField, int x, int y, int z, int * n) {
    int i, coord_dest[3];
    double * cell_ptr;

    /* for each lattice */
    for (i = 0; i < Q; i++) {
        /* compute a cell where lattice is pointing*/
        coord_dest[0] = x + LATTICEVELOCITIES[i][0];
        coord_dest[1] = y + LATTICEVELOCITIES[i][1];
        coord_dest[2] = z + LATTICEVELOCITIES[i][2];
        if (coord_dest[0]<n[0] && coord_dest[1]<n[1] && coord_dest[2]<n[2] && coord_dest[0]>=0 && coord_dest[1]>=0 && coord_dest[2]>=0){
            // TODO printf("%d %d %d \n",coord_dest[0],coord_dest[1],coord_dest[2] );
            /** TODO do we need to consider FLUID-OBSTACLE boundaries? */
            /* if pointed cell is FLUID */
            if (*getFlag(flagField, coord_dest[0], coord_dest[1], coord_dest[2], n) == FLUID) {
                /* get pointer to the i-th lattice of boundary cell */
                cell_ptr = getEl(collideField, x, y, z, i, n);
    
                /* NOSLIP */
                /* set i-th lattice to inverse lattice of the computed inner cell */
                *cell_ptr= *getEl(collideField, coord_dest[0], coord_dest[1], coord_dest[2], 18-i, n);
            }
        }
    }
}

/**
 * Set MOVING_WALL condition
 */
void setMovingWall(double * collideField, int * flagField,  const double * const wallVelocity, int x, int y, int z, int * n) {
    int i, coord_dest[3];
    double * cell_ptr;
    double dotProd;
    double density;
    
/* for each lattice */
    for (i = 0; i < Q; i++) {
        /* compute a cell where lattice is pointing*/
        coord_dest[0] = x + LATTICEVELOCITIES[i][0];
        coord_dest[1] = y + LATTICEVELOCITIES[i][1];
        coord_dest[2] = z + LATTICEVELOCITIES[i][2];

        if (coord_dest[0]<n[0] && coord_dest[1]<n[1] && coord_dest[2]<n[2] && coord_dest[0]>=0 && coord_dest[1]>=0 && coord_dest[2]>=0){
            /** TODO do we need to consider FLUID-OBSTACLE boundaries? */
            /* if pointed cell is FLUID */
            if (*getFlag(flagField, coord_dest[0], coord_dest[1], coord_dest[2], n) == FLUID) {
                /* get pointer to the i-th lattice of boundary cell */
                cell_ptr = getEl(collideField, x, y, z, i, n);

                /* NOSLIP */
                /* set i-th lattice to inverse lattice of the computed inner cell */
                *cell_ptr= *getEl(collideField, coord_dest[0], coord_dest[1], coord_dest[2], 18-i, n);

                dotProd = 0;

                /* compute inner product of wall velocity and i-th lattice velocity */
                for (int j = 0; j < 3; ++j)	{
                    dotProd += LATTICEVELOCITIES[i][j] * wallVelocity[j];
                }

                computeDensity(getEl(collideField, coord_dest[0], coord_dest[1], coord_dest[2], 0, n), &density);

                /* Set boundary i-th lattice with respect to the formula for MOVING_WALL */
                *cell_ptr += 2.0 * LATTICEWEIGHTS[i] * dotProd * density / (C_S*C_S);
            }
        }
    }
}

/**
 * Set OUTFLOW condition
 */
void setOutflow(double * collideField, int * flagField, const double * const ro_ref, int x, int y, int z, int * n) {
    int i, coord_dest[3];
    double * cell_ptr;
    double feq[Q];
    double velocity[D];
    double * fluidCell;

    /* for each lattice */
    for (i = 0; i < Q; i++) {
        /* compute a cell where lattice is pointing*/
        /** TODO Am I sure that we need to use those coordinates? */
        coord_dest[0] = x + LATTICEVELOCITIES[i][0];
        coord_dest[1] = y + LATTICEVELOCITIES[i][1];
        coord_dest[2] = z + LATTICEVELOCITIES[i][2];

        if (coord_dest[0]<n[0] && coord_dest[1]<n[1] && coord_dest[2]<n[2] && coord_dest[0]>=0 && coord_dest[1]>=0 && coord_dest[2]>=0){
            if (*getFlag(flagField, coord_dest[0], coord_dest[1], coord_dest[2], n) == FLUID) {
                /* get pointer to the fluid cell */
                fluidCell = getEl(collideField, coord_dest[0], coord_dest[1], coord_dest[2], 0, n);

                /* compute velocity of the fluid cell */
                computeVelocity(fluidCell, ro_ref, velocity);

                /* compute f-equilibrium of the fluid cell */
                computeFeq(ro_ref, velocity, feq);

                /* pointer to the i-th lattice of the boundary cell */
                cell_ptr = getEl(collideField, x, y, z, i, n);

                /* set boundary */
                *cell_ptr = feq[Q - i -1] + feq[i] - fluidCell[Q - 1 - i];
            }
        }
    }
}

/**
 * Set INFLOW condition
 */
void setInflow(double * collideField, int * flagField, const double * const ro_ref, const double * const inVelocity, int x, int y, int z, int * n) {
    int i, coord_dest[3];
    double * cell_ptr;
    double feq[Q];
    double * fluidCell;

    /* for each lattice */
    for (i = 0; i < Q; i++) {
        /* compute a cell where lattice is pointing*/
        /** TODO */
        coord_dest[0] = x + LATTICEVELOCITIES[i][0];
        coord_dest[1] = y + LATTICEVELOCITIES[i][1];
        coord_dest[2] = z + LATTICEVELOCITIES[i][2];

       if (coord_dest[0]<n[0] && coord_dest[1]<n[1] && coord_dest[2]<n[2] && coord_dest[0]>=0 && coord_dest[1]>=0 && coord_dest[2]>=0){
            if (*getFlag(flagField, coord_dest[0], coord_dest[1], coord_dest[2], n) == FLUID) {
                fluidCell = getEl(collideField, coord_dest[0], coord_dest[1], coord_dest[2], 0, n);
            
                computeFeq(ro_ref, inVelocity, feq);

                cell_ptr = getEl(collideField, x, y, z, i, n);

                *cell_ptr = feq[i];
            }
        }
    }
}

void boundaryCell(double * collideField,
                  int * flagField,
                  const double * const ro_ref,
                  const double * const inVelocity,
                  const double * const wallVelocity,
                  int x, int y, int z,
                  int * n) {
    /** TODO maybe it is more convenient to use length[0] for x and length[2] for z, and not vice versa */

    /* Type of boundary cell */
    int flag = *getFlag(flagField, x, y, z, n);

    if (flag == NOSLIP) {
        setNoSlip(collideField, flagField, x, y, z, n);
    } else if (flag == MOVING_WALL) {
        setMovingWall(collideField, flagField, wallVelocity, x, y, z, n);
    } else if (flag == INFLOW) {
        setInflow(collideField, flagField, ro_ref, inVelocity, x, y, z, n);
    } else if (flag == OUTFLOW) {
        setOutflow(collideField, flagField, ro_ref, x, y, z, n);
    }
}

void treatBoundary(double *collideField, int *flagField, const double * const ro_ref, const double * const inVelocity, const double * const wallVelocity, int * length){
    int x, y, z, flag;
    int n[3] = { length[0] + 2, length[1] + 2, length[2] + 2 };

    for (z = 0; z < n[0]; z++) {
        for (y = 0; y < n[1]; y++) {
            for (x = 0; x < n[2]; x++) {
                flag = *getFlag(flagField, x, y, z, n);
                if (flag != FLUID && flag != OBSTACLE) {
                    //printf("Debug: started boundaryCell %d %d %d\n", z,y,x);
                    boundaryCell(collideField, flagField, ro_ref, inVelocity, wallVelocity, x, y, z, n);
                }
            }
        }
    }
    printf("Debug: ended treatBoundary\n");

//    
///*
//  /** TODO maybe it is more convenient to use length[0] for x and length[2] for z, and not vice versa */
//
//    /* Left boundary for Y */
//    z = 0;
//    /* All cells on XY plane */
//    for (y = 0; y < n[1]; ++y){
//        for (x = 0; x < n[2]; ++x) {
//            boundaryCell(collideField, flagField, wallVelocity, length, x, y, z);
//        }
//    }
//
//    /* Inner cells for Z-axis */
//    for (z = 1; z <= length[0]; ++z) {
//        /* Left boundary for Y */
//        y = 0;
//        /* All cells on X axis */
//        for (x = 0; x < n[2]; ++x) {
//            boundaryCell(collideField, flagField, wallVelocity, length, x, y, z);
//        }
//
//        /* Inner cells for Y-axis */
//        for (y = 1; y <= length[1]; ++y){
//            /* Left boundary for X */
//            x = 0;
//            boundaryCell(collideField, flagField, wallVelocity, length, x, y, z);
//            
//            /* Right boundary for X */
//            x = length+1;
//            boundaryCell(collideField, flagField, wallVelocity, length, x, y, z);
//        }
//        
//        /* Right boundary for Y */
//        y = length+1;
//        /* All cells on X axis */
//        for (x = 0; x < n[2]; ++x){
//            boundaryCell(collideField, flagField, wallVelocity, length, x, y, z);
//        }
//    }
//    
//    /* Right boundary for Z */
//    z = length+1;
//    /* All cells on XY plane */
//    for (y = 0; y < n[1]; ++y){
//        for (x = 0; x < n[2]; ++x){
//            boundaryCell(collideField, flagField, wallVelocity, length, x, y, z);
//        }
//    }
    
}
