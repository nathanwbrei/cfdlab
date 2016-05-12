#include "boundary.h"
#include "LBDefinitions.h"
#include "helper.h"

double density(double * cell){
	double density = 0;
	for (int i = 0; i < 19; ++i){
		density += cell[i];
	}
	return density;
}

void boundaryCell(double * collideField, int* flagField, const double * const wallVelocity, int length_tot, int x, int y, int z){
	int i, coord_dest[3];
	double * cell_ptr;
	double DotProd;

	for (i=0;i<19;i++){
		coord_dest[0] = x + LATTICEVELOCITIES[i][0];
		coord_dest[1] = y + LATTICEVELOCITIES[i][1];
		coord_dest[2] = z + LATTICEVELOCITIES[i][2];
		if (coord_dest[0] >0 && coord_dest[0] <length_tot && coord_dest[1] >0 && coord_dest[1] <length_tot && coord_dest[2] >0 && coord_dest[2] <length_tot ){
			cell_ptr = getEl(collideField, x, y, z, i, length_tot);
			* cell_ptr= *getEl(collideField, coord_dest[0], coord_dest[1], coord_dest[2], 18-i, length_tot);
			if (flagField[length_tot*length_tot*z+length_tot*y+x] ==2 ){			
				for (int j = 0; j < 3; ++j)	{
						DotProd += (double) LATTICEVELOCITIES[i][j] * wallVelocity[j];
					}	
				*cell_ptr += 2.0 * LATTICEWEIGHTS[i] * DotProd/ (C_S*C_S);
				*cell_ptr += density(getEl(collideField, coord_dest[0], coord_dest[1], coord_dest[2], 0, length_tot));
			}
		}
	}
}

void treatBoundary(double *collideField, int *flagField, const double * const wallVelocity, int xlength){
	int x, y, z, length_tot=xlength +2;

	z = 0;
	for (y = 0; y <= xlength+1; ++y){
	    for (x = 0; x <= xlength+1; ++x){           	        
	        boundaryCell(collideField, flagField, wallVelocity, length_tot, x, y, z);
	    }
	}
	for (z = 1; z <= xlength; ++z){
	    y = 0;
	    for (x = 0; x <= xlength +1; ++x){
	        boundaryCell(collideField, flagField, wallVelocity, length_tot, x, y, z);
	    }       
	    for (y = 1; y <= xlength; ++y){
	        x = 0;
	        boundaryCell(collideField, flagField, wallVelocity, length_tot, x, y, z);
	        x = xlength+1;
	        boundaryCell(collideField, flagField, wallVelocity, length_tot, x, y, z);
	    }
	    y = xlength+1;
	    for (x = 0; x <= xlength+1; ++x){
	        boundaryCell(collideField, flagField, wallVelocity, length_tot, x, y, z);
	    }
	}
	z = xlength+1;
	for (y = 0; y <= xlength+1; ++y){
	    for (x = 0; x <= xlength+1; ++x){
	        boundaryCell(collideField, flagField, wallVelocity, length_tot, x, y, z);
	    }
	}
}