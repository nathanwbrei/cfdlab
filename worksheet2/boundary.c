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

  for (i = 0; i < 19; i++){
    coord_dest[0] = x + LATTICEVELOCITIES[i][0];
    coord_dest[1] = y + LATTICEVELOCITIES[i][1];
    coord_dest[2] = z + LATTICEVELOCITIES[i][2];
    
    if (coord_dest[0] > 0 && coord_dest[0] <= xlength &&
        coord_dest[1] > 0 && coord_dest[1] <= xlength &&
        coord_dest[2] > 0 && coord_dest[2] <= xlength) {
      
      cell_ptr = getEl(collideField, x, y, z, i, length_tot);
      *cell_ptr= *getEl(collideField, coord_dest[0], coord_dest[1], coord_dest[2], 18-i, length_tot);

      if (flagField[length_tot*length_tot*z+length_tot*y+x] == 2 ) {
        dotProd = 0;
        for (int j = 0; j < 3; ++j)	{
          dotProd += LATTICEVELOCITIES[i][j] * wallVelocity[j];
        }	
        
        computeDensity(getEl(collideField, coord_dest[0], coord_dest[1], coord_dest[2], 0, length_tot), &density);
        *cell_ptr += 2.0 * LATTICEWEIGHTS[i] * dotProd / (C_S*C_S);
      }
    }
  }
}

void treatBoundary(double *collideField, int *flagField, const double * const wallVelocity, int xlength){
	int x, y, z;

	z = 0;
	for (y = 1; y <= xlength; ++y){
	    for (x = 1; x <= xlength; ++x){           	        
	        boundaryCell(collideField, flagField, wallVelocity, xlength, x, y, z);
	    }
	}
	for (z = 1; z <= xlength; ++z){
	    y = 0;
	    for (x = 1; x <= xlength; ++x){
	        boundaryCell(collideField, flagField, wallVelocity, xlength, x, y, z);
	    }       
	    for (y = 1; y <= xlength; ++y){
	        x = 0;
	        boundaryCell(collideField, flagField, wallVelocity, xlength, x, y, z);
	        x = xlength+1;
	        boundaryCell(collideField, flagField, wallVelocity, xlength, x, y, z);
	    }
	    y = xlength+1;
	    for (x = 1; x <= xlength; ++x){
	        boundaryCell(collideField, flagField, wallVelocity, xlength, x, y, z);
	    }
	}
	z = xlength+1;
	for (y = 1; y <= xlength; ++y){
	    for (x = 1; x <= xlength; ++x){
	        boundaryCell(collideField, flagField, wallVelocity, xlength, x, y, z);
	    }
	}
}
