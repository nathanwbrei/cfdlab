#include "initLB.h"

int readParameters(int *xlength, double *tau, double *velocityWall, int *timesteps, int *timestepsPerPlotting, int argc, char *argv[]){
  /* TODO */



  return 0;
}


/* TODO: chech if colidefield or streamfield or both */
/* TODO: define flux for boundary cells*/
/* TODO: put all inside initialiceCell and inline it*/
void initialiseCell(double *collideField, double *streamField, int *flagField, int length_tot, int x, int y, int z, int flag){
	int pos_flag = length_tot*length_tot*z+length_tot*y + x;
	int pos_field = 19 * pos_flag;
	const double w0=12.0/36.0, w1=2.0/36.0, w2=1.0/36.0;
	flagField[pos_flag] = flag;

	collideField[pos_field+0] = w2; 
	collideField[pos_field+1] = w2; 
	collideField[pos_field+2] = w1; 
	collideField[pos_field+3] = w2; 
	collideField[pos_field+4] = w2; 
	collideField[pos_field+5] = w2; 
	collideField[pos_field+6] = w1; 
	collideField[pos_field+7] = w2; 
	collideField[pos_field+8] = w1; 
	collideField[pos_field+9] = w0; 
	collideField[pos_field+10] = w1; 
	collideField[pos_field+11] = w2; 
	collideField[pos_field+12] = w1; 
	collideField[pos_field+13] = w2; 
	collideField[pos_field+14] = w2; 
	collideField[pos_field+15] = w2; 
	collideField[pos_field+16] = w1; 
	collideField[pos_field+17] = w2; 
	collideField[pos_field+18] = w2; 


}

void initialiseFields(double *collideField, double *streamField, int *flagField, int xlength){
  
  	int x, y, z, length_tot;
  	length_tot = xlength +2 ;
  	/* Definition of the fields */
  	/* TODO: There are some overkills or cells not written in the code, think carefully in the boundaries*/
  	z = 0;
	for (y = 0; y <= xlength+1; ++y){
		for (x = 0; x <= xlength+1; ++x){			
			initialiseCell(collideField, streamField, flagField, length_tot, x, y, z, 1);
		}
	}
	for (z = 1; z <= xlength; ++z){
		y = 0;
		for (x = 0; x <= xlength +1; ++x){
			initialiseCell(collideField, streamField, flagField, length_tot, x, y, z, 1);
		}		
		for (y = 1; y <= xlength; ++y){
			x = 0;
			initialiseCell(collideField, streamField, flagField, length_tot, x, y, z, 1);
			for (x = 1; x <= xlength; ++x){
				flagField[length_tot*length_tot*z+length_tot*y + x] = 0;
				initialiseCell(collideField, streamField, flagField, length_tot, x, y, z, 0);
			}
			x = xlength+1;
			initialiseCell(collideField, streamField, flagField, length_tot, x, y, z, 1);
		}
		y = xlength+1;
		for (x = 0; x <= xlength+1; ++x){
			initialiseCell(collideField, streamField, flagField, length_tot, x, y, z, 1);
		}
	}
	z = xlength+1;
	for (y = 0; y <= xlength+1; ++y){
		for (x = 0; x <= xlength+1; ++x){
			initialiseCell(collideField, streamField, flagField, length_tot, x, y, z, 2);
		}
	}


}

