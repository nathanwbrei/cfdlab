#include "initLB.h"

int readParameters(int *xlength, double *tau, double *velocityWall, int *timesteps, int *timestepsPerPlotting, int argc, char *argv[]){

	char * FileName;
	
	if (argc == 2) {
	    FileName = argv[1];	    
	}
	else{
	    printf("Usage: ./lbsim inputfile \n");
	    return -1;
	}

	READ_INT   ( FileName, *xlength );

	READ_DOUBLE( FileName,  *tau);
	READ_DOUBLE( FileName,  *velocityWall);

	READ_INT   ( FileName, *timesteps );
	READ_INT   ( FileName, *timestepsPerPlotting );

	return 0;
}


void initialiseFields(double *collideField, double *streamField, int *flagField, int xlength){
  /* TODO */
}

