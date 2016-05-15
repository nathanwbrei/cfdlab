#ifndef _MAIN_C_
#define _MAIN_C_

#include "collision.h"
#include "streaming.h"
#include "initLB.h"
#include "visualLB.h"
#include "boundary.h"

int main(int argc, char *argv[]){
  /* TODO */
    int xlength, timesteps, timestepsPerPlotting;
    double tau, velocityWall[3];    

    readParameters(&xlength, &tau, velocityWall, &timesteps, &timestepsPerPlotting, argc, argv);

    /* TODO: Add errors to failed allocation */
    double  *collideField = (double *)  malloc((size_t)( 19*(xlength+2)*(xlength+2)*(xlength+2) ) * sizeof( double ));
    double  *streamField = (double *)  malloc((size_t)( 19*(xlength+2)*(xlength+2)*(xlength+2) ) * sizeof( double ));
    int  *flagField = (int *)  malloc((size_t)( (xlength+2)*(xlength+2)*(xlength+2) ) * sizeof( int ));

    initialiseFields( collideField, streamField, flagField, xlength);

    //treatBoundary(collideField,flagfield,velocityWall,xlength);

    free(collideField);
    free(streamField);
    free(flagField);
    
    return 0;
}

#endif

