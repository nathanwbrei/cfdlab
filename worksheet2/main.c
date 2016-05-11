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
    double tau, velocityWall;

    readParameters(&xlength, &tau, &velocityWall, &timesteps, &timestepsPerPlotting, argc, argv);

    return 0;
}

#endif

