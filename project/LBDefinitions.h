#include <math.h>

#ifndef _LBDEFINITIONS_H_
#define _LBDEFINITIONS_H_

#define Q 19
#define q 5
#define D 3

static const int LATTICEVELOCITIES[Q][D] = {
    {0, -1, -1}, {-1, 0, -1}, {0, 0, -1}, {1, 0, -1}, {0, 1, -1},
    {-1, -1, 0}, {0, -1, 0}, {1, -1, 0}, {-1, 0, 0}, {0, 0, 0}, {1, 0, 0}, {-1, 1, 0}, {0, 1, 0}, {1, 1, 0},
    {0, -1, 1}, {-1, 0, 1}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}
};

static const double LATTICEWEIGHTS[Q] = {
    1.0/36, 1.0/36, 2.0/36, 1.0/36, 1.0/36,
    1.0/36, 2.0/36, 1.0/36, 2.0/36, 12.0/36, 2.0/36, 1.0/36, 2.0/36, 1.0/36, 
    1.0/36, 1.0/36, 2.0/36, 1.0/36, 1.0/36
};

/*
  Lattices which will be sent to corresponding neighbor
*/
static const int right_lattices[q] = { 3, 7, 10, 13, 17};
static const int left_lattices[q] =  { 1, 5, 8, 11, 15 };
static const int top_lattices[q] = { 14, 15, 16, 17, 18 };
static const int bottom_lattices[q] = { 0, 1, 2, 3, 4 };
static const int front_lattices[q] = { 0, 5, 6, 7, 14 };
static const int back_lattices[q] = { 4, 11, 12, 13, 18 };

/* Math.sqrt(3.0) = 1.732050807568877 */
static const double C_S = 1.0 / 1.732050807568877;
static const double C_S2 = 1.0 / 3.0; 

static const int FLUID = 0;
static const int NOSLIP = 1;
static const int MOVING_WALL = 2;
static const int INFLOW = 3;
static const int OUTFLOW = 4;
static const int PRESSURE_IN = 5;
static const int OBSTACLE = 6;
static const int FREESLIP = 7;
static const int INTERFACE = 8;
static const int GAS = 9;
static const int PARALLEL = 10;

static const char * const PARABOLIC_SCENARIO = "parabolic";
static const char * const CONSTANT_SCENARIO = "constant";

#endif

