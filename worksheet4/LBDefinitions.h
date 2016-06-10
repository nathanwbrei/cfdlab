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

static const int right_lattices[q] = { 1, 7, 9, 11, 13 };
static const int left_lattices[q] = { 2, 8, 10, 12, 14 };
static const int top_lattices[q] = { 3, 5, 7, 8, 17 };
static const int bottom_lattices[q] = { 4, 9, 10, 16, 18 };
static const int front_lattices[q] = { 6, 13, 14, 17, 18 };
static const int back_lattices[q] = { 5, 11, 12, 15, 16 };

/* Math.sqrt(3.0) = 1.732050807568877 */
static const double C_S = 1.0 / 1.732050807568877;
/* C_S * C_S */
static const double C_S_2 = 1.0 / 3;

static const int FLUID = 0;
static const int NOSLIP = 1;
static const int MOVING_WALL = 2;
static const int PARALLEL = 10;

#endif

