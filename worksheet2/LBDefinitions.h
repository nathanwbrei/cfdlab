#include <math.h>

#ifndef _LBDEFINITIONS_H_
#define _LBDEFINITIONS_H_

#define Q 19
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

/* Math.sqrt(3.0) = 1.732050807568877 */
static const double C_S = 1.0 / 1.732050807568877;

#endif

