#ifndef _MASS_H_
#define _MASS_H_

typedef enum { STANDART_CELL, NO_FLUID_NB, NO_EMPTY_NB } CellType;

void calculateMassExchange(double * collideField, double * massField, double * fractionField, int * length); 

#endif
