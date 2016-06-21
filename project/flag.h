#ifndef _FLAG_H_
#define _FLAG_H_

void updateFlagField(double * collideField, int * flagField, double * fractionField, int * length);
/* TODO change emptiedCells type to pointer */
void updateEmptiedNeighbors(double * collideField, int * flagField, int emptiedCells[][3], int nEmptied, int * n);

#endif
