#include "boundary_val.h"

void boundaryvalues(
	  int imax,
	  int jmax,
	  double **U,
	  double **V){

	for (int i = 1; i <= imax ; ++i)	{
		V[i][0] = 0;
		V[i][jmax] = 0;
		U[i][0]=-U[i][1];
		U[i][jmax+1] = -U[i][jmax];
		//P[i][0] = P[i][1];
		//P[i][jmax+1] = P[i][jmax];
	}
	for (int j = 0; j <= jmax; ++j)
	{
		U[0][j] = 0;
		U[imax][j] = 0;
		V[0][j] = -V[1][j];
		V[imax +1] [j] = -V[imax][j];
		//P[0][j] = P[1][j];
		//P[imax+1][j] = P[imax][j];
	}

	//return void;
}