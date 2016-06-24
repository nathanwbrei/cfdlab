#include "helper.h"
#include "computeCellValues.h"
#include "LBDefinitions.h"
#include "checks.h"

void check_in_rank(double *collideField, int *flagField, int * length, int t){
	int node[D], i, x, y, z;
	int n[D] = { length[0] + 2, length[1] + 2, length[2] + 2 };
	double velocity[D], density, norm_v=0, *currentCell;

	for (z = 1; z < length[2]; ++z){
		node[2] = z;
		for (y = 1; y < length[1]; ++y){
			node[1] = y;
			for (x = 1; x < length[0]; ++x){
				node[0] = x;
				if (*getFlag(flagField, node, n) == FLUID){
					currentCell = getEl(collideField, node, 0, n);
					computeDensity(currentCell, &density);
					computeVelocity(currentCell, &density, velocity);
					
					for (i = 0; i < D; ++i)
						norm_v += velocity[i]*velocity[i];

					if(density>1.1 || density<0.9)
						printf("Warning: In timestep %d position %d %d %d is an anormal density of %f \n", t, node[0],node[1],node[2],density);

					if(norm_v>(3 * C_S * C_S))
						printf("Warning: In timestep %d position %d %d %d is an anormal velocity of %f \n", t, node[0],node[1],node[2],norm_v);
				}
			}
		}
	}
}

void check_flags(int * flagField, int* length, int flag1, int flag2, int t){
	int node[D], node2[D], i, j, x, y, z;
	int n[D] = { length[0] + 2, length[1] + 2, length[2] + 2 };

	for (z = 1; z < length[2]; ++z){
		node[2] = z;
		for (y = 1; y < length[1]; ++y){
			node[1] = y;
			for (x = 1; x < length[0]; ++x){
				node[0] = x;
				if (*getFlag(flagField, node, n) == flag1){
					for (i = 0; i < Q; ++i){
						for (j = 0; j < D; ++j)
							node2[j]=node[j]+LATTICEVELOCITIES[i][j];
						if (*getFlag(flagField, node, n) == flag2)
							printf("Warning: In timestep %d position %d %d %d a cell %d is adjacent to a forbiden cell %d in %d %d %d \n", t, node[0],node[1],node[2],flag1, node2[0],node2[1],node2[2], flag2);
					}
				}
			}
		}
	}
}

void check_mass(double *massField, int* flagField, int* length, int t){
	double tot_mass = 0;
	int node[D], x, y, z;
	int n[D] = { length[0] + 2, length[1] + 2, length[2] + 2 };

	for (z = 1; z < length[2]; ++z){
		node[2] = z;
		for (y = 1; y < length[1]; ++y){
			node[1] = y;
			for (x = 1; x < length[0]; ++x){
				node[0] = x;
				if (*getFlag(flagField, node, n) == FLUID){
					tot_mass += *getMass(massField, node, n);
				}
			}
		}
	}
	printf("On timestep %d : Total mass: %f \n",t,tot_mass);
}

void run_checks(double *collideField, double *massField, int *flagField, int * length, int t ){
//	check_in_rank(collideField, flagField, length, t);
//	check_flags(flagField, length, FLUID, GAS, 	 t);
	check_mass( massField, flagField, length, t);
}
