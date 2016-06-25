#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>

#ifndef PARALLEL_H
#define PARALLEL_H

void Program_Message(char *txt);
/* produces a stderr text output  */



void Programm_Sync(char *txt);
/* produces a stderr textoutput and synchronize all processes */



void Programm_Stop(char *txt);
/* all processes will produce a text output, be synchronized and finished */

void initializeMPI(int * my_rank,int *number_of_ranks,int argc, char *argv[]);
/* Initializes the MPI Porcess, also obtains the total number of ranks and the process assigned */

typedef enum {LEFT, RIGHT, TOP, BOTTOM, FRONT, BACK} face_t;

void exchange(face_t face,
              double * field,
              double * sendBuffer,
              double * readBuffer,
              int * length,
              int * my_pos,
              int * Proc,
              int my_rank
);

/* Obtain the position in the D dimentional space my_pos, given the process number rank amd the number of process per dimention (Proc)*/
void get_rank_pos(int * my_pos, int rank, int *Proc);

/* Obtain the size of the prtion assigned to the procces my_lengths, using the length of the cavity xlength, the position on each dimetion y the total number of processes per dimention*/
void get_my_lengths(int* my_pos, int * length, int * Proc, int * my_lengths, int * my_origin);

int get_rank(int px, int py, int pz,  int * Proc);

/* Debugging printing for buffer */
void printBuffer(double * buffer, int n1, int n2);

void printXZvelocities(double * field, int y, int * n);
#endif
