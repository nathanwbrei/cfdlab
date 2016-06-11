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
              double ** sendBuffer,
              double ** readBuffer,
              int * length,
              int * my_pos,
              int * Proc,
              int my_rank
);
#endif
