#include "parallel.h"
#include <mpi.h>

#include "LBDefinitions.h"
#include "initLB.h"
#include "computeCellValues.h"

/* Debugging printing for velocities */
/* TODO move to helper */
void printXZvelocities(double * collideField, int y, int * n) {
    int x, z, node[3];
    double density, velocity[3], * el;

    node[1] = y;

    for(z = 0; z < n[2]; z++) {
        node[2] = z;
        for(x = 0; x < n[0]; x++) {
            node[0] = x;
            el = getEl(collideField, node, 0, n);
            computeDensity(el, &density);
            computeVelocity(el, &density, velocity);
            printf("(%f %f %f) ", velocity[0], velocity[1], velocity[2]);
            //printf("%f ", el[2]);
        }
        printf("\n");
    }
    printf("\n");
}

/* Debugging printing for buffer */
void printBuffer(double * buffer, int n1, int n2) {
    int i, j;

    for (int k = 0; k < q; k++) {
    for (j = 0; j < n2; j++) {
        for (i = 0; i < n1; i++) {
            printf("%f ", *getBufferEl(buffer, i, j, k, n1));
        }
        printf("\n");
    }
    printf("\n");
    }
    printf("\n");
}

void Program_Message(char *txt)
/* produces a stderr text output  */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   fprintf(stderr,"-MESSAGE- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
}


void Programm_Sync(char *txt)
/* produces a stderr textoutput and synchronize all processes */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Barrier(MPI_COMM_WORLD);                             /* synchronize output */  
   fprintf(stderr,"-MESSAGE- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
   MPI_Barrier(MPI_COMM_WORLD);
}

void Programm_Stop(char *txt)
/* all processes will produce a text output, be synchronized and finished */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Barrier(MPI_COMM_WORLD);                           /* synchronize output */
   fprintf(stderr,"-STOP- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();
   exit(1);
}

/* Initializes the MPI Porcess, also obtains the total number of ranks and the process assigned */
void initializeMPI(int * my_rank,int * number_of_ranks,int argc, char * argv[]){
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, number_of_ranks );
    MPI_Comm_rank( MPI_COMM_WORLD, my_rank );
}

/*
  Extraction for TOP & BOTTOM
  node[2] = 0 || length[2] + 1
*/
/* TODO Extract_<direction> is also too copy-paste */
void extract_z(double * field,
               double * sendBuffer,
               const int * lattices, /* lattices which will be sent */
               int * node,
               int * n) {
    int x, y, i, index;

    for (y = 1; y <= n[1] - 2; y++) {
        node[1] = y;
        for (x = 0; x < n[0]; x++) {
            node[0] = x;
            for (i = 0; i < q; i++) {
                index = lattices[i];
                *getBufferEl(sendBuffer, x, y, i, n[0]) = *getEl(field, node, index, n);
            }
        }
    }
}

/*
  Extraction for FRONT & BACK
  node[1] = 0 || length[1] + 1
*/
void extract_y(double * field,
               double * sendBuffer,
               const int * lattices, /* lattices which will be sent */
               int * node,
               int * n) {
    int x, z, i, index;

    for (z = 0; z < n[2]; z++) {
        node[2] = z;
        for (x = 0; x < n[0]; x++) {
            node[0] = x;
            for (i = 0; i < q; i++) {
                index = lattices[i];
                *getBufferEl(sendBuffer, x, z, i, n[0]) = *getEl(field, node, index, n);
            }
        }
    }
}

/*
  Extraction for LEFT & RIGHT
  node[0] = 0 || length[0] + 1
*/
void extract_x(double * field,
               double * sendBuffer,
               const int * lattices, /* lattices which will be sent */
               int * node,
               int * n
    ) {
    int y, z, i, index;
    int yLength = n[1] - 2;

    for (z = 1; z <= n[2] - 2; z++) {
        node[2] = z;
        for (y = 1; y <= yLength; y++) {
            node[1] = y;
            for (i = 0; i < q; i++) {
                index = lattices[i];
                *getBufferEl(sendBuffer, y, z, i, n[1]) = *getEl(field, node, index, n);
            }
        }
    }
}

void swap(face_t face, double * sendBuffer, double * readBuffer, int count, int destination, int my_rank) {
    // MPI send to neighboring process, block until receives 
    MPI_Status status;

    MPI_Sendrecv(sendBuffer, count, MPI_DOUBLE, destination, 0,
                 readBuffer, count, MPI_DOUBLE, destination, 0, MPI_COMM_WORLD, &status);
}

/*
  Injection for TOP & BOTTOM
  node[2] = 1 || length[2]
*/
void inject_z(double * field,  double * readBuffer, const int * lattices, int * node, int * n, int my_rank) {
    int x, y, i, index;

    for (y = 1; y <= n[1] - 2; y++) {
        node[1] = y;
        for (x = 0; x < n[0]; x++) {
            node[0] = x;
            for (i = 0; i < q; i++) {
                index = lattices[i];
                *getEl(field, node, index, n) = *getBufferEl(readBuffer, x, y, i, n[0]);
            }
        }
    }
}

/*
  Injection for FRONT & BACK
  node[1] = 1 || length[1]
*/
void inject_y(double * field,  double * readBuffer, const int * lattices, int * node, int * n) {
    int x, z, i, index;

    for (z = 0; z < n[2]; z++) {
        node[2] = z;
        for (x = 0; x < n[0]; x++) {
            node[0] = x;
            for (i = 0; i < q; i++) {
                index = lattices[i];
                *getEl(field, node, index, n) = *getBufferEl(readBuffer, x, z, i, n[0]);
            }
        }
    }
}

/*
  Injection for LEFT & RIGHT
  node[0] = 1 || length[0]
*/
void inject_x(double * field,  double * readBuffer, const int * lattices, int * node, int * n) {
    int y, z, i, index;
    int yLength = n[1] - 2;

    for (z = 1; z <= n[2] - 2; z++) {
        node[2] = z;
        for (y = 1; y <= yLength; y++) {
            node[1] = y;
            for (i = 0; i < q; i++) {
                index = lattices[i];
                *getEl(field, node, index, n) = *getBufferEl(readBuffer, y, z, i, n[1]);
            }
        }
    }
}

void exchange(face_t face,
              double * field,
              double * sendBuffer,
              double * readBuffer,
              int * length,
              int * my_pos,
              int * Proc,
              int my_rank
    ) {
    int n[3] = { length[0] + 2, length[1] + 2, length[2] + 2 };
    int node[3];
    int count;
    int destination;

    /* TODO It would be really cool to get rig of this copy-paste */
    if (face == LEFT) {
        /* check existance of a neighbor */
        if (my_pos[0] != 0) {
            /* Number of elements in a buffer */
            count = n[1] * n[2] * q;
            /* rank of neighbor process */
            destination = get_rank(my_pos[0] - 1, my_pos[1], my_pos[2], Proc);
 
            node[0] = 1;
            extract_x(field, sendBuffer, left_lattices, node, n);
           
            swap(face, sendBuffer, readBuffer, count, destination, my_rank);
 
            node[0] = 0;
            inject_x(field, readBuffer, right_lattices, node, n);
            //printXZvelocities(field, 2, n);
       }
    } else if (face == RIGHT) {
        /* check existance of a neighbor */
        if (my_pos[0] != Proc[0] - 1) {
            /* Number of elements in a buffer */
            count = n[1] * n[2] * q;
            /* rank of neighbor process */
            destination = get_rank(my_pos[0] + 1, my_pos[1], my_pos[2], Proc);
 
            node[0] = length[0];
            extract_x(field, sendBuffer, right_lattices, node, n);

            swap(face, sendBuffer, readBuffer, count, destination, my_rank);

            node[0] = length[0] + 1;
            inject_x(field, readBuffer, left_lattices, node, n);
        }
    } else if (face == TOP) {
        /* check existance of a neighbor */
        if (my_pos[2] != Proc[2] - 1) {
            /* Number of elements in a buffer */
            count = n[0] * n[1] * q;
            /* rank of neighbor process */
            destination = get_rank(my_pos[0], my_pos[1], my_pos[2] + 1, Proc);

            node[2] = length[2];
            extract_z(field, sendBuffer, top_lattices, node, n);
 
            swap(face, sendBuffer, readBuffer, count, destination, my_rank);
            
            node[2] = length[2] + 1;
            inject_z(field, readBuffer, bottom_lattices, node, n, my_rank);
        }
    } else if (face == BOTTOM) {
        /* check existance of a neighbor */
        if (my_pos[2] != 0) {
            /* Number of elements in a buffer */
            count = n[0] * n[1] * q;
            /* rank of neighbor process */
            destination = get_rank(my_pos[0], my_pos[1], my_pos[2] - 1, Proc);
            
            node[2] = 1;
            extract_z(field, sendBuffer, bottom_lattices, node, n);
 
            swap(face, sendBuffer, readBuffer, count, destination, my_rank);

            node[2] = 0;
            inject_z(field, readBuffer, top_lattices, node, n, my_rank);
        }
    } else if (face == FRONT) {
        /* check existance of a neighbor */
        if (my_pos[1] != 0) {
            /* Number of elements in a buffer */
            count = n[0] * n[2] * q;
            /* rank of neighbor process */
            destination = get_rank(my_pos[0], my_pos[1] - 1, my_pos[2], Proc);
 
            node[1] = 1;
            extract_y(field, sendBuffer, front_lattices, node, n);
 
            swap(face, sendBuffer, readBuffer, count, destination, my_rank);

            node[1] = 0;
            inject_y(field, readBuffer, back_lattices, node, n);
        }
    } else if (face == BACK) {
        /* check existance of a neighbor */
        if (my_pos[1] != Proc[1] - 1) {
            /* Number of elements in a buffer */
            count = n[0] * n[2] * q;
            /* rank of neighbor process */
            destination = get_rank(my_pos[0], my_pos[1] + 1, my_pos[2], Proc);
 
            node[1] = length[1];
            extract_y(field, sendBuffer, back_lattices, node, n);
 
            swap(face, sendBuffer, readBuffer, count, destination, my_rank);

            node[1] = length[1] + 1;
            inject_y(field, readBuffer, front_lattices, node, n);
        }
    }
}

