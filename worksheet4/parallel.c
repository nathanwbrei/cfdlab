#include "parallel.h"
#include <mpi.h>

#include "LBDefinitions.h"
#include "initLB.h"


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

void initializeMPI(int * my_rank,int * number_of_ranks,int argc, char * argv[]){
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, number_of_ranks );
    MPI_Comm_rank( MPI_COMM_WORLD, my_rank );
}

void extract_z(double * field,  double * sendBuffer, const int * lattices, int * bufferLength, int * node, int * n) {
    int x, y, i, index;

    for (y = 1; y <= bufferLength[1]; y++) {
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

void extract_y(double * field,  double * sendBuffer, const int * lattices, int * bufferLength, int * node, int * n) {
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

void extract_x(double * field,  double * sendBuffer, const int * lattices, int * bufferLength, int * node, int * n) {
    int y, z, i, index;

    for (z = 1; z <= bufferLength[1]; z++) {
        node[2] = z;
        for (y = 1; y <= bufferLength[0]; y++) {
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
    // TODO: Needs to know who its neighbors are / whether they exist


//    MPI_Send(sendBuffer, count, MPI_DOUBLE, destination, 0, MPI_COMM_WORLD);
//    MPI_Recv(readBuffer, count, MPI_DOUBLE, destination, 0, MPI_COMM_WORLD, &status);
}

void inject_z(double * field,  double * readBuffer, const int * lattices, int * bufferLength, int * node, int * n) {
    int x, y, i, index;

    for (y = 1; y <= bufferLength[1]; y++) {
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

void inject_y(double * field,  double * readBuffer, const int * lattices, int * bufferLength, int * node, int * n) {
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

void inject_x(double * field,  double * readBuffer, const int * lattices, int * bufferLength, int * node, int * n) {
    int y, z, i, index;

    for (z = 1; z <= bufferLength[1]; z++) {
        node[2] = z;
        for (y = 1; y <= bufferLength[0]; y++) {
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
              double ** sendBuffer,
              double ** readBuffer,
              int * length,
              int * my_pos,
              int * Proc,
              int my_rank
    ) {
    int n[3] = { length[0] + 2, length[1] + 2, length[2] + 2 };
    int node[3];
    int bufferLength[2];
    int count;
    int destination;

    if (face == LEFT) {
        if (my_pos[0] != 0) {
            printf("%i Left\n", my_rank);
            bufferLength[0] = length[1];
            bufferLength[1] = length[2];
            count = n[1] * n[2] * q;
            destination = get_rank(my_pos[0] - 1, my_pos[1], my_pos[2], Proc);
 
            node[0] = 1;
            extract_x(field, sendBuffer[0], left_lattices, bufferLength, node, n);
            
            swap(face, sendBuffer[0], readBuffer[0], count, destination, my_rank);

            node[0] = 0;
            inject_x(field, readBuffer[0], left_lattices, bufferLength, node, n);
       }
    } else if (face == RIGHT) {
        if (my_pos[0] != Proc[0] - 1) {
            printf("%i Right\n", my_rank);
            bufferLength[0] = length[1];
            bufferLength[1] = length[2];
            count = n[1] * n[2] * q;
            destination = get_rank(my_pos[0] + 1, my_pos[1], my_pos[2], Proc);
 
            node[0] = length[0];
            extract_x(field, sendBuffer[1], right_lattices, bufferLength, node, n);

            swap(face, sendBuffer[1], readBuffer[1], count, destination, my_rank);

            node[0] = length[0] + 1;
            inject_x(field, readBuffer[1], right_lattices, bufferLength, node, n);
        }
    } else if (face == TOP) {
        if (my_pos[2] != Proc[2] - 1) {
            bufferLength[0] = length[0];
            bufferLength[1] = length[1];
            count = n[0] * n[1];
            destination = get_rank(my_pos[0], my_pos[1], my_pos[2] + 1, Proc);

            node[2] = length[2];
            extract_z(field, sendBuffer[2], top_lattices, bufferLength, node, n);
 
            swap(face, sendBuffer[2], readBuffer[2], count, destination, my_rank);
            
            node[2] = length[2] + 1;
            inject_z(field, readBuffer[2], top_lattices, bufferLength, node, n);
        }
    } else if (face == BOTTOM) {
        if (my_pos[2] != 0) {
            bufferLength[0] = length[0];
            bufferLength[1] = length[1];
            count = n[0] * n[1];
            destination = get_rank(my_pos[0], my_pos[1], my_pos[2] - 1, Proc);
            
            node[2] = 1;
            extract_z(field, sendBuffer[3], bottom_lattices, bufferLength, node, n);
 
            swap(face, sendBuffer[3], readBuffer[3], count, destination, my_rank);

            node[2] = 0;
            inject_z(field, readBuffer[3], bottom_lattices, bufferLength, node, n);
        }
    } else if (face == FRONT) {
        if (my_pos[1] != 0) {
            bufferLength[0] = length[0];
            bufferLength[1] = length[2];
            count = n[0] * n[2];
            destination = get_rank(my_pos[0], my_pos[1] - 1, my_pos[2], Proc);
 
            node[1] = 1;
            extract_y(field, sendBuffer[4], front_lattices, bufferLength, node, n);
 
            swap(face, sendBuffer[4], readBuffer[4], count, destination, my_rank);

            node[1] = 0;
            inject_y(field, readBuffer[4], front_lattices, bufferLength, node, n);
        }
    } else if (face == BACK) {
        if (my_pos[1] != Proc[1] - 1) {
            bufferLength[0] = length[0];
            bufferLength[1] = length[2];
            count = n[0] * n[2];
            destination = get_rank(my_pos[0], my_pos[1] + 1, my_pos[2], Proc);
 
            node[1] = length[1];
            extract_y(field, sendBuffer[5], back_lattices, bufferLength, node, n);
 
            swap(face, sendBuffer[5], readBuffer[5], count, destination, my_rank);

            node[1] = length[1] + 1;
            inject_y(field, readBuffer[5], back_lattices, bufferLength, node, n);
        }
    }
}

