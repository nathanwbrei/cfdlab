#include "visualLB.h"
#include "computeCellValues.h"
#include "helper.h"

void write_vtkFile(const char *szProblem, int t, int * length, double * collideField, int * my_pos, int my_rank) {
    char szFileName[80];
    FILE *fp=NULL;
 
    int x, y, z;
    int node[3];
    int n[3] = { length[0] + 2, length[1] + 2, length[2] + 2 };
    double velocity[3];
    double density;
    double * el = NULL;

    sprintf( szFileName, "%s_%i.%i.vtk", szProblem, my_rank, t);
    fp = fopen( szFileName, "w");
    if( fp == NULL )		       
    {
        char szBuff[80];
        sprintf( szBuff, "Failed to open %s", szFileName );
        ERROR( szBuff );
        return;
    }

    write_vtkHeader(fp, length);
    write_vtkPointCoordinates(fp, length, my_pos);

    fprintf(fp,"POINT_DATA %i \n", length[0] * length[1] * length[2]);
  
    fprintf(fp,"\n");
    fprintf(fp, "VECTORS velocity float\n");
    
    for(z = 1; z  <= length[2]; z++) {
        node[2] = z;
        for(y = 1; y  <= length[1]; y++) {
            node[1] = y;
            for(x = 1; x  <= length[0]; x++) {
                node[0] = x;
                el = getEl(collideField, node, 0, n);
                computeDensity(el, &density);
                computeVelocity(el, &density, velocity);
                fprintf(fp, "%f %f %f\n", velocity[0], velocity[1], velocity[2]);
            }
        }
    }

    fprintf(fp,"\n");

    fprintf(fp, "SCALARS density float 1 \n"); 
    fprintf(fp, "LOOKUP_TABLE default \n");

    for(z = 1; z <= length[2]; z++) {
        node[2] = z;
        for(y = 1; y <= length[1]; y++) {
            node[1] = y;
            for(x = 1; x <= length[0]; x++) {
                node[0] = x;
                computeDensity(getEl(collideField, node, 0, n), &density);
                fprintf(fp, "%f\n", density);
            }
        }
    }

    if( fclose(fp) )
    {
        char szBuff[80];
        sprintf( szBuff, "Failed to close %s", szFileName );
        ERROR( szBuff );
    }
}


void write_vtkHeader( FILE *fp, int * length) {
    if( fp == NULL ) {
        char szBuff[80];
        sprintf( szBuff, "Null pointer in write_vtkHeader" );
    
        ERROR( szBuff );
    
        return;
    }

    fprintf(fp,"# vtk DataFile Version 2.0\n");
    fprintf(fp,"generated by CFD-lab course output (written by Tobias Neckel) \n");
    fprintf(fp,"ASCII\n");
    fprintf(fp,"\n");	
    fprintf(fp,"DATASET STRUCTURED_GRID\n");
    fprintf(fp,"DIMENSIONS  %i %i %i \n", length[0], length[1], length[2]);
    fprintf(fp,"POINTS %i float\n", length[0] * length[1] * length[2]);
    fprintf(fp,"\n");
}


void write_vtkPointCoordinates(FILE *fp, int * length, int * my_pos) {
    int originX = 0;
    int originY = 0;
    int originZ = 0;

    int x = my_pos[0] * length[0];
    int y = my_pos[1] * length[1];
    int z = my_pos[2] * length[2];

    for(z = 1; z  <= length[2]; z++) {
        for(y = 1; y  <= length[1]; y++) {
            for(x = 1; x  <= length[0]; x++) {
                /* dx = dy = dz = 1 */
                fprintf(fp, "%d %d %d\n", originX + x, originY + y, originZ + z);
            }
        }
    }
}

void writeVtkOutput(double * collideField,
                    const int * const flagField,
                    const char * filename,
                    unsigned int t,
                    int * length,
                    int * my_pos,
                    int my_rank) {
    write_vtkFile(filename, t, length, collideField, my_pos, my_rank);
}

