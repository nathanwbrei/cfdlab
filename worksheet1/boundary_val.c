


/*
 * Apply the boundary values to the velocity matrices U, V in-place. 
 * TODO: This constrains our code to only solving cavity problems!
 */
void boundaryvalues(int imax, int jmax, double **U, double **V) {
 
    // Right and left boundaries
    for (int j=1; j<=jmax; j++) {
        // No-slip along both
        // u(0,j) = 0 for all j
        // u(imax,j) = 0 for all j
        U[0][j] = 0;
        V[0][j] = -V[1][j];

        //TODO: Are these indices correct?
        U[imax+1][j] = 0;
        V[imax+1][j] = -V[imax][j];

    }

    // Upper and lower boundaries
    for (int i=1; i<=imax; i++) {
        // Moving-wall along upper boundary: u = u_wall; v = 0
        // u(i,jmax+1) = -u(i,jmax) + 2*u_wall, for all i<-1:imax, where u_wall=1 (See p. 16)
        U[i][jmax+1] = -U[i][jmax] + 2;
        V[i][jmax+1] = 0;

        // No-slip along lower boundary: u = v = 0
        U[i][0] = -U[i][1];
        V[i][0] = 0;
    }

    


}

