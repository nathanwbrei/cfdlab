/*
 * Apply the boundary values to the velocity matrices U, V in-place. 
 * TODO: This constrains our code to only solving cavity problems!
 */
void boundaryvalues(int imax, int jmax, double **U, double **V) {

    /*
     * Suppose our full staggered grid has shape 0:imax+1 x 0:jmax+1.
     * Then U has shape 0:imax x 0:jmax+1
     *      V has shape 0:imax+1 x 0:jmax
     *      P has shape 0:imax+1 x 0:jmax+1
    
     * To populate the full boundary of U, 
     *    we first move along the sides i<-{0,imax} _excluding_ corner points j<-{0,jmax+1}
     *    and then move along the top/bottom, _including_ corner points
     *    This is necessary because the corner values depend on their
     *    nearest side neighbors.
    
     * We do something similar for V.
    
     * Note that the assignment sheet incorrectly specifies the range for the U[i][1]=-U[i][0], etc.
     *    i needs to range over 0:imax, NOT 1:imax, otherwise we miss two corner values. 
     *    These values are constrained to be 0 anyway and might not even be used, 
     *    but this will come back to bite us if we ever try to make this program more general.
    
     * Set u=0 along sides
     */

    int i, j;

    for (j=1; j<=jmax; j++) {
        U[0][j] = 0;
        U[imax][j] = 0;
    }
    /* Set v=0 along top, bottom */
    for (i=1; i<=imax; i++) {
        V[i][0] = 0;
        V[i][jmax] = 0;
    }
    /* Set v along sides */
    for (j=0; j<=jmax; j++) {   /* Note the 0 */
        V[0][j] = -V[1][j];
        V[imax+1][j] = -V[imax][j];
    }
    /* Set u along top and bottom */
    for (i=0; i<=imax; i++) {  /* Note the 0 */
        U[i][0] = -U[i][1];
        U[i][jmax+1] = 2.0 - U[i][jmax];
    }
}
