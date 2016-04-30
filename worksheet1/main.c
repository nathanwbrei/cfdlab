#include "helper.h"
#include "visual.h"
#include "init.h"
#include "boundary_val.h"
#include "uvp.h"
#include "sor.h"
#include <stdio.h>



/**
 * The main operation reads the configuration file, initializes the scenario and
 * contains the main loop. So here are the individual steps of the algorithm:
 *
 * - read the program configuration file using read_parameters()
 * - set up the matrices (arrays) needed using the matrix() command
 * - create the initial setup init_uvp(), init_flag(), output_uvp()
 * - perform the main loop
 * - trailer: destroy memory allocated and do some statistics
 *
 * The layout of the grid is decribed by the first figure below, the enumeration
 * of the whole grid is given by the second figure. All the unknowns corresond
 * to a two dimensional degree of freedom layout, so they are not stored in
 * arrays, but in a matrix.
 *
 * @image html grid.jpg
 *
 * @image html whole-grid.jpg
 *
 * Within the main loop the following big steps are done (for some of the 
 * operations a definition is defined already within uvp.h):
 *
 * - calculate_dt() Determine the maximal time step size.
 * - boundaryvalues() Set the boundary values for the next time step.
 * - calculate_fg() Determine the values of F and G (diffusion and confection).
 *   This is the right hand side of the pressure equation and used later on for
 *   the time step transition.
 * - calculate_rs()
 * - Iterate the pressure poisson equation until the residual becomes smaller
 *   than eps or the maximal number of iterations is performed. Within the
 *   iteration loop the operation sor() is used.
 * - calculate_uv() Calculate the velocity at the next time step.
 */
int main(int argn, char** args){

    /*
     * Define the problem parameters
     */
    /*
     * Geometry parameters
     */
    double xlength, ylength, dx, dy;
    int imax, jmax;
    /*
     * Timestep parameters
     */
    double t, t_end, dt, tau, dt_value;
    int n;
    /*
     * Pressure iteration data
     */
    double eps, omg, alpha, res;
    int itermax, it;
    /*
     * Problem-dependent quantities
     */
    double Re, GX, GY, UI, VI, PI;
    /*
     * Data structures
     */
    double **U, **V, **P, **RS, **F, **G;

    /*
      Read the problem parameters
    */
    char *filename = "cavity100.dat"; 
    if (argn == 2) {
        filename = args[1];
    }
           
    printf("Loading parameters from %s\n",filename);
    read_parameters(filename, &Re, &UI, &VI, &PI, &GX, &GY, &t_end, &xlength, 
                    &ylength, &dt, &dx, &dy, &imax, &jmax, &alpha, &omg, &tau, 
                    &itermax, &eps, &dt_value); 

    n = 0;
    t = 0;
    /*
     * Assign initial values to u,v,p
     */
    U = matrix(0, imax+1, 0, jmax+1);
    V = matrix(0, imax+1, 0, jmax+1);
    P = matrix(0, imax+1, 0, jmax+1);
    F = matrix(0, imax+1, 0, jmax+1); /*TODO: Definitely not the right indices*/
    G = matrix(0, imax+1, 0, jmax+1);
    RS = matrix(0, imax+1, 0, jmax+1);

    print_matrix("Created U", 0, imax+1, 0, jmax+1, U);

    init_uvp(UI, VI, PI, imax, jmax, U, V, P);
    print_matrix("Initialized U", 0, imax+1, 0, jmax+1, U);

    while (t < t_end) {
        calculate_dt(Re,tau,&dt,dx,dy,imax,jmax,U,V); 
        printf("@t=%f, dt=%f\n", t, dt);
        boundaryvalues(imax,jmax,U,V);
        /*print_matrix("Boundaryvals for U", 0, imax+1, 0, jmax+1, U);
        print_matrix("Boundaryvals for V", 0, imax+1, 0, jmax+1, V); */
        calculate_fg(Re,GX,GY,alpha,dt,dx,dy,imax,jmax,U,V,F,G);
        calculate_rs(dt,dx,dy,imax,jmax,F,G,RS);
        it = 0;
        res = eps+1; /* Residual starts off > eps */
        while (it < itermax && res > eps) {
            sor(omg,dx,dy,imax,jmax,P,RS,&res);
            it += 1;
        }
        calculate_uv(dt,dx,dy,imax,jmax,U,V,F,G,P);
        /*TODO: Condition this on some parameter*/
        write_vtkFile(args[0], n, xlength, ylength, imax, jmax, dx, dy, U, V, P);
        t += dt;
        n += 1;
    }
    /*
      Output u,v,p for visualization
      TODO: Do we want this?
    */
    /*
     * Free the matrices we've allocated
     */
    free_matrix(U, 0, imax+1, 0, jmax+1);
    free_matrix(V, 0, imax+1, 0, jmax+1);
    free_matrix(P, 0, imax+1, 0, jmax+1);
    free_matrix(F, 0, imax+1, 0, jmax+1); /*TODO: Definitely the wrong indices*/
    free_matrix(G, 0, imax+1, 0, jmax+1);
    free_matrix(RS, 0, imax+1, 0, jmax+1);
     
    return -1;
}
