#include "helper.h"

/**
 * Determines the value of U and G according to the formula
 *
 * @f$ F_{i,j} := u_{i,j} + \delta t \left( \frac{1}{Re} \left( \left[
    \frac{\partial^2 u}{\partial x^2} \right]_{i,j} + \left[
    \frac{\partial^2 u}{\partial y^2} \right]_{i,j} \right) - \left[
    \frac{\partial (u^2)}{\partial x} \right]_{i,j} - \left[
    \frac{\partial (uv)}{\partial y} \right]_{i,j} + g_x \right) @f$
 *
 * @f$ i=1,\ldots,imax-1, \quad j=1,\ldots,jmax @f$
 *
 * @f$ G_{i,j} := v_{i,j} + \delta t \left( \frac{1}{Re} \left(
   \left[ \frac{\partial^2 v}{\partial x^2}\right]_{i,j} + \left[ \frac{\partial^2 v}{\partial
                   y^2} \right]_{i,j} \right) - \left[ \frac{\partial
                   (uv)}{\partial x} \right]_{i,j} - \left[
                 \frac{\partial (v^2)}{\partial y} \right]_{i,j} + g_y
               \right) @f$
 *
 * @f$ i=1,\ldots,imax, \quad j=1,\ldots,jmax-1 @f$
 *
 */
void calculate_fg(
    double Re, double GX, double GY, double alpha, double dt, double dx, double dy,
    int imax, int jmax, double **U, double **V, double **F, double **G) {
    
    int i;
    int j;
    
    double u_frd_sum_i;
    double u_bkd_sum_i;
    double v_frd_sum_i;
    double v_bkd_sum_i;
    
    double u_frd_sum_j;
    double u_bkd_sum_j;
    double v_frd_sum_j;
    double v_bkd_sum_j;
    
    double d2ux;
    double d2uy;
    double du2x;
    double duvy;
    double d2vx;
    double d2vy;
    double duvx;
    double dv2y;
    
    for (j = 1; j < jmax + 1; j++) {
        F[0][j] = U[0][j];
        F[imax][j] = U[imax][j];
    }
    
    for (i = 1; i < imax + 1; i++) {
        G[i][0] = V[i][0];
        G[i][jmax] = V[i][jmax];
    }
    
    for (i = 1; i < imax; i++) {
        for (j = 1; j < jmax; j++) {
            u_frd_sum_i = (U[i][j] + U[i + 1][j]) / 2;
            u_frd_sum_j = (U[i][j] + U[i][j + 1]) / 2;
            
            u_bkd_sum_i = (U[i][j] + U[i - 1][j]) / 2;
            u_bkd_sum_j = (U[i][j] + U[i][j - 1]) / 2;
            
            v_frd_sum_i = (V[i][j] + V[i + 1][j]) / 2;
            v_frd_sum_j = (V[i][j] + V[i][j + 1]) / 2;
            
            v_bkd_sum_i = (V[i][j] + V[i - 1][j]) / 2;
            v_bkd_sum_j = (V[i][j] + V[i][j - 1]) / 2;
            
            d2ux = (U[i - 1][j] - 2 * U[i][j] + U[i + 1][j]) / (dx * dx);
            d2uy = (U[i][j - 1] - 2 * U[i][j] + U[i][j + 1]) / (dy * dy);
            
            du2x = (u_frd_sum_i * u_frd_sum_i - u_bkd_sum_i * u_bkd_sum_i) / dx +
            alpha * (abs(u_frd_sum_i) * (U[i][j] - U[i + 1][j]) / 2 -
                     abs(u_bkd_sum_i) * (U[i - 1][j] - U[i][j]) / 2) / dx;
            
            duvy = (v_frd_sum_i * u_frd_sum_j - (V[i][j - 1] + V[i + 1][j - 1]) / 2 * u_bkd_sum_j) / dy +
            alpha * (abs(v_frd_sum_i) * (U[i][j] - U[i][j - 1]) / 2 -
                     abs(V[i][j - 1] + V[i + 1][j - 1]) * (U[i][j - 1] - U[i][j]) / 4) / dy;
            
            d2vx = (V[i - 1][j] - 2 * V[i][j] + V[i + 1][j]) / (dx * dx);
            d2vy = (V[i][j - 1] - 2 * V[i][j] + V[i][j + 1]) / (dy * dy);
            
            dv2y = (v_frd_sum_j * v_frd_sum_j - v_bkd_sum_j * v_bkd_sum_j) / dy +
            alpha * (abs(v_frd_sum_j) * (V[i][j] - V[i][j + 1]) / 2 -
                     abs(v_bkd_sum_j) * (V[i][j - 1] - V[i][j]) / 2) / dy;
            
            duvx = (u_frd_sum_j * v_frd_sum_i - (U[i - 1][j] + U[i - 1][j + 1]) / 2 * v_bkd_sum_i) / dx +
            alpha * (abs(u_frd_sum_j) * (V[i][j] - V[i + 1][j]) / 2 -
                     abs(U[i - 1][j] + U[i - 1][j + 1]) * (V[i - 1][j] - V[i][j]) / 4) / dx;
            
            F[i][j] = U[i][j] + dt * ((d2ux + d2uy) / Re - du2x - duvy + GX);
            G[i][j] = V[i][j] + dt * ((d2vx + d2vy) / Re - duvx - dv2y + GY);
        }
        
        j = jmax;
        
        u_frd_sum_i = (U[i][j] + U[i + 1][j]) / 2;
        u_frd_sum_j = (U[i][j] + U[i][j + 1]) / 2;
        u_bkd_sum_i = (U[i][j] + U[i - 1][j]) / 2;
        u_bkd_sum_j = (U[i][j] + U[i][j - 1]) / 2;
        v_frd_sum_i = (V[i][j] + V[i + 1][j]) / 2;
        
        d2ux = (U[i - 1][j] - 2 * U[i][j] + U[i + 1][j]) / (dx * dx);
        d2uy = (U[i][j - 1] - 2 * U[i][j] + U[i][j + 1]) / (dy * dy);
        
        du2x = (u_frd_sum_i * u_frd_sum_i - u_bkd_sum_i * u_bkd_sum_i) / dx +
        alpha * (abs(u_frd_sum_i) * (U[i][j] - U[i + 1][j]) / 2 -
                 abs(u_bkd_sum_i) * (U[i - 1][j] - U[i][j]) / 2) / dx;
        
        duvy = (v_frd_sum_i * u_frd_sum_j - (V[i][j - 1] + V[i + 1][j - 1]) / 2 * u_bkd_sum_j) / dy +
        alpha * (abs(v_frd_sum_i) * (U[i][j] - U[i][j - 1]) / 2 -
                 abs(V[i][j - 1] + V[i + 1][j - 1]) * (U[i][j - 1] - U[i][j]) / 4) / dy;
        
        F[i][j] = U[i][j] + dt * ((d2ux + d2uy) / Re - du2x - duvy + GX);
    }
    
    i = imax;
    
    for (j = 1; j < jmax; j++) {
        u_frd_sum_j = (U[i][j] + U[i][j + 1]) / 2;
        
        v_frd_sum_i = (V[i][j] + V[i + 1][j]) / 2;
        v_frd_sum_j = (V[i][j] + V[i][j + 1]) / 2;
        
        v_bkd_sum_i = (V[i][j] + V[i - 1][j]) / 2;
        v_bkd_sum_j = (V[i][j] + V[i][j - 1]) / 2;
        
        d2vx = (V[i - 1][j] - 2 * V[i][j] + V[i + 1][j]) / (dx * dx);
        d2vy = (V[i][j - 1] - 2 * V[i][j] + V[i][j + 1]) / (dy * dy);
        
        dv2y = (v_frd_sum_j * v_frd_sum_j - v_bkd_sum_j * v_bkd_sum_j) / dy +
        alpha * (abs(v_frd_sum_j) * (V[i][j] - V[i][j + 1]) / 2 -
                 abs(v_bkd_sum_j) * (V[i][j - 1] - V[i][j]) / 2) / dy;
        
        duvx = (u_frd_sum_j * v_frd_sum_i - (U[i - 1][j] + U[i - 1][j + 1]) / 2 * v_bkd_sum_i) / dx +
        alpha * (abs(u_frd_sum_j) * (V[i][j] - V[i + 1][j]) / 2 -
                 abs(U[i - 1][j] + U[i - 1][j + 1]) * (V[i - 1][j] - V[i][j]) / 4) / dx;
        
        G[i][j] = V[i][j] + dt * ((d2vx + d2vy) / Re - duvx - dv2y + GY);
    }


}


/**
 * This operation computes the right hand side of the pressure poisson equation.
 * The right hand side is computed according to the formula
 *
 * @f$ rs = \frac{1}{\delta t} \left( \frac{F^{(n)}_{i,j}-F^{(n)}_{i-1,j}}{\delta x} + \frac{G^{(n)}_{i,j}-G^{(n)}_{i,j-1}}{\delta y} \right)  @f$
 *
 */
void calculate_rs(
    double dt, double dx, double dy, int imax, int jmax, double **F, double **G, double **RS) {
    
    int i;
    int j;
    
    for (i = 1; i < imax + 1; i++) {
        for (j = 1; j < jmax + 1; j++) {
            RS[i][j] = ((F[i][j] - F[i - 1][j]) / dx + (G[i][j] - G[i][j - 1]) / dy) / dt;
        }
    }
}


/**
 * Determines the maximal time step size. The time step size is restricted
 * accordin to the CFL theorem. So the final time step size formula is given
 * by
 *
 * @f$ {\delta t} := \tau \, \min\left( \frac{Re}{2}\left(\frac{1}{{\delta x}^2} + \frac{1}{{\delta y}^2}\right)^{-1},  \frac{{\delta x}}{|u_{max}|},\frac{{\delta y}}{|v_{max}|} \right) @f$
 *
 */
void calculate_dt(double Re, double tau, double *dt, double dx, double dy, 
                  int imax, int jmax, double **U, double **V) {

    // let u_max = foldr (max . abs) 0 U
    double u_max = 0;
    double v_max = 0;

    for (int i = 0; i <= imax+1; i++) {
        for (int j = 0; j <= jmax+1; j++) {
            u_max = fmax(u_max, fabs(U[i][j]));
            v_max = fmax(v_max, fabs(V[i][j]));
        }
    }

    double dx2 = dx * dx;
    double dy2 = dy * dy;
    
    *dt = Re / 2 * (dx2 * dy2) / (dx2 + dy2);
    *dt = fmin(*dt, dx / u_max);
    *dt = fmin(*dt, dy / v_max);
    *dt *= tau;
}


/**
 * Calculates the new velocity values according to the formula
 *
 * @f$ u_{i,j}^{(n+1)}  =  F_{i,j}^{(n)} - \frac{\delta t}{\delta x} (p_{i+1,j}^{(n+1)} - p_{i,j}^{(n+1)}) @f$
 * @f$ v_{i,j}^{(n+1)}  =  G_{i,j}^{(n)} - \frac{\delta t}{\delta y} (p_{i,j+1}^{(n+1)} - p_{i,j}^{(n+1)}) @f$
 *
 * As always the index range is
 *
 * @f$ i=1,\ldots,imax-1, \quad j=1,\ldots,jmax @f$
 * @f$ i=1,\ldots,imax, \quad j=1,\ldots,jmax-1 @f$
 *
 * @image html calculate_uv.jpg
 */
void calculate_uv(
    double dt, double dx, double dy, int imax, int jmax, 
    double **U, double **V, double **F, double **G, double **P) {
    
    int i;
    int j;
    
    for (i = 1; i < imax; i++) {
        for (j = 1; j < jmax; j++) {
            U[i][j] = F[i][j] - dt / dx * (P[i + 1][j] - P[i][j]);
            V[i][j] = G[i][j] - dt / dy * (P[i][j + 1] - P[i][j]);
        }
        
        j = jmax;
        U[i][j] = F[i][j] - dt / dx * (P[i + 1][j] - P[i][j]);
        V[i][j] = G[i][j] - dt / dy * (P[i][j + 1] - P[i][j]);
    }
    
    i = imax;
    for (j = 1; j < jmax; j++) {
        U[i][j] = F[i][j] - dt / dx * (P[i + 1][j] - P[i][j]);
        V[i][j] = G[i][j] - dt / dy * (P[i][j + 1] - P[i][j]);
    }
}



