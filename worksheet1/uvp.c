#include "helper.h"
#include "uvp.h"

void calculate_fg(
    double Re,
    double GX,
    double GY,
    double alpha,
    double dt,
    double dx,
    double dy,
    int imax,
    int jmax,
    double **U,
    double **V,
    double **F,
    double **G){

  double parU2x, parUVy, parUVx, par2Ux, par2Uy, par2Vx, par2Vy, parV2y;

  for (int i = 1; i <= imax; ++i)  {
    for (int j = 1; j <= jmax; ++j)    {
      parU2x = (pow((U[i][j]+U[i+1][j])/2,2)-pow((U[i][j]+U[i-1][j])/2,2))/dx 
        + alpha * (abs(U[i][j]+U[i+1][j])*(U[i][j]-U[i+1][j])  - abs(U[i-1][j]+U[i][j])*(U[i-1][j]-U[i][j]))/(4*dx);
      parV2y = (pow((V[i][j]+V[i][j+1])/2,2)-pow((V[i][j]+V[i][j-1])/2,2))/dy 
        + alpha * (abs(V[i][j]+V[i][j+1])*(V[i][j]-V[i][j+1])  - abs(V[i][j-1]+V[i][j])*(V[i][j-1]-V[i][j]))/(4*dy);

      parUVy = (      (V[i][j]+V[i+1][j])*(U[i][j]+U[i][j+1])  -    (V[i][j-1]+V[i+1][j-1])*(U[i][j-1]-U[i][j]))/(4*dy)
        + alpha * (abs(V[i][j]+V[i+1][j])*(U[i][j]-U[i][j+1])  - abs(V[i][j-1]+V[i+1][j-1])*(U[i][j-1]-U[i][j]))/(4*dy);        
      parUVx = (      (U[i][j]+U[i][j+1])*(V[i][j]+V[i+1][j])  -    (U[i-1][j]+U[i-1][j+1])*(V[i-1][j]+V[i][j]))/(4*dx)
        + alpha * (abs(U[i][j]+U[i][j+1])*(V[i][j]-V[i+1][j])  - abs(U[i-1][j]+U[i-1][j+1])*(V[i-1][j]-V[i][j]))/(4*dx);              

      par2Ux = (U[i+1][j]-2*U[i][j]+U[i-1][j])/(dx*dx);
      par2Uy = (U[i][j+1]-2*U[i][j]+U[i][j-1])/(dy*dy);

      par2Vx = (V[i+1][j]-2*V[i][j]+V[i-1][j])/(dx*dx);
      par2Vy = (V[i][j+1]-2*V[i][j]+V[i][j-1])/(dy*dy);

      if (i<imax){        
        F[i][j] = U[i][j] + dt * ((par2Ux+par2Uy)/Re - parU2x - parUVy + GX);  
      }
      if (j<jmax)      {
        G[i][j] = V[i][j] + dt * ((par2Vx+par2Vy)/Re - parUVx - parV2y + GY);
      }
      
    }
  }
}


void calculate_rs(
    double dt,
    double dx,
    double dy,
    int imax,
    int jmax,
    double **F,
    double **G,
    double **RS){

  for (int i = 1; i < imax; ++i)  {
    for (int j = 1; j < jmax; ++j)    {
      RS [i][j] = ((F[i][j]-F[i-1][j])/dx + (G[i][j]-G[i][j-1])/dy )/dt;
    }
  }  
}



void calculate_dt(
    double Re,
    double tau,
    double *dt,
    double dx,
    double dy,
    int imax,
    int jmax,
    double **U,
    double **V){

  double Umax=0, Vmax=0;
  for (int i = 0; i < imax; ++i)
  {
    for (int j = 0; j < jmax; ++j)
    {
      Umax = fmax(Umax,abs(U[i][j]));
      Vmax = fmax(Vmax,abs(V[i][j]));
    }
  }
  if (Umax == 0){
    if (Vmax == 0){
      *dt = tau * ((Re/2)/((1/dx)+(1/dy)));
    }
    else{
      *dt = tau * fmin ((Re/2)/((1/dx)+(1/dy)) , (dx / Umax));
    }
  }
  else{
    if (Vmax == 0){
      *dt = tau * fmin ((Re/2)/((1/dx)+(1/dy)) , (dx / Umax));
    }
    else{
      *dt = tau * fmin ((Re/2)/((1/dx)+(1/dy)) , fmin ((dx / Umax),(dy / Vmax)));
    }
  }
  
}



void calculate_uv(
    double dt,
    double dx,
    double dy,
    int imax,
    int jmax,
    double **U,
    double **V,
    double **F,
    double **G,
    double **P){

  for (int i = 0; i <= imax; ++i)  {
    for (int j = 0; j <= jmax; ++j)    {
      if (i<imax)      {
        U[i][j] = F[i][j] - dt * (P[i+1][j]-P[i][j]) / dx;
      }
      if (j<jmax)      {
        V[i][j] = G[i][j] - dt * (P[i][j+1]-P[i][j]) / dy;
      }
    }
  }  
}
