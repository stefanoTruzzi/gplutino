#include "plutino.hpp"


/* ********************************************************************* */
int flux_single(double *v, double *u, double *f)
/*
 *
 *********************************************************************** */
{
  double ptot; 
  // HD  = P
  // MHD = P + Pmag

  f[RHO] = u[MXn];            // MHD = HD  
  f[MX1] = v[VXn]*u[MX1];   // MHD = HD -= Bn * V
  f[MX2] = v[VXn]*u[MX2];   
  f[MX3] = v[VXn]*u[MX3];

  f[MXn] += v[PRS];
  ptot = v[PRS];

  #if PHYSICS == HD 
  f[ENG] = (u[ENG] + ptot) * v[VXn];
  #endif

  #if PHYSICS == IDEALMHD
  double Bmag2; // square of B 
  double vB;    // B*speed
    
  Bmag2 = v[BX1]*v[BX1] + v[BX2]*v[BX2] + v[BX3]*v[BX3];
  vB    = v[VX1]*v[BX1] + v[VX2]*v[BX2] + v[VX3]*v[BX3];
  
  f[MXn] += 0.5*Bmag2;
  f[MX1] -= v[BXn]*v[BX1];
  f[MX2] -= v[BXn]*v[BX2];
  f[MX3] -= v[BXn]*v[BX3];

  f[BXn] = 0.0;
  f[BXt] = v[VXn]*v[BXt] - v[BXn]*v[VXt];
  f[BXb] = v[VXn]*v[BXb] - v[BXn]*v[VXb];
  
  ptot  += 0.5*Bmag2;
  f[ENG] = (u[ENG] + ptot)*v[VXn] - v[BXn]*vB;
  return (0);
  #endif

  return 0;
}




int ustar_single(double *u, double *v, double *uk, double sk, double sstar, int dir)
{
  double coeff;

  coeff = v[RHO]*(sk - v[VXn])/(sk - sstar);
  uk[RHO] = coeff * 1.0; 
  uk[MX1] = coeff * v[VX1];
  uk[MX2] = coeff * v[VX2];
  uk[MX3] = coeff * v[VX3];
  uk[MXn] = coeff * sstar;
  uk[ENG] = coeff * (u[ENG]/v[RHO]) 
          + (sstar - v[VXn]) * ( sstar + v[PRS]/(v[RHO]*(sk-v[VXn])));
  
  return(0);
}