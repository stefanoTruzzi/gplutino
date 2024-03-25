#include "plutino.hpp"


/* ********************************************************************* */
int flux_single(double *v, double *u, double *f, Indices &indices)
/*
 *
 *********************************************************************** */
{
	double ptot; 
  // HD  = P
  // MHD = P + Pmag

  f[RHO] = u[indices.MXn];            // MHD = HD  
  f[MX1] = v[indices.VXn]*u[MX1];   // MHD = HD -= Bn * V
  f[MX2] = v[indices.VXn]*u[MX2];   
  f[MX3] = v[indices.VXn]*u[MX3];

  f[indices.MXn] += v[PRS];
  ptot = v[PRS];

  #if PHYSICS == HD 
  f[ENG] = (u[ENG] + ptot) * v[indices.VXn];
  #endif

  #if PHYSICS == IDEALMHD
  double Bmag2; // square of B 
  double vB;    // B*speed
    
  Bmag2 = v[BX1]*v[BX1] + v[BX2]*v[BX2] + v[BX3]*v[BX3];
  vB    = v[VX1]*v[BX1] + v[VX2]*v[BX2] + v[VX3]*v[BX3];
  
  f[indices.MXn] += 0.5*Bmag2;
  f[MX1] -= v[indices.BXn]*v[BX1];
  f[MX2] -= v[indices.BXn]*v[BX2];
  f[MX3] -= v[indices.BXn]*v[BX3];

  f[indices.BXn] = 0.0;
  f[indices.BXt] = v[indices.VXn]*v[indices.BXt] - v[indices.BXn]*v[indices.VXt];
  f[indices.BXb] = v[indices.VXn]*v[indices.BXb] - v[indices.BXn]*v[indices.VXb];
  
  ptot  += 0.5*Bmag2;
  f[ENG] = (u[ENG] + ptot)*v[indices.VXn] - v[indices.BXn]*vB;
  return (0);
  #endif

  return 0;
}




int ustar_single(double *u, double *v, double *uk, double sk, double sstar, int dir)
{
  double coeff;
#if 0

  coeff = v[RHO]*(sk - v[VXn])/(sk - sstar);
  uk[RHO] = coeff * 1.0; 
  uk[MX1] = coeff * v[VX1];
  uk[MX2] = coeff * v[VX2];
  uk[MX3] = coeff * v[VX3];
  uk[MXn] = coeff * sstar;
  uk[ENG] = coeff * (u[ENG]/v[RHO]) 
          + (sstar - v[VXn]) * ( sstar + v[PRS]/(v[RHO]*(sk-v[VXn])));
#endif  
  return(0);
}
