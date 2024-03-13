#include "plutino.hpp"

        


int prim_to_cons(double *v, double *u)
{
  //rho vx P (vettori V)
  double Ekin;

  u[RHO] = v[RHO];
  
  u[MX1] = v[RHO]*v[VX1];
  u[MX2] = v[RHO]*v[VX2];
  u[MX3] = v[RHO]*v[VX3];

  Ekin   = v[RHO]*(v[VX1]*v[VX1] + v[VX2]*v[VX2] + v[VX3]*v[VX3]);
  
  #if PHYSICS == IDEALMHD
    u[BX1] = v[BX1];
    u[BX2] = v[BX2];
    u[BX3] = v[BX3];
    Ekin  +=  v[BX1]*v[BX1] + v[BX2]*v[BX2] + v[BX3]*v[BX3];
  #endif

  Ekin  *= 0.5;
  u[ENG] = Ekin + v[PRS]/(GAMMA_EOS-1.0);
  
  return 0;
}

int cons_to_prim(double *u, double *v){
  double m2, b2, kinb2;
  double tau, rho;
  
  m2 = u[MX1]*u[MX1] + u[MX2]*u[MX2] + u[MX3]*u[MX3];
  
  v[RHO] = rho = u[RHO];
  tau = 1.0/u[RHO];
  v[VX1] = u[MX1]*tau;
  v[VX2] = u[MX2]*tau;
  v[VX3] = u[MX3]*tau;

  kinb2 = m2*tau;    

  #if PHYSICS == IDEALMHD
    b2 = u[BX1]*u[BX1] + u[BX2]*u[BX2] + u[BX3]*u[BX3];
    v[BX1] = u[BX1];
    v[BX2] = u[BX2];
    v[BX3] = u[BX3];
    kinb2 += b2;
  #endif
  
  kinb2 *= 0.5;
  v[PRS] = (GAMMA_EOS-1.0)*(u[ENG] - kinb2);  
  
  return 0;
}

