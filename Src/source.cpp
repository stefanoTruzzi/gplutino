#include "plutino.hpp"


int powell_source(double *v, double *s, double divB)
{
  // Correction divb
  s[RHO] = 0.0;
  s[MX1] = -v[BX1] * divB; 
  s[MX2] = -v[BX2] * divB;
  s[MX3] = -v[BX3] * divB; 
  
  s[BX1] = -v[VX1] * divB;
  s[BX2] = -v[VX2] * divB;
  s[BX3] = -v[VX3] * divB;
  
  s[ENG] = -divB* (v[VX1] * v[BX1] + v[VX2] * v[BX2] + v[VX3] * v[BX3]);

  return(0);
}

