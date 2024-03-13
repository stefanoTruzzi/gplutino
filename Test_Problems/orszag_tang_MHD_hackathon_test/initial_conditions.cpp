#include "plutino.hpp"

void Initial_Condition(double x1, double x2, double *v)
/*
 *  http://plutocode.ph.unito.it/Doxygen/Test_Problems/_m_h_d_2_orszag___tang_2init_8c.html
 *
 *
 *********************************************************************** */
{  
  double x = x1;
  double y = x2;

  v[RHO] = 25.0/9.0; 
  v[PRS] = 5.0/3.0;

  v[VX1] = -sin(2.0*M_PI*y);
  v[VX2] =  sin(2.0*M_PI*x);
  v[VX3] =  0.0;

  #if PHYSICS == IDEALMHD 
    v[BX1] = -sin(2.0*M_PI*y);
    v[BX2] =  sin(4.0*M_PI*x);
    v[BX3] = 0.0;
  #endif
}
