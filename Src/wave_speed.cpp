#include "plutino.hpp"


double MaxSignalSpeedlr(double *vl, double *vr)
{
  int i; 
  double cl, cr, lambda_max;
  double lambda_l, lambda_r;
  cl = sqrt(GAMMA_EOS*vl[PRS]/vl[RHO]); 
  cr = sqrt(GAMMA_EOS*vr[PRS]/vr[RHO]);  
      
  #if PHYSICS == HD
    lambda_l = fabs(vl[VXn]);
    lambda_r = fabs(vr[VXn]);
    lambda_max = MAX(lambda_l + cl, lambda_r + cr);
  #endif

  #if PHYSICS == IDEALMHD
    double b2l, b2r,val,var, alfcl, alfcr, normcr, normcl;
    b2l = vl[BX1]*vl[BX1] + vl[BX2]*vl[BX2] + vl[BX3]*vl[BX3];
    b2r = vr[BX1]*vr[BX1] + vr[BX2]*vr[BX2] + vr[BX3]*vr[BX3];
    val = (b2l/vl[RHO]);
    var = (b2r/vr[RHO]);
    alfcl  = cl*cl + val;
    alfcr  = cr*cr + var; 
    normcl = 4.0*cl*cl*vl[BXn]*vl[BXn]/(vl[RHO]);
    normcr = 4.0*cr*cr*vr[BXn]*vr[BXn]/(vr[RHO]);
    lambda_l = sqrt(0.5*(alfcl + sqrt((alfcl*alfcl) - normcl)));
    lambda_r = sqrt(0.5*(alfcr + sqrt((alfcr*alfcr) - normcr)));
    lambda_max=MAX(fabs(vl[VXn])+lambda_l,fabs(vr[VXn])+lambda_r);

  #endif
  return (lambda_max);
}